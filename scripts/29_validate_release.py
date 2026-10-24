#!/usr/bin/env python
"""Phase 8 / step 29 — local validation of dmg_atlas_v1.h5ad against the CELLxGENE schema 7.1.0 text.

Stands in for `cellxgene-schema validate` (not in the multiverse env; no installs). Rules are transcribed from
https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/7.1.0/schema.md (read 2026-10-18) and restricted to what
applies to this dataset: human, non-spatial, UMI scRNA/snRNA/multiome-RNA, raw counts in X (no normalised matrix, no raw), no
perturbations, not a pre-analysis collection.

Inputs : data/processed/dmg_atlas_v1.h5ad (h5py; X streamed in row blocks), results/tables/release_manifest.csv (expected md5),
         results/tables/release_ontology_map.csv (term labels hard-coded in step 27), GENCODE v48 GTF (data/raw/gencode/, from step 27),
         OLS4 REST API (EBI) for every ontology term used: existence, obsolescence, label, hierarchical ancestors, served version.
Outputs (never overwritten if present, unless --force):
  results/tables/release_validation.csv        one row per check: section, check, level (MUST / SHOULD / info), status (pass / fail /
                                               warn / info / skip), n_bad, detail
  results/tables/release_validation_terms.csv  every ontology term used: field, term_id, n_cells, OLS label, step-27 label, obsolete,
                                               required ancestor found, ontology version served by OLS4
  docs/release_validation.md                   generated summary (input md5, counts by status, failures / warnings, limits)
Run: nice -n 10 python scripts/29_validate_release.py [--force] [--offline]
  --offline skips the OLS4 checks (status 'skip'). Exit status 1 if any MUST check fails.

Checks
  structure : AnnData groups, X csr encoding, no layers / raw, non-empty obsp; no '__' names, unique column names, no deprecated or
              portal-annotated (curator MUST NOT) names in obs / var / uns.
  X         : float32, finite, non-zero values are positive integers, column indices in range and strictly increasing per row,
              every cell has >= 1 non-zero, no two cells with identical values for all features (blake2b hash per row).
  obs       : required fields present, categorical with str categories (is_primary_data bool); value vocabularies and formats
              (sex, tissue_type, suspension_type, 'unknown' rules, disease ' || ' ordering, donor_id not 'na'); suspension_type
              allowed for the assay (EFO:0030080 10x -> cell|nucleus, EFO:0007045 ATAC-seq -> nucleus); development stage label
              consistent with age_years (floor years / months; group stages only where age is missing or 0).
  var       : versionless Ensembl IDs, unique, all present as gene_id in GENCODE v48 primary assembly; feature_is_filtered bool
              and all False (no raw matrix); gene_symbol vs v48 gene_name differences reported (info; the portal re-derives names).
  obsm      : X_<suffix> names, numeric ndarray, n_obs rows, >= 2 columns, no +-inf, not all NaN (NaN rows reported as info).
  uns       : organism_ontology_term_id NCBITaxon:9606, title, batch_condition keys in obs, default_embedding in obsm, no empty
              values (recursive), X_approximate_distribution vocabulary if present.
  ontology  : each term exists in OLS4, is not obsolete, has the required ancestor (assay: EFO:0002772 | EFO:0010183; cell type:
              CL term not in CL:0000255 / CL:0000257 / CL:0000548; stage: HsapDv:0000001; disease: MONDO:0000001 | MONDO:0021178;
              tissue: UBERON:0001062), and its OLS label matches the label hard-coded in step 27 (case-insensitive).

Limits (why a pass here is not certification)
  OLS4 serves current ontology releases, not the versions pinned by schema 7.1.0 (CL 2026-03-26, EFO 3.90.0, MONDO 2026-05-05,
  UBERON 2026-04-01, HsapDv 2025-01-23); a term obsoleted or re-parented after the pin is judged against the newer release.
  Ancestry uses OLS hierarchicalAncestors (is_a + part_of); the official validator may use a narrower relation set.
  "Most accurate term" requirements cannot be checked by software.
"""
import sys, re, gzip, time, hashlib, argparse, urllib.request, urllib.parse, urllib.error, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import dmg_atlas.utils as U  # noqa: E402  (thread caps first)
import h5py, numpy as np, pandas as pd  # noqa: E402
from anndata.io import read_elem  # noqa: E402
from dmg_atlas import ROOT, RAW, PROCESSED, RESULTS, LOGS  # noqa: E402

log = U.StepLogger(LOGS / "phase8_release.log", context="29_validate_release")
H5 = PROCESSED / "dmg_atlas_v1.h5ad"
TAB = RESULTS / "tables"
OUT_CSV, OUT_TERMS, OUT_MD = TAB / "release_validation.csv", TAB / "release_validation_terms.csv", ROOT / "docs/release_validation.md"
GTF = RAW / "gencode" / "gencode.v48.primary_assembly.annotation.gtf.gz"
OLS = "https://www.ebi.ac.uk/ols4/api"
SCHEMA_URL = "https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/7.1.0/schema.md"
ROW_BLOCK = 20_000

REQUIRED_OBS = ["assay_ontology_term_id", "cell_type_ontology_term_id", "development_stage_ontology_term_id", "disease_ontology_term_id",
                "donor_id", "self_reported_ethnicity_ontology_term_id", "sex_ontology_term_id", "suspension_type", "tissue_ontology_term_id",
                "tissue_type"]
OBS_PORTAL = {"assay", "cell_type", "development_stage", "disease", "experimental_condition", "observation_joinid", "perturbation_types",
              "self_reported_ethnicity", "sex", "tissue"}
OBS_DEPRECATED = {"organism", "organism_ontology_term_id", "ethnicity", "ethnicity_ontology_term_id"}
OBS_PERTURB = {"experimental_condition_ontology_term_id", "genetic_perturbation_id", "genetic_perturbation_strategy"}
VAR_PORTAL = {"feature_biotype", "feature_length", "feature_name", "feature_reference", "feature_type"}
UNS_PORTAL = {"citation", "is_pre_analysis", "organism", "schema_reference", "schema_version"}
UNS_DEPRECATED = {"X_normalization", "default_field", "layer_descriptions", "tags", "version", "contributors", "preprint_doi",
                  "project_description", "project_links", "project_name", "publication_doi"}
CL_FORBIDDEN = {"CL:0000255", "CL:0000257", "CL:0000548"}
ONTO = {  # field -> (id regex, required ancestors (any; the root itself also accepted))
    "assay_ontology_term_id": (r"EFO:\d{7}", {"EFO:0002772", "EFO:0010183"}),
    "cell_type_ontology_term_id": (r"CL:\d{7}", {"CL:0000000"}),
    "development_stage_ontology_term_id": (r"HsapDv:\d{7}", {"HsapDv:0000001"}),
    "disease_ontology_term_id": (r"MONDO:\d{7}", {"MONDO:0000001", "MONDO:0021178"}),
    "self_reported_ethnicity_ontology_term_id": (r"(HANCESTRO|AfPO):\d{4,7}", {"HANCESTRO:0601", "HANCESTRO:0602"}),
    "tissue_ontology_term_id": (r"UBERON:\d{7}", {"UBERON:0001062"}),
    "sex_ontology_term_id": (r"PATO:\d{7}", set()),
}
UNKNOWN_OK = {"cell_type_ontology_term_id", "development_stage_ontology_term_id", "self_reported_ethnicity_ontology_term_id",
              "sex_ontology_term_id"}
SEX_OK = {"PATO:0000383", "PATO:0000384", "PATO:0001340", "unknown"}
SUSPENSION_RULES = {"EFO:0030080": {"cell", "nucleus"}, "EFO:0007045": {"nucleus"}}  # assay (and descendants) -> allowed

rows = []


def chk(section, name, ok, level="MUST", n_bad=0, detail=""):
    status = "info" if level == "info" else "pass" if ok else ("fail" if level == "MUST" else "warn")
    rows.append(dict(section=section, check=name, level=level, status=status, n_bad=int(n_bad), detail=str(detail)))
    if status in ("fail", "warn"):
        log.write(f"{status.upper()} [{section}] {name}: n_bad {n_bad} {detail}")


def skip(section, name, why):
    rows.append(dict(section=section, check=name, level="MUST", status="skip", n_bad=0, detail=why))


def md5(p: Path) -> str:
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()


def bad_names(names) -> list:
    return [n for n in names if str(n).startswith("__")]


def get_json(url: str):
    for attempt in range(4):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"Accept": "application/json"}), timeout=60) as r:
                return json.load(r)
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return {}
            err = e
        except OSError as e:
            err = e
        time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"OLS4 unreachable: {url} ({err})")


def ols_term(tid: str) -> dict:
    ont = tid.split(":")[0].lower()
    d = get_json(f"{OLS}/ontologies/{ont}/terms?obo_id={urllib.parse.quote(tid)}")
    ts = [t for t in d.get("_embedded", {}).get("terms", []) if t.get("obo_id") == tid]
    if not ts:
        return dict(found=False, label="", obsolete=False, ancestors=set())
    t = next((t for t in ts if t.get("is_defining_ontology")), ts[0])
    anc, url = set(), t["_links"]["hierarchicalAncestors"]["href"].replace("http://", "https://") + "?size=1000"
    while url:
        a = get_json(url)
        anc |= {x.get("obo_id") for x in a.get("_embedded", {}).get("terms", [])}
        url = a.get("_links", {}).get("next", {}).get("href", "").replace("http://", "https://")
    return dict(found=True, label=t.get("label", ""), obsolete=bool(t.get("is_obsolete")), ancestors=anc)


def ols_version(ont: str) -> str:
    c = get_json(f"{OLS}/ontologies/{ont}").get("config", {})
    return str(c.get("version") or c.get("versionIri") or "unknown")


def is_str_cat(s: pd.Series) -> bool:
    return isinstance(s.dtype, pd.CategoricalDtype) and all(isinstance(c, str) for c in s.cat.categories)


def check_structure(f):
    chk("structure", "required groups X, obs, var, obsm, uns present", {"X", "obs", "var", "obsm", "uns"} <= set(f.keys()),
        detail=", ".join(f.keys()))
    enc = f["X"].attrs.get("encoding-type", "dense" if isinstance(f["X"], h5py.Dataset) else "?")
    chk("structure", "X encoded as csr_matrix (>= 50% zeros)", enc == "csr_matrix", detail=enc)
    chk("structure", "no raw (raw counts already in X)", "raw" not in f, detail="raw present" if "raw" in f else "")
    nl = len(f["layers"]) if "layers" in f else 0
    chk("structure", "no extra layers (all layers must share cells/genes; none expected)", nl == 0, level="SHOULD", n_bad=nl)
    for g in ("obsp", "varm", "varp"):
        empty = [k for k in f.get(g, {}) if np.prod(f[g][k].attrs.get("shape", getattr(f[g][k], "shape", (1,)))) == 0]
        chk("structure", f"{g}: no zero-size entries", not empty, n_bad=len(empty), detail=", ".join(empty))


def check_obs(obs: pd.DataFrame, col_order) -> pd.DataFrame:
    sec = "obs"
    chk(sec, "index unique", obs.index.is_unique, n_bad=int(obs.index.duplicated().sum()))
    chk(sec, "index is str", all(isinstance(i, str) for i in obs.index[:1000]))
    dup = pd.Index(col_order)[pd.Index(col_order).duplicated()].tolist()
    chk(sec, "column names unique", not dup, n_bad=len(dup), detail=", ".join(dup))
    b = bad_names(obs.columns); chk(sec, "no names starting with '__'", not b, n_bad=len(b), detail=", ".join(b))
    for label, names in [("deprecated reserved names absent", OBS_DEPRECATED), ("portal-annotated names absent (curator MUST NOT)", OBS_PORTAL)]:
        hit = sorted(names & set(obs.columns)); chk(sec, label, not hit, n_bad=len(hit), detail=", ".join(hit))
    hit = sorted(OBS_PERTURB & set(obs.columns))
    chk(sec, "7.1.0 perturbation fields absent (dataset has no perturbations)", not hit, level="SHOULD", n_bad=len(hit), detail=", ".join(hit))
    miss = [c for c in REQUIRED_OBS + ["is_primary_data"] if c not in obs]
    chk(sec, "required fields present", not miss, n_bad=len(miss), detail=", ".join(miss))
    for c in [c for c in REQUIRED_OBS if c in obs]:
        chk(sec, f"{c}: categorical with str categories", is_str_cat(obs[c]), detail=str(obs[c].dtype))
        chk(sec, f"{c}: no missing values", obs[c].notna().all(), n_bad=int(obs[c].isna().sum()))
        unused = set(obs[c].cat.categories) - set(obs[c].unique()) if isinstance(obs[c].dtype, pd.CategoricalDtype) else set()
        chk(sec, f"{c}: no unused categories", not unused, level="SHOULD", n_bad=len(unused), detail=", ".join(sorted(unused)))
    if "is_primary_data" in obs:
        chk(sec, "is_primary_data: bool", obs["is_primary_data"].dtype == bool, detail=str(obs["is_primary_data"].dtype))
        chk(sec, "is_primary_data: value counts", True, level="info", detail=obs["is_primary_data"].value_counts().to_dict())
    s = lambda c: obs[c].astype(str)  # noqa: E731
    for c, (rx, _) in ONTO.items():
        v = s(c)
        ok = v.str.fullmatch(rx) | (v.isin(["unknown"]) if c in UNKNOWN_OK else False)
        if c == "disease_ontology_term_id":
            parts = v.str.split(" || ", regex=False)
            ok = (v == "PATO:0000461") | parts.map(lambda p: all(re.fullmatch(rx, t) for t in p) and p == sorted(set(p)))
        chk(sec, f"{c}: id format / allowed literal", ok.all(), n_bad=int((~ok).sum()), detail=", ".join(sorted(set(v[~ok]))[:10]))
    v = s("sex_ontology_term_id"); chk(sec, "sex_ontology_term_id in female/male/hermaphrodite/unknown", v.isin(SEX_OK).all(), n_bad=int((~v.isin(SEX_OK)).sum()))
    v = s("tissue_type"); ok = v.isin(["tissue", "primary cell culture", "cell line", "organoid"])
    chk(sec, "tissue_type vocabulary", ok.all(), n_bad=int((~ok).sum()), detail=", ".join(sorted(set(v))))
    cl = s("tissue_type") == "cell line"
    for c in ["cell_type_ontology_term_id", "development_stage_ontology_term_id", "donor_id", "sex_ontology_term_id",
              "self_reported_ethnicity_ontology_term_id"]:
        na = s(c) == "na"; bad = (na & ~cl) | (~na & cl & (c != "cell_type_ontology_term_id"))  # cell type: 'na' MAY for cell lines
        chk(sec, f"{c}: 'na' only (and then required) for cell lines", not bad.any(), n_bad=int(bad.sum()))
    v = s("suspension_type"); ok = v.isin(["cell", "nucleus", "na"])
    chk(sec, "suspension_type vocabulary", ok.all(), n_bad=int((~ok).sum()), detail=", ".join(sorted(set(v))))
    d = s("donor_id"); chk(sec, "donor_id non-empty", (d.str.strip() != "").all(), n_bad=int((d.str.strip() == "").sum()), detail=f"{d.nunique()} donors")
    ct = s("cell_type_ontology_term_id"); bad = ct.isin(CL_FORBIDDEN)
    chk(sec, "cell_type_ontology_term_id not a forbidden CL root", not bad.any(), n_bad=int(bad.sum()))
    chk(sec, "cell_type_ontology_term_id 'unknown' cells", True, level="info", n_bad=int((ct == "unknown").sum()),
        detail="allowed when no appropriate term exists")
    return obs


def check_stage_vs_age(obs: pd.DataFrame, labels: dict):
    if "age_years" not in obs:
        return skip("obs", "development stage consistent with age_years", "no age_years column")
    t = pd.DataFrame({"age": obs["age_years"].astype(float).values, "term": obs["development_stage_ontology_term_id"].astype(str).values})
    bad = 0; ex = []
    for (age, term), n in t.value_counts(dropna=False).items():
        lab = labels.get(term, ""); m = re.fullmatch(r"(\d+)-(year|month)-old stage", lab)
        if m:
            exp = int(age) if m.group(2) == "year" else int(age * 12) if np.isfinite(age) else -1
            ok = np.isfinite(age) and age > 0 and int(m.group(1)) == exp
        else:
            ok = not (np.isfinite(age) and age > 0) or term == "unknown"
        if not ok:
            bad += n; ex.append(f"age {age:g} -> {term} '{lab}' ({n})")
    chk("obs", "development stage consistent with age_years (floor years / months)", bad == 0, n_bad=bad, detail="; ".join(ex[:8]))


def check_var(var: pd.DataFrame, col_order):
    sec = "var"
    idx = var.index.astype(str)
    chk(sec, "index unique", idx.is_unique, n_bad=int(idx.duplicated().sum()))
    ok = idx.str.fullmatch(r"ENSG\d{11}")
    chk(sec, "index = versionless human Ensembl gene IDs", ok.all(), n_bad=int((~ok).sum()), detail=", ".join(idx[~ok][:5]))
    dup = pd.Index(col_order)[pd.Index(col_order).duplicated()].tolist()
    chk(sec, "column names unique", not dup, n_bad=len(dup))
    b = bad_names(var.columns); chk(sec, "no names starting with '__'", not b, n_bad=len(b))
    hit = sorted(VAR_PORTAL & set(var.columns)); chk(sec, "portal-annotated names absent (curator MUST NOT)", not hit, n_bad=len(hit), detail=", ".join(hit))
    fif = var.get("feature_is_filtered")
    chk(sec, "feature_is_filtered present, bool", fif is not None and fif.dtype == bool, detail=str(getattr(fif, "dtype", "missing")))
    if fif is not None:
        chk(sec, "feature_is_filtered all False (no raw matrix)", not fif.any(), n_bad=int(fif.sum()))
    if not GTF.exists():
        return skip(sec, "gene IDs in GENCODE v48", f"{GTF} missing (run step 27)")
    with log.step("GENCODE v48 gene IDs"):
        g = {}
        with gzip.open(GTF, "rt") as fh:
            for line in fh:
                if line[0] != "#":
                    p = line.split("\t", 9)
                    if p[2] == "gene":
                        gid = p[8].split('gene_id "', 1)[1].split('"', 1)[0].split(".", 1)[0]
                        g[gid] = p[8].split('gene_name "', 1)[1].split('"', 1)[0] if 'gene_name "' in p[8] else gid
        g = pd.Series(g)
    miss = ~idx.isin(g.index)
    chk(sec, "every gene ID is a GENCODE v48 gene_id", not miss.any(), n_bad=int(miss.sum()), detail=", ".join(idx[miss][:10]))
    if "gene_symbol" in var:
        sym = var["gene_symbol"].astype(str).values; ref = g.reindex(idx).fillna("").values
        diff = sym != ref
        chk(sec, "gene_symbol differs from v48 gene_name (portal uses gene_name)", True, level="info", n_bad=int(diff.sum()),
            detail=", ".join(f"{a}->{b}" for a, b in list(zip(sym[diff], ref[diff]))[:10]))


def check_obsm(f, n_obs):
    sec = "obsm"
    keys = list(f["obsm"].keys())
    chk(sec, "at least one X_ embedding (non-spatial assay)", any(k.startswith("X_") for k in keys), detail=", ".join(keys))
    for k in keys:
        okname = re.fullmatch(r"X_[a-zA-Z][a-zA-Z0-9_.-]*", k) is not None and k != "X_spatial"
        chk(sec, f"{k}: name X_<suffix> pattern", okname or k == "spatial", detail=k)
        a = f["obsm"][k]
        if not isinstance(a, h5py.Dataset):
            chk(sec, f"{k}: numpy.ndarray", False, detail=a.attrs.get("encoding-type", "group")); continue
        a = a[:]
        chk(sec, f"{k}: numeric dtype, n_obs rows, >= 2 columns", a.dtype.kind in "fiu" and a.shape[0] == n_obs and a.ndim == 2 and a.shape[1] >= 2,
            detail=f"{a.dtype} {a.shape}")
        if a.dtype.kind == "f":
            chk(sec, f"{k}: no +-inf", not np.isinf(a).any(), n_bad=int(np.isinf(a).any(axis=1).sum()))
            nanrow = np.isnan(a).any(axis=1)
            chk(sec, f"{k}: not all NaN", not np.isnan(a).all())
            if nanrow.any():
                # audit P7-21: the official CELLxGENE validator WARNS on NaN embedding rows; record it at the same level so this
                # validator cannot look cleaner than the one that gates submission. chk() maps ok=True to 'pass' whatever the
                # level, so the condition is stated negatively ("no rows with NaN") and fails here at SHOULD level.
                chk(sec, f"{k}: no rows with NaN", False, level="SHOULD", n_bad=int(nanrow.sum()),
                    detail="the official cellxgene-schema validator warns on these rows (not all NaN, so not a MUST failure)")


def check_uns(uns: dict, obs_cols, obsm_keys):
    sec = "uns"
    for label, names in [("deprecated reserved keys absent", UNS_DEPRECATED), ("portal-annotated keys absent (curator MUST NOT)", UNS_PORTAL)]:
        hit = sorted(names & set(uns)); chk(sec, label, not hit, n_bad=len(hit), detail=", ".join(hit))
    b = bad_names(uns); chk(sec, "no keys starting with '__'", not b, n_bad=len(b))
    org = uns.get("organism_ontology_term_id")
    chk(sec, "organism_ontology_term_id == NCBITaxon:9606", str(org) == "NCBITaxon:9606", detail=str(org))
    t = uns.get("title"); chk(sec, "title is a non-empty str", isinstance(t, str) and t.strip() != "", detail=str(t)[:80])
    bc = uns.get("batch_condition")
    if bc is not None:
        bc = [str(x) for x in np.atleast_1d(bc)]; miss = [x for x in bc if x not in obs_cols]
        chk(sec, "batch_condition keys are obs columns", not miss, n_bad=len(miss), detail=", ".join(bc))
    de = uns.get("default_embedding")
    if de is not None:
        chk(sec, "default_embedding is an obsm key", str(de) in obsm_keys, detail=str(de))
    xd = uns.get("X_approximate_distribution")
    chk(sec, "X_approximate_distribution in {count, normal} if present", xd is None or str(xd) in ("count", "normal"), detail=str(xd))
    empty = []

    def walk(d, path):
        for k, v in d.items():
            p = f"{path}/{k}"
            if isinstance(v, dict):
                (empty.append(p) if len(v) == 0 else walk(v, p))
            elif v is None or isinstance(v, (bool, np.bool_)):
                continue
            elif (len(v) if isinstance(v, (str, list, pd.DataFrame)) else np.size(v)) == 0:
                empty.append(p)
    walk(uns, "uns")
    chk(sec, "no empty values (True / False / None allowed)", not empty, n_bad=len(empty), detail=", ".join(empty[:10]))
    for k in [k for k in uns if k.endswith("_colors")]:
        col = k[:-7]
        chk(sec, f"{k}: refers to an obs column", col in obs_cols, level="MUST", detail=col)


def check_ontology(obs: pd.DataFrame, offline: bool) -> tuple:
    sec = "ontology"
    m27 = pd.read_csv(TAB / "release_ontology_map.csv").drop_duplicates(["field", "term_id"]).set_index(["field", "term_id"])["term_label"]
    terms = []
    for c in ONTO:
        vc = obs[c].astype(str).value_counts()
        for v, n in vc.items():
            for t in ([v] if c != "disease_ontology_term_id" else v.split(" || ")):
                if t not in ("unknown", "na", "PATO:0000461"):
                    terms.append(dict(field=c, term_id=t, n_cells=int(n), label_step27=m27.get((c, t), "")))
    terms.append(dict(field="uns.organism_ontology_term_id", term_id="NCBITaxon:9606", n_cells=int(len(obs)), label_step27="Homo sapiens"))
    terms = pd.DataFrame(terms).groupby(["field", "term_id", "label_step27"], as_index=False)["n_cells"].sum()
    if offline:
        skip(sec, "ontology terms exist / not obsolete / required ancestor / label", "--offline")
        return terms, {}
    with log.step(f"OLS4 lookups ({terms.term_id.nunique()} terms)"):
        info = {t: ols_term(t) for t in terms.term_id.unique()}
        versions = {o: ols_version(o) for o in sorted({t.split(':')[0].lower() for t in terms.term_id})}
    terms["ols_label"] = terms.term_id.map(lambda t: info[t]["label"])
    terms["found"] = terms.term_id.map(lambda t: info[t]["found"])
    terms["obsolete"] = terms.term_id.map(lambda t: info[t]["obsolete"])
    roots = {**{k: v[1] for k, v in ONTO.items()}, "uns.organism_ontology_term_id": set()}
    terms["ancestor_ok"] = [not roots[f] or t in roots[f] or bool(roots[f] & info[t]["ancestors"]) for f, t in zip(terms.field, terms.term_id)]
    terms["label_match"] = [a.strip().lower() == b.strip().lower() or not a for a, b in zip(terms.label_step27, terms.ols_label)]
    terms["ols_version"] = terms.term_id.map(lambda t: versions[t.split(":")[0].lower()])
    for name, col in [("term exists in OLS4", "found"), ("term not obsolete", "obsolete"), ("term has required ancestor", "ancestor_ok")]:
        bad = ~terms[col] if col != "obsolete" else terms[col]
        bad &= terms.found if col != "found" else True
        chk(sec, name, not bad.any(), n_bad=int(terms.n_cells[bad].sum()),
            detail="; ".join(f"{f}={t}" for f, t in zip(terms.field[bad], terms.term_id[bad])))
    bad = ~terms.label_match & terms.found
    chk(sec, "OLS label matches label hard-coded in step 27", not bad.any(), level="SHOULD", n_bad=int(terms.n_cells[bad].sum()),
        detail="; ".join(f"{t}: step27 '{a}' vs OLS '{b}'" for t, a, b in zip(terms.term_id[bad], terms.label_step27[bad], terms.ols_label[bad])))
    s = obs["suspension_type"].astype(str); a = obs["assay_ontology_term_id"].astype(str)
    bad_n, ex = 0, []
    for (assay, susp), n in pd.crosstab(a, s).stack().items():
        if n == 0:
            continue
        anc = info[assay]["ancestors"] | {assay}
        allowed = [v for k, v in SUSPENSION_RULES.items() if k in anc]
        if not allowed:
            chk(sec, f"suspension_type rule for {assay}", False, level="SHOULD", n_bad=int(n), detail="assay not in the schema table")
        elif not all(susp in v for v in allowed):
            bad_n += int(n); ex.append(f"{assay}/{susp} ({n})")
    chk(sec, "suspension_type allowed for assay", bad_n == 0, n_bad=bad_n, detail="; ".join(ex))
    return terms, versions


def check_X(f, n_obs, n_var):
    sec = "X"
    x = f["X"]; shape = tuple(int(i) for i in x.attrs["shape"])
    chk(sec, "shape = (n_obs, n_var)", shape == (n_obs, n_var), detail=str(shape))
    chk(sec, "data dtype float32", x["data"].dtype == np.float32, detail=str(x["data"].dtype))
    with log.step("stream X (values, indices, empty and duplicate rows)"):
        indptr = x["indptr"][:]
        chk(sec, "indptr monotone, starts at 0, ends at nnz", indptr[0] == 0 and (np.diff(indptr) >= 0).all() and indptr[-1] == x["data"].shape[0])
        nonfinite = nonpos = nonint = badix = unsorted = 0
        empty = int((np.diff(indptr) == 0).sum()); hashes = {}
        for r0 in range(0, n_obs, ROW_BLOCK):
            r1 = min(r0 + ROW_BLOCK, n_obs); a, b = indptr[r0], indptr[r1]
            d = x["data"][a:b]; ix = x["indices"][a:b]
            fin = np.isfinite(d); nonfinite += int((~fin).sum())
            nonpos += int((d[fin] <= 0).sum()); nonint += int((d[fin] != np.round(d[fin])).sum())
            badix += int(((ix < 0) | (ix >= n_var)).sum())
            row = np.repeat(np.arange(r1 - r0), np.diff(indptr[r0:r1 + 1]))
            unsorted += int(((np.diff(ix) <= 0) & (row[1:] == row[:-1])).sum())
            p = indptr[r0:r1 + 1] - a
            for i in range(r1 - r0):
                h = hashlib.blake2b(ix[p[i]:p[i + 1]].tobytes() + d[p[i]:p[i + 1]].tobytes(), digest_size=16).digest()
                hashes.setdefault(h, []).append(r0 + i)
    nnz = int(indptr[-1])
    chk(sec, "density < 50% (sparse encoding required)", True, level="info", detail=f"nnz {nnz:,}, density {nnz / (n_obs * n_var):.4f}")
    chk(sec, "values finite", nonfinite == 0, n_bad=nonfinite)
    chk(sec, "non-zero values positive", nonpos == 0, n_bad=nonpos)
    chk(sec, "non-zero values integers (de-duplicated UMI counts)", nonint == 0, n_bad=nonint)
    chk(sec, "column indices in range", badix == 0, n_bad=badix)
    chk(sec, "column indices strictly increasing within rows (canonical CSR)", unsorted == 0, n_bad=unsorted)
    chk(sec, "every cell has >= 1 non-zero value", empty == 0, n_bad=empty)
    dups = [v for v in hashes.values() if len(v) > 1]
    return dups


def write_report(h5_md5, man_md5, terms, versions, dups, obs_index):
    r = pd.DataFrame(rows)
    r.to_csv(OUT_CSV, index=False)
    terms.to_csv(OUT_TERMS, index=False)
    n = r.status.value_counts()
    lines = ["# Release v1 — local schema validation", "",
             f"Generated by `scripts/29_validate_release.py` on {time.strftime('%Y-%m-%d')}. Do not edit by hand.", "",
             f"Target: CELLxGENE schema 7.1.0 ({SCHEMA_URL}), transcribed into local checks because `cellxgene-schema` is not installed.", "",
             f"- Input: `data/processed/dmg_atlas_v1.h5ad` md5 `{h5_md5}` "
             + ("(matches release_manifest.csv)." if h5_md5 == man_md5 else f"(**differs** from release_manifest.csv `{man_md5}`)."),
             f"- Checks: {len(r)} — " + ", ".join(f"{k} {v}" for k, v in n.items()) + ".",
             f"- Ontology terms: {terms.term_id.nunique()} (`results/tables/release_validation_terms.csv`); OLS4 versions: "
             + (", ".join(f"{k} {v}" for k, v in versions.items()) if versions else "not queried (--offline)") + ".", ""]
    for st, head in [("fail", "Failures (MUST)"), ("warn", "Warnings (SHOULD)"), ("skip", "Skipped")]:
        sub = r[r.status == st]
        lines += [f"## {head}", ""] + (["None.", ""] if sub.empty else
                                        ["| section | check | n_bad | detail |", "|---|---|---|---|"]
                                        + [f"| {x.section} | {x.check} | {x.n_bad:,} | {x.detail[:300]} |" for x in sub.itertuples()] + [""])
    if dups:
        lines += ["Duplicate-row examples: " + "; ".join(" = ".join(obs_index[i] for i in g[:3]) for g in dups[:5]), ""]
    info = r[r.status == "info"]
    lines += ["## Notes (info)", "", "| section | check | n | detail |", "|---|---|---|---|"]
    lines += [f"| {x.section} | {x.check} | {x.n_bad:,} | {x.detail[:300]} |" for x in info.itertuples()]
    lines += ["", "## Limits", "",
              "- OLS4 serves current ontology releases, not the versions pinned by schema 7.1.0 (CL 2026-03-26, EFO 3.90.0, MONDO 2026-05-05, "
              "UBERON 2026-04-01, HsapDv 2025-01-23); obsolescence and ancestry are judged against the served release.",
              "- Ancestry uses OLS hierarchicalAncestors (is_a + part_of); the official validator may use a narrower relation set.",
              "- 'Most accurate term' requirements are not machine-checkable. A pass here is strong evidence, not certification: run "
              "`cellxgene-schema validate` (official validator) before submission.", ""]
    OUT_MD.write_text("\n".join(lines) + "\n")
    log.write(f"wrote {OUT_CSV.name}, {OUT_TERMS.name}, {OUT_MD.relative_to(ROOT)}; status counts {n.to_dict()}")
    return r


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--offline", action="store_true")
    args = ap.parse_args()
    if OUT_CSV.exists() and not args.force:
        log.write(f"{OUT_CSV.name} exists, skipping (use --force)"); sys.exit(0)
    with log.step("input md5"):
        h5_md5 = md5(H5)
        man = pd.read_csv(TAB / "release_manifest.csv").set_index("path")["md5"]
        man_md5 = man.get("data/processed/dmg_atlas_v1.h5ad", "")
        chk("structure", "h5ad md5 matches release_manifest.csv", h5_md5 == man_md5, detail=f"{h5_md5} vs {man_md5}")
        log.write(f"{H5.name} md5 {h5_md5}")
    with h5py.File(H5, "r") as f:
        with log.step("structure, obs, var, obsm, uns"):
            check_structure(f)
            obs = read_elem(f["obs"]); var = read_elem(f["var"]); uns = read_elem(f["uns"])
            check_obs(obs, list(f["obs"].attrs["column-order"]))
            check_var(var, list(f["var"].attrs["column-order"]))
            check_obsm(f, len(obs))
            check_uns(uns, set(obs.columns), set(f["obsm"].keys()))
        terms, versions = check_ontology(obs, args.offline)
        check_stage_vs_age(obs, dict(zip(terms.term_id, terms.get("ols_label", terms.label_step27))))
        dups = check_X(f, len(obs), len(var))
    chk("X", "no two cells identical across all features", not dups, n_bad=sum(len(g) for g in dups),
        detail=f"{len(dups)} groups" if dups else "")
    r = write_report(h5_md5, man_md5, terms, versions, dups, obs.index)
    fails = int((r.status == "fail").sum())
    log.write(f"done | {fails} MUST failures | peak RSS {U.peak_rss_gb():.2f} GB")
    sys.exit(1 if fails else 0)
