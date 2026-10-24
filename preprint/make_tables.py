"""Generate the preprint's LaTeX tables from results/tables/*.csv (run from the repo root)."""
import pandas as pd, pathlib, re
R = pathlib.Path(__file__).resolve().parents[1]
T = lambda n: pd.read_csv(R / "results" / "tables" / n)
OUT = R / "preprint" / "tables"; OUT.mkdir(exist_ok=True)
esc = lambda s: str(s).replace("_", r"\_").replace("&", r"\&").replace("%", r"\%")

def write(name, body):
    (OUT / name).write_text(body + "\n"); print("wrote", name)

# Deposited data: ONLY the deposits that actually enter the released results. The selection CSV's
# recommended_role is a selection-time proposal, not the inclusion record (it still carries Xie 2026, the Jessa
# mouse reference and ScPCA, all of which were dropped), so the list is taken from the pipeline outputs:
# atlas_v1_summary.csv (reference + projected queries) and results/qc/spatial_sections_used.csv (spatial).
USED = {  # pipeline study key -> (display name, accession, source)
    "okada":          ("Okada et al., post-mortem DMG (snRNA)", "EGAS50000001447 (raw, controlled); processed matrix from the study's UCSC Cell Browser instance", "EGA / UCSC Cell Browser"),
    "jessa_scrna":    ("Jessa et al., DMG (sc/snRNA)", "GSE210567", "GEO"),
    "jessa_multiome": ("Jessa et al., DMG (multiome RNA)", "GSE210566", "GEO"),
    "collot":         ("Collot et al., DMG (snRNA)", "Zenodo: 10.5281/zenodo.17922773 (scRNA.tar.gz)", "Zenodo"),
    "zang":           ("Zang et al., DMG (scRNA)", "GSE245674 (GSM7846769)", "GEO"),
    "desisto":        ("DeSisto et al., paediatric HGG (scRNA)", "GSE231859", "GEO"),
    "filbin":         ("Filbin et al., H3K27M glioma (Smart-seq2, TPM)", "GSE102130", "GEO"),
    "liu":            ("Liu et al., H3-K27M DMG (Smart-seq2, TPM)", "GSE184357", "GEO"),
    "labelle":        ("DMG immune cohort (SS2 + 10x, TPM)", "GSE227983", "GEO"),
    "braun2023":      ("Braun et al., first-trimester human brain", "CZ CELLxGENE Discover (collection 4d8fed08-2d6d-4692-b5ea-464f1d072077)", "CZ CELLxGENE"),
    "siletti2023":    ("Siletti et al., adult human brain (pons dissections)", "CZ CELLxGENE Discover (collection 283d65eb-dd53-496d-adb7-7570c7caa443)", "CZ CELLxGENE"),
    "huchede":        ("Huchede et al., DMG Visium", "GSE268577", "GEO"),
    "ren":            ("Ren et al., DMG and glioblastoma Visium", "GSE194329", "GEO"),
    "collot_st":      ("Collot et al., DMG Visium", "Zenodo: 10.5281/zenodo.17922773 (ST.tar.gz)", "Zenodo"),
    "brainstem":      ("Whole-brainstem DIPG autopsy Visium", "GSE280990", "GEO"),
}

# ---- Table 1: reference cohort ----
a = T("atlas_v1_summary.csv"); ref = a[(a.role == "reference") & (a.study != "TOTAL reference")].copy()
names = {"okada": "Okada 2026", "jessa_scrna": "Jessa 2022 (sc/snRNA)", "jessa_multiome": "Jessa 2022 (multiome)",
         "collot": "Collot 2026", "zang": "Zang 2024", "desisto": "DeSisto 2024"}
rows = []
for _, r in ref.sort_values("cells", ascending=False).iterrows():
    rows.append(f'{names.get(r.study, r.study)} & {int(r.samples)} & {int(r.patients)} & {int(r.cells):,} & '
                f'{esc(r.assay)} & {esc(r.site)} & {esc(r.h3_variants)} & {r.malignant_fraction_consensus:.3f} \\\\')
tot = ref[["samples", "patients", "cells"]].sum()
write("table1_cohort.tex", r"""\begin{table}[htbp]
\centering\footnotesize
\caption{\textbf{Reference cohort.} Samples, patients and cells per contributing study after quality control, with
assay, anatomical site and reported H3 variant (sample counts), and the fraction of cells called malignant by the
consensus rule. Counts are cells for whole-cell and nuclei for single-nucleus and multiome assays.}
\label{tab:cohort}
\resizebox{\textwidth}{!}{%
\begin{tabular}{lrrrlllr}
\toprule
Study & Samples & Patients & Cells & Assay (cells) & Site (samples) & H3 variant (samples) & Malignant \\
\midrule
""" + "\n".join(rows) + r"""
\midrule
Total & """ + f"{int(tot.samples)} & {int(tot.patients)} & {int(tot.cells):,}" + r""" & \multicolumn{4}{l}{18,882 genes shared by all studies} \\
\bottomrule
\end{tabular}}
\end{table}""")

# ---- Table 2: malignant-call validation ----
v = T("malignant_call_validation_summary.csv")
lbl = dict(names, liu="Liu 2022 (query)", labelle="LaBelle (query)", filbin="Filbin 2018 (query)", all=r"\textbf{All}")
rows = []
for _, r in v.iterrows():
    b = r.study == "all"
    f = (lambda s: r"\textbf{" + s + "}") if b else (lambda s: s)
    rows.append(" & ".join([f(esc(lbl.get(r.study, r.study))), f(f"{int(r.n_labelled):,}"), f(str(int(r.n_samples))),
        f(f"{r.pooled_agreement:.3f}"), f(f"{r.pooled_precision:.3f}"), f(f"{r.pooled_recall:.3f}"),
        f(f"{r.sample_median:.3f} ({r.sample_q25:.3f}--{r.sample_q75:.3f})"), f(f"{r.sample_mean:.3f}")]) + r" \\")
write("table2_cnv.tex", r"""\begin{table}[htbp]
\centering\footnotesize
\caption{\textbf{Agreement of malignant calls with published annotations.} Tier-2 comparison of the copy-number call
against each study's own annotation, for every annotated cell. Pooled values weight cells equally; the sample median
(with interquartile range) and unweighted sample mean weight samples equally.}
\label{tab:cnv}
\resizebox{\textwidth}{!}{%
\begin{tabular}{lrrrrrlr}
\toprule
Study & Cells & Samples & Agreement & Precision & Recall & Sample median (IQR) & Sample mean \\
\midrule
""" + "\n".join(rows) + r"""
\bottomrule
\end{tabular}}
\end{table}""")

# ---- Table 3: metaprograms ----
mp = T("metaprograms.csv"); fam = {"OPC_like": "OPC-like", "AC_like": "AC-like", "OC_like": "OC-like",
    "MES_like": "MES-like", "cycling": "cycling", "hypoxia": "hypoxia", "novel_MP9": "ciliated (MP9)"}
tech = {"technical:ribosomal": "ribosomal", "technical:mitochondrial": "mitochondrial",
        "technical:IEG_dissociation_stress": "IEG / dissociation stress"}
gcol = [c for c in mp.columns if "gene" in c.lower()][0]
rows = []
for _, r in mp.iterrows():
    f = fam.get(r.family) if isinstance(r.family, str) else r"\textit{" + tech.get(r.technical_class, "technical") + "}"
    genes = ", ".join(str(r[gcol]).replace(";", ",").split(",")[:6])
    rows.append(f'{r.metaprogram} & {f} & {int(r.n_samples)} & {"yes" if r.replicated else "no"} & '
                r"\textit{" + esc(genes) + r"} \\")
write("table3_metaprograms.tex", r"""\begin{table}[htbp]
\centering\footnotesize
\caption{\textbf{Consensus metaprograms.} Metaprograms at the primary dendrogram cut, their assigned lineage family or
technical class, the number of samples in which a contributing programme was found, whether they satisfy the
pre-registered replication rule (at least three studies and both assay classes), and their top genes.}
\label{tab:mp}
\resizebox{\textwidth}{!}{%
\begin{tabular}{llrll}
\toprule
Metaprogram & Family & Samples & Replicated & Top genes \\
\midrule
""" + "\n".join(rows) + r"""
\bottomrule
\end{tabular}}
\end{table}""")

# ---- Table 4: spatial co-localisation ----
s = T("spatial_colocalisation.csv")
h = s[(s.metric == "proportion") & (s.subset == "all") & (s.pair_type == "state_vs_TME")].sort_values("q_bh_state_vs_TME")
nm = lambda c: (c.replace("malignant_", "").replace("novel_MP9", "ciliated (MP9)").replace("_like", "-like")
                 .replace("_", " ").replace("T NK", "T/NK"))
rows = [f'{nm(r.class_a)} & {nm(r.class_b)} & {r.rho_re:.3f} & {r.rho_re_lo:.3f} to {r.rho_re_hi:.3f} & '
        f'{r.q_bh_state_vs_TME:.3g} & {int(r.n_sections_positive)}/{int(r.n_sections)} & {r.I2:.3f} \\\\'
        for _, r in h.head(10).iterrows()]
write("table4_spatial.tex", r"""\begin{table}[htbp]
\centering\footnotesize
\caption{\textbf{Spatial co-localisation of malignant states with microenvironment lineages.} Spot-level Spearman
correlations of cell2location proportions, Fisher-$z$ averaged within patient and pooled by random effects over 12 DMG
patients and 21 sections, Benjamini--Hochberg corrected across 42 state-by-microenvironment pairs. The ten strongest
pairs are shown; $I^2$ quantifies between-patient heterogeneity.}
\label{tab:spatial}
\resizebox{\textwidth}{!}{%
\begin{tabular}{llrlrrr}
\toprule
Malignant state & Lineage & $\rho$ & 95\% CI & $q$ & Sections $>0$ & $I^2$ \\
\midrule
""" + "\n".join(rows) + r"""
\bottomrule
\end{tabular}}
\end{table}""")

# ---- Supplementary Table 1: deposit inventory ----
d = pd.read_csv(R / "docs" / "datasets" / "dmg_dataset_selection.csv")
d = d[d.recommended_role.notna()]
def wrap(x, n):
    x = str(x)
    return esc(x[:n] + ("\\ldots" if len(x) > n else ""))
# Which screened deposits actually entered the released results (USED, built above from the pipeline outputs).
ACC_KEYS = [tok for _, acc, _ in USED.values()
            for tok in re.findall(r"GSE\d+|EGAS\d+|zenodo\.\d+", acc)]
def used_mark(r):
    acc = str(r.accession)
    return "yes" if (any(k in acc for k in ACC_KEYS) or "cellxgene" in str(r.repository).lower()) else "no"
rows = [f'{wrap(r.dataset, 30)} & {wrap(r.accession, 26)} & {wrap(r.repository, 14)} & {wrap(r.access, 16)} & '
        f'{wrap(r.n_cells_or_spots, 11)} & {used_mark(r)} & {wrap(r.recommended_role, 22)} \\\\'
        for _, r in d.iterrows()]
write("tableS1_deposits.tex", r"""\begin{table}[htbp]
\centering\scriptsize
\caption{\textbf{Supplementary Table 1 $\vert$ Public deposits screened and their role in this study.} Accession,
repository, access status and deposited cell or spot count for every deposit screened, with an explicit column
for whether it entered the released results. The final column is the role proposed at screening, which is not the
inclusion record: deposits marked ``no'' were dropped before or during quality control (the mouse developmental
reference and the adult Paired-Tag reference were never used, and the ScPCA samples failed quality control), and
the Phase~7 bulk cohort is outside the scope of this paper.}
\label{tab:deposits}
\setlength{\tabcolsep}{3pt}
\resizebox{\textwidth}{!}{%
\begin{tabular}{p{2.6cm}p{2.4cm}p{1.5cm}p{1.5cm}p{1.2cm}p{0.9cm}p{2.1cm}}
\toprule
Dataset & Accession & Repository & Access & Cells/spots & Used here & Role at screening \\
\midrule
""" + "\n".join(rows) + r"""
\bottomrule
\end{tabular}}
\end{table}""")

# ---- Key resources table (Cell Press STAR Methods requirement) ----
ver = pd.read_csv(R / "env" / "versions.csv")
SOFT = {"scanpy": "Wolf et al.", "scvi-tools": "Lopez et al.; Xu et al.", "cell2location": "Kleshchevnikov et al.",
        "infercnvpy": "Sturm et al. (inferCNV port)", "umap-learn": "McInnes et al.", "leidenalg": "Traag et al.",
        "scikit-learn": "Pedregosa et al.", "statsmodels": "Seabold and Perktold", "torch": "Paszke et al.",
        "anndata": "Virshup et al.", "pyro-ppl": "Bingham et al.", "python": "Python Software Foundation",
        "numpy": "Harris et al.", "scipy": "Virtanen et al.", "pandas": "McKinney", "matplotlib": "Hunter",
        "lifelines": "Davidson-Pilon"}
soft_rows = []
for _, r in ver.iterrows():
    if r["package"] not in SOFT:
        continue
    v = str(r["local"]) if pd.notna(r["local"]) and str(r["local"]) != "nan" else str(r["modal"])
    soft_rows.append(f'{esc(r["package"])} v{esc(v)} & {esc(SOFT[r["package"]])} & \\url{{https://pypi.org/project/{r["package"]}/}} \\\\')

atlas = T("atlas_v1_summary.csv"); atlas = atlas[atlas.study != "TOTAL reference"]
spatial = pd.read_csv(R / "results" / "qc" / "spatial_sections_used.csv")
used_keys = list(atlas.study) + sorted(set(spatial.study))
missing = [k for k in used_keys if k not in USED]
assert not missing, f"no accession recorded for: {missing}"
dep_rows = [f'{esc(USED[k][0])} & {esc(USED[k][1])} & {esc(USED[k][2])} \\\\' for k in used_keys]
print(f"key resources: {len(dep_rows)} deposits actually used ({len(atlas)} atlas + {len(set(spatial.study))} spatial)")

write("table_key_resources.tex", r"""\begin{table}[htbp]
\centering\scriptsize
\caption{\textbf{Key resources table.} Deposited data reused in this study and the software that produced the
results. Version numbers are those recorded by \texttt{scripts/28\_record\_environment.py}; the full environment
(557 conda and 167 pip entries for the local environment, plus the GPU image pins) is in \texttt{env/}.}
\label{tab:keyresources}
\setlength{\tabcolsep}{3pt}
\begin{tabular}{p{5.0cm}p{6.0cm}p{3.6cm}}
\toprule
\multicolumn{3}{l}{\textbf{Deposited data}} \\
Dataset & Identifier & Source \\
\midrule
""" + "\n".join(dep_rows) + r"""
\midrule
\multicolumn{3}{l}{\textbf{Software and algorithms}} \\
Software & Reference & Identifier \\
\midrule
""" + "\n".join(soft_rows) + r"""
\bottomrule
\end{tabular}
\end{table}""")
