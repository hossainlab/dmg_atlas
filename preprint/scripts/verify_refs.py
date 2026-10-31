#!/usr/bin/env python
"""Verify every references/references.bib entry against its Crossref record (fetched by DOI).

For each entry: GET https://api.crossref.org/works/<doi> (cached in the system temp dir, not in the repo) and compare
title (normalised), year (print or online), volume, first-author family name and journal. Prints mismatches, uncited
entries and unresolved \\cite keys.
Run: python preprint/scripts/verify_refs.py preprint/references.bib preprint/dmg_preprint.tex
"""
import sys, re, json, time, tempfile, unicodedata, urllib.request, urllib.parse
from pathlib import Path

BIB = Path(sys.argv[1]); TEX = Path(sys.argv[2]) if len(sys.argv) > 2 else None
OUT = Path(tempfile.gettempdir()) / "dmg_atlas_crossref"; OUT.mkdir(parents=True, exist_ok=True)


def norm(s):
    s = re.sub("[\u2010-\u2015\u2212]", "-", s or "")
    s = unicodedata.normalize("NFKD", re.sub(r"<[^>]+>", "", s)).encode("ascii", "ignore").decode().lower()
    return re.sub(r"[^a-z0-9]+", " ", s).strip()


def field(body, name):
    m = re.search(rf"\b{name}\s*=\s*\{{(.*?)\}}\s*,?\s*(?:\n|$)", body, re.S | re.I)
    return re.sub(r"[{}]", "", m.group(1)).strip() if m else ""


entries = re.findall(r"@(\w+)\{([^,]+),(.*?)\n\}", BIB.read_text(), re.S)
rows = []
for typ, key, body in entries:
    doi = field(body, "doi"); f = OUT / f"{key}.json"
    if not doi:
        rows.append(dict(key=key, doi="", status="NO_DOI")); continue
    if not f.exists():
        req = urllib.request.Request(f"https://api.crossref.org/works/{urllib.parse.quote(doi)}",
                                     headers={"User-Agent": "dmg-atlas-preprint/1.0 (reference verification)"})
        try:
            f.write_text(urllib.request.urlopen(req, timeout=30).read().decode()); time.sleep(0.2)
        except Exception as e:
            rows.append(dict(key=key, doi=doi, status=f"FETCH_FAIL {e}")); continue
    m = json.loads(f.read_text())["message"]
    cr_title = (m.get("title") or [""])[0]
    years = {str(m[k]["date-parts"][0][0]) for k in ("published-print", "published-online", "issued") if m.get(k)}
    cr_first = norm((m.get("author") or [{}])[0].get("family", ""))
    cr_journal = norm((m.get("container-title") or [""])[0]); cr_short = norm((m.get("short-container-title") or [""])[0])
    bib_author = field(body, "author").split(" and ")[0]; bib_first = norm(bib_author.split(",")[0])
    issues = []
    if norm(field(body, "title")) != norm(cr_title):
        issues.append(f"title: bib '{field(body, 'title')[:70]}' vs crossref '{cr_title[:70]}'")
    if field(body, "year") not in years:
        issues.append(f"year: bib {field(body, 'year')} vs crossref {sorted(years)}")
    if field(body, "volume") and m.get("volume") and field(body, "volume") != m.get("volume"):
        issues.append(f"volume: bib {field(body, 'volume')} vs crossref {m.get('volume')}")
    if bib_first and cr_first and bib_first != cr_first:
        issues.append(f"first author: bib {bib_first} vs crossref {cr_first}")
    bj = norm(field(body, "journal") or field(body, "booktitle"))
    if typ.lower() == "article" and bj and bj not in (cr_journal, cr_short) and not cr_journal.startswith(bj) and not bj.startswith(cr_short or "@"):
        issues.append(f"journal: bib '{bj}' vs crossref '{cr_journal}'")
    rows.append(dict(key=key, doi=doi, status="OK" if not issues else "CHECK", issues=" | ".join(issues), crossref_title=cr_title))

cited = set()
if TEX:
    for g in re.findall(r"\\cite[pt]?\*?\{([^}]+)\}", TEX.read_text()):
        cited |= {k.strip() for k in g.split(",")}
keys = {r["key"] for r in rows}
ok = sum(r["status"] == "OK" for r in rows)
print(f"{len(rows)} entries: {ok} OK, {len(rows) - ok} need attention")
for r in rows:
    if r["status"] != "OK":
        print(f"  {r['key']}: {r['status']} {r.get('issues', '')}")
if TEX:
    print("uncited bib entries:", sorted(keys - cited) or "none")
    print("unresolved cite keys:", sorted(cited - keys) or "none")
