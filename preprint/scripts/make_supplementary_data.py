"""Write the preprint's Supplementary Data workbooks from existing tables (run: python preprint/scripts/make_supplementary_data.py).

  Supplementary_Data_1.xlsx  per-sample metadata of the reference and query cohorts (data/processed/sample_manifest.csv,
                             local file paths dropped)
  Supplementary_Data_2.xlsx  metaprograms: summary, top-50 genes, member programmes, replication matrix
  Supplementary_Data_3.xlsx  malignant-call agreement with author labels, per sample and per study
  Supplementary_Data_4.xlsx  spatial co-localisation, all pairs and metrics
  Supplementary_Data_5.xlsx  robustness analyses (leave-one-publication-out classifier, metaprogram recovery, composition
                             label sensitivity)
  Supplementary_Data_6.xlsx  composition models: every category x contrast, estimability, per-sample proportions
Each workbook starts with a README sheet naming the source table of every sheet.
"""
import pathlib
import pandas as pd

R = pathlib.Path(__file__).resolve().parents[2]
T, Q = R / "results/tables", R / "results/tables/quality"
OUT = R / "preprint/supplementary_data"; OUT.mkdir(exist_ok=True)
BOOKS = {
    1: ("Per-sample metadata", {"samples": R / "data/processed/sample_manifest.csv"}),
    2: ("Consensus metaprograms", {"metaprograms": T / "metaprograms.csv", "metaprogram_genes": T / "metaprogram_genes.csv",
                                   "members": T / "metaprogram_members.csv", "replication_matrix": T / "replication_matrix.csv",
                                   "cut_sensitivity": T / "metaprogram_cut_sensitivity.csv"}),
    3: ("Malignant-call agreement with author labels", {"per_sample": T / "malignant_call_validation.csv",
                                                        "per_study": T / "malignant_call_validation_summary.csv"}),
    4: ("Spatial co-localisation", {"colocalisation": T / "spatial_colocalisation.csv",
                                    "per_section": T / "spatial_colocalisation_per_section.csv"}),
    5: ("Robustness analyses", {"loso_summary": Q / "loso_classifier_summary.csv", "loso_per_class": Q / "loso_classifier_per_class.csv",
                                "loso_per_patient": Q / "loso_classifier_per_patient.csv", "loso_confusion": Q / "loso_classifier_confusion.csv",
                                "mp_recovery_summary": Q / "metaprogram_lopo_summary.csv", "mp_scenarios": Q / "metaprogram_lopo_scenarios.csv",
                                "composition_compare": Q / "composition_label_sensitivity_compare.csv"}),
    6: ("Composition models", {"results": T / "composition_results.csv", "estimability": T / "composition_estimability.csv",
                               "per_sample": T / "composition_per_sample.csv"}),
}
for n, (title, sheets) in BOOKS.items():
    readme = pd.DataFrame({"sheet": list(sheets), "source": [str(p.relative_to(R)) for p in sheets.values()]})
    with pd.ExcelWriter(OUT / f"Supplementary_Data_{n}.xlsx") as xw:
        pd.DataFrame({"Supplementary Data": [f"{n}: {title}"]}).to_excel(xw, sheet_name="README", index=False)
        readme.to_excel(xw, sheet_name="README", index=False, startrow=3)
        for name, p in sheets.items():
            d = pd.read_csv(p, low_memory=False)
            if name == "samples":
                d = d.drop(columns=["file_paths"], errors="ignore")
            d.to_excel(xw, sheet_name=name[:31], index=False)
    print("wrote", OUT / f"Supplementary_Data_{n}.xlsx")
