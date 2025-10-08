# MI-Checklist (condensed)

Source: https://www.nature.com/articles/s41591-020-1041-y/tables/1

| Item | Answer | Where |
| --- | --- | --- |
| Study design stated | Yes | REPORT.md: Introduction |
| Data origin described | Yes | REPORT.md: Methods (Data) |
| Inclusion/exclusion criteria | Yes | REPORT.md: Methods (Preprocessing filters) |
| Outcome/targets defined | Yes | REPORT.md: Introduction, Results |
| Cohort splits | N/A | Single dataset; no train/test split |
| Sample size justification | N/A | Unsupervised clustering |
| Handling of missing data | N/A | Expression matrix; cells/genes filtered |
| Preprocessing steps | Yes | REPORT.md: Methods (Preprocessing) |
| Model description | Yes | Graph clustering (Louvain/Leiden), k-means |
| Hyperparameters | Yes | REPORT.md: Methods (resolutions, neighbors, PCs) |
| Evaluation metrics | Yes | ARI, Silhouette (Results) |
| Statistical uncertainty | Partial | Variation across algorithms reported |
| External validation | N/A | No external dataset |
| Code availability | Yes | REPORT.md: Code availability |
| Data availability | Partial | `.h5ad` not included; path provided in README |
| Hardware/software | Yes | README: Python 3.10 + libs |
| Interpretability | Yes | Marker genes and dot/violin plots |
| Reproducibility | Yes | run_all.py, build_report.py |
| Limitations | Partial | Discussed in Results/Conclusion |
| Ethical considerations | N/A | Non-human, secondary analysis |

Date: 2025-10-08
