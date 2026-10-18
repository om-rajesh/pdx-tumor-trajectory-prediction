# Depth and durability of response in patient-derived xenografts

Analysis and figure code for a re-analysis of the Novartis PDX Encyclopedia (Gao et al., 2015). Preclinical studies usually score a treated tumor's growth curve by its depth of response (how far it shrank). This project also measures durability (how quickly it regrows after its smallest volume, the nadir), and tests whether the early curve's shape forecasts later volume and chooses between treatments better than depth alone.

## Repository layout

```
├── Analysis_Notebook.ipynb    all analyses (Sections 1–15)
├── Notebook_Figures.ipynb     all figures, drawn from outputs/
├── pdxe_utils.py              shared endpoint functions, with a self-test
├── prespecification/          analysis plans fixed before testing
├── data/                      inputs (not in the repository)
├── outputs/                   results register, manifest and tables
├── environment.yml            conda environment (requirements.txt for pip)
├── CITATION.cff
└── LICENSE
```

## Data

Only the Gao et al. (2015) workbook is required: 67,445 tumor-volume readings in 4,706 arms (one PDX model given one treatment, in one mouse) from 277 models of six tumor types. A section whose optional input is missing prints a skip line.

| Input | Path under `data/` | Access |
|---|---|---|
| Gao et al. (2015) workbook (required) | `gao2015/41591_2015_BFnm3954_MOESM10_ESM.xlsx` | Public |
| PDXNet tumor-volume export (Section 13) | `pdxnet/PDXNet-PDTC-Export-*.csv` | Public ([portal](https://pdxnet.shinyapps.io/portal/)) |
| Oberg et al. (2021) ovarian cohort (Section 13) | `oberg2021/supp2_data.csv` | From the authors on request |
| MSigDB hallmark gene sets (Section 9.3) | `msigdb/*.gmt` | Public |
| Trials A6181087, CA031, EFC10261 (Section 14) | `project_data_sphere/<trial>/` | Restricted ([Project Data Sphere](https://data.projectdatasphere.org)) |

The trial files and the patient-level files Section 14 writes are covered by a data use agreement and are not in this repository.

## Main results

Intervals are 95%; every value is in `outputs/results_register.json`. Forecasting predicts day-24 volume from measurements up to day 10 with a depth set (baseline volume, best average response) or a trajectory set (depth plus eight shape descriptors); percentages are relative reductions in root-mean-square error.

| Analysis | Result | Section |
|---|---|---|
| Depth against durability, Spearman ρ (1,832 post-nadir arms) | 0.355 | 4 |
| Depth-matched treatment pairs that differ in durability | 27.7% | 4 |
| Mechanism discrimination accuracy, unpaired / paired | 0.548 / 0.557 | 5 |
| Forecasting, trajectory over depth (prespecified; margin 10%) | 19.6% (13.1 to 26.0) | 7 |
| Treatment selection, reduction in regret (prespecified; margin 15%) | 44.3% (20.7 to 61.9) | 7 |
| Pairwise ranking accuracy, depth / trajectory | 0.639 / 0.706 | 7 |
| PDX models where the two sets choose differently | 50.8% | 7 |
| Trajectory over expanded depth | 10.2% (3.1 to 17.1) | 8 |
| Trajectory over depth plus current volume | 3.5% (−3.0 to 10.1) | 8 |
| Noise-padded depth over depth (control) | 1.9% (−2.9 to 5.9) | 8 |
| Without the post-nadir restriction | 14.5% (10.5 to 18.3) | 8 |
| Untreated arms, day 14 to 21 (negative control) | 16.8% (7.2 to 25.4) | 8 |
| Ranking accuracy, pairs more than 8 noise SDs apart, depth / trajectory | 0.714 / 0.906 | 8 |
| Molecular dissociation statistic D | −0.054 (−0.125 to 0.032) | 9 |
| Combination evaluations better than the best component on both endpoints | 26.7% | 10 |
| Cytostatic minus comparator families, durability (per day) | 0.0082 (0.0009 to 0.0158) | 11 |
| Regrowth-model separation, fold change from 28 to 84 days (half-life 0.25 d) | 6.90 (4.24 to 9.95) | 12 |
| PDXNet between-animal SD / source-study between-model SD | 0.495 | 13 |
| Patients, change in C-index, depth plus shape over depth | −0.004 (−0.040 to 0.030) | 14 |

## Reproducing the analysis

1. `conda env create -f environment.yml`, `conda activate pdxe`, then `python -m ipykernel install --user --name pdxe` (or `pip install -r requirements.txt` on Python 3.11).
2. Run `Analysis_Notebook.ipynb` with *Restart and run all* (about 30 minutes on a CPU), then `Notebook_Figures.ipynb` (about 6 minutes).

## Citation

Code metadata are in [`CITATION.cff`](CITATION.cff). Please also cite the data sources:

- Gao, H. *et al.* High-throughput screening using patient-derived tumor xenografts to predict clinical trial drug response. *Nat. Med.* **21**, 1318–1325 (2015). https://doi.org/10.1038/nm.3954
- Koc, S. *et al.* PDXNet portal: patient-derived xenograft model, data, workflow and tool discovery. *NAR Cancer* **4**, zcac014 (2022). https://doi.org/10.1093/narcan/zcac014
- Oberg, A. L. *et al.* Statistical analysis of comparative tumor growth repeated measures experiments in the ovarian cancer patient derived xenograft (PDX) setting. *Sci. Rep.* **11**, 8076 (2021).
- Liberzon, A. *et al.* The Molecular Signatures Database hallmark gene set collection. *Cell Syst.* **1**, 417–425 (2015).

This work uses information obtained from Project Data Sphere. Neither Project Data Sphere nor the owners of that information contributed to, approved or are responsible for this work.

## License

MIT. The input data remain under their sources' terms.
