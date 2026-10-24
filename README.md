# Depth and durability of response in patient-derived xenografts

Code for a re-analysis of the Novartis PDX Encyclopedia (Gao et al., 2015). It measures durability (regrowth after a tumor's smallest volume) alongside depth of response, and tests whether the early growth curve forecasts later volume and selects treatments better than depth.

## Repository layout

```
├── Analysis_Notebook.ipynb   every analysis (Sections 1–15); writes outputs/
├── Figures_Notebook.ipynb    every figure, drawn only from outputs/
├── pdxe_utils.py             shared endpoint functions, with a self-test
├── prespecification/         the three analysis plans fixed before testing
├── environment.yml           conda environment (requirements.txt for pip)
├── CITATION.cff
└── LICENSE
```

`data/` and `outputs/` are created locally and not tracked.

## Data

| Input | Path under `data/` | Access |
|---|---|---|
| Gao et al. (2015) Supplementary Table 1 (required) | `gao2015/41591_2015_BFnm3954_MOESM10_ESM.xlsx` | Public ([doi:10.1038/nm.3954](https://doi.org/10.1038/nm.3954)) |
| PDXNet tumor volumes (Section 13) | `pdxnet/PDXNet-PDTC-Export-*.csv` | Public ([portal](https://pdxnet.shinyapps.io/portal/)) |
| Oberg et al. (2021) ovarian cohort (Section 13) | `oberg2021/supp2_data.csv` | From the authors on request |
| MSigDB hallmark gene sets (Section 9.3) | `msigdb/` | Public; downloaded if absent |
| Trials A6181087, CA031, EFC10261 (Section 14) | `project_data_sphere/<trial>/` | Restricted ([Project Data Sphere](https://data.projectdatasphere.org)) |

Sections whose optional input is absent print a skip line. Trial files and the patient-level files Section 14 writes fall under a data use agreement and must not be committed.

## Prespecified plans

| File | Analysis | SHA-256 prefix |
|---|---|---|
| `mechanism_family_plan.json` | Mechanism class (Section 11) | `7a013fd2991e9848` |
| `clinical_single_trial_plan.json` | A6181087 clinical plan (Section 14.3) | `971249c51d66f110` |
| `clinical_pooled_plan.json` | Pooled clinical plan (Section 14.4), written after the single-trial primary result | `2a34fec1d3aef86d` |

The notebook stops if a plan is missing or its checksum (line endings normalized to LF) differs from the recorded prefix. The forecasting endpoints and margins (Section 7) were set before that analysis first ran; later additions are labelled as added in the notebook.

## Reproducing the analysis

1. `conda env create -f environment.yml`, `conda activate pdxe`, `python -m ipykernel install --user --name pdxe`.
2. Place the inputs under `data/`.
3. Run `Analysis_Notebook.ipynb` with *Restart and run all* (about 60 minutes on a CPU; 282 consistency checks), then `Figures_Notebook.ipynb` (about 6 minutes). Every reported value is written to `outputs/results_register.json`. All seeds are fixed (project seed 20260905).

## Main results

Intervals are 95% and resample PDX models. Day-24 volume is forecast from readings up to day 10; the depth set is baseline volume and best average response, the trajectory set adds eight shape descriptors.

| Analysis | Result | Section |
|---|---|---|
| Depth vs durability, Spearman ρ (1,832 post-nadir arms) | 0.355 (0.301 to 0.406) | 4 |
| Depth-matched pairs differing in durability (noise alone: 4.2%) | 27.6% | 4 |
| Forecasting error, trajectory vs depth (prespecified; margin 10%) | 19.6% (13.1 to 26.0) | 7 |
| Selection regret, trajectory vs depth (prespecified; margin 15%) | 44.3% (20.7 to 61.9) | 7 |
| Forecasting, trajectory vs expanded depth | 10.2% (3.5 to 17.0) | 8 |
| Forecasting, trajectory vs baseline plus slope | 4.2% (−1.2 to 9.3) | 8 |
| Forecasting, trajectory vs baseline plus slope by least squares | −9.7% (−19.0 to −1.5) | 8 |
| Regret, trajectory vs depth by least squares | 19.5% (−22.8 to 47.4) | 8 |
| Forecasting, all treated arms | 14.5% (10.5 to 18.3) | 8 |
| Combinations better than best component on both endpoints | 21.6% | 10 |
| Cytostatic minus cytotoxic durability, per day (planned permutation P = 0.11) | 0.0082 (0.0009 to 0.0158) | 11 |
| Patients, C-index change, depth plus shape vs depth | −0.004 (−0.040 to 0.030) | 14 |

## Citation

Code metadata are in [`CITATION.cff`](CITATION.cff). Please also cite the data sources: Gao et al., *Nat. Med.* 2015; Koc et al., *NAR Cancer* 2022; Oberg et al., *Sci. Rep.* 2021; Liberzon et al., *Cell Syst.* 2015.

This work uses information obtained from Project Data Sphere. Neither Project Data Sphere nor the owners of that information contributed to, approved or are responsible for this work.

## License

MIT. Input data remain under their sources' terms.
