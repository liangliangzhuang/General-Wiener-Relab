
# Code and Data Description

This repository contains the real-data analysis code for the manuscript submitted to IISE Transactions:

**Modeling Multivariate Degradation with Time-Varying Mean–Variance Dynamics for Reliability Assessment: A Hierarchical Bayesian Approach**

The `Real-data` folder includes two real case studies used in the manuscript: fatigue crack growth data and permanent magnet brake data. Each case study contains the raw data, exploratory data analysis scripts, main analysis scripts, and generated results.

---

# Folder Structure

```text
Real-data/
├── crack/
│   ├── eda_crack_growth.R
│   ├── main_crack_growth.R
│   ├── crack_growth_data.csv
│   └── result/
│
└── PMB/
    ├── eda_PMB.R
    ├── main_PMB_linear.R
    ├── utility_PMB.R
    ├── PMB_data.csv
    └── result/
```

---

# Case Study 1: Fatigue Crack Growth Data

Folder: `Real-data/crack/`

- `eda_crack_growth.R`  
  Performs exploratory data analysis for the fatigue crack growth data.

- `main_crack_growth.R`  
  Runs the main model-fitting procedures and generates the results reported in the manuscript.

- `crack_growth_data.csv`  
  Raw fatigue crack growth data.

- `result/`  
  Stores generated numerical results, diagnostic outputs, and figures.

---

# Case Study 2: Permanent Magnet Brake Data

Folder: `Real-data/PMB/`

- `eda_PMB.R`  
  Performs exploratory data analysis for the permanent magnet brake data.

- `main_PMB_linear.R`  
  Runs the main linear-model analysis and generates the results reported in the manuscript.

- `utility_PMB.R`  
  Contains utility functions for:
  
  - data loading,
  - model specification,
  - Bayesian inference,
  - MCMC diagnostics,
  - posterior analysis,
  - plotting utilities.

- `PMB_data.csv`  
  Raw permanent magnet brake degradation data.

- `result/`  
  Stores generated numerical results, diagnostic outputs, and figures.

---

# Notes

The scripts are organized so that exploratory analysis and main model-fitting procedures are separated for clarity.

All generated figures, tables, posterior summaries, and diagnostic outputs are saved in the corresponding `result/` folders.

## Reproducibility Workflow

The following table summarizes the data files, code files, expected outputs, and approximate runtime for reproducing the main empirical results.

| Which results to reproduce | Data File | Code File | Expected output | Run time at the above-specified computer conditions |
|---|---|---|---|---|
| Table 1 | `PMB_data.csv` | `eda_PMB.R` | `1-pmb_path.pdf`<br>`2-pmb_ratio.pdf`<br>`3-pmb_corr.pdf`<br>`4-pmb_hetero.pdf` | 15 seconds |
| Table 4, Table 5, Figure 7 | `PMB_data.csv` | `main_PMB_linear.R` | `PMB_GoF_QQ.pdf`<br>`PMB_Path_Fitting.pdf`<br>`Final_PMB_Model_all.RData` | 2 minutes |
| Figure 9 | Embedded in script | `eda_crack_growth.R` | `1-crack.pdf`<br>`2-ratio.pdf` | 10 seconds |
| Table 6, Figure 10 | `crack_growth_data.csv` | `main_crack_growth.R` | `Goodness_of_Fit_QQ.pdf`<br>`Path_Fitting_Best_Model.pdf` | 30 seconds |

