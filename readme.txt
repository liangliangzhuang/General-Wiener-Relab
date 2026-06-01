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

