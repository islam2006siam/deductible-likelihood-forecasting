# Modelling the deductible inside the likelihood: Density forecasting of dependent claim frequency and severity

**Eslam Abdelhakim Seyam**<br>
Department of Insurance and Risk Management, College of Business,<br>
Imam Mohammad Ibn Saud Islamic University (IMSIU), Riyadh 13318, Saudi Arabia<br>
Email: [isiam@imamu.edu.sa](mailto:isiam@imamu.edu.sa)

This repository contains the R code and data needed to reproduce the empirical
results, simulation study and figures of the paper.

## Contents

| Path | Description |
|---|---|
| `analysis.R` | Complete analysis in a single script |
| `data/Insample.csv` | LGPIF policy-year file, 2006-2010 (5,639 records) |
| `data/CLAIMLEVEL.csv` | LGPIF claim-level file, building-and-contents claims, 2006-2010 (6,258 claims) |

## Data

The data come from the Wisconsin Local Government Property Insurance Fund
(LGPIF) and are distributed with the *Loss Data Analytics* open textbook:

- `Insample.csv`: https://raw.githubusercontent.com/OpenActTexts/Loss-Data-Analytics/master/Data/Insample.csv
- `CLAIMLEVEL.csv`: https://raw.githubusercontent.com/OpenActTexts/Loss-Data-Analytics/master/Data/CLAIMLEVEL.csv

| File | MD5 |
|---|---|
| `Insample.csv` | `14ced05b9af8e70045f2796c57528bf6` |
| `CLAIMLEVEL.csv` | `bcfd1179a285a3cfe83a53419df1862c` |

The 2011 policy year is not publicly available, so the models are trained on
2006-2009, 2010 is the held-out year, and 2009 is the calibration year for the
conformal intervals.

## Requirements

- R 4.3.3
- rstan 2.32.6 (Stan 2.32.2) with a working C++ toolchain (Rtools 4.3 on Windows)
- loo 2.7.0, data.table 1.18.4, MASS 7.3-60, scoringRules 1.1.3,
  ggplot2 3.5.1, patchwork 1.3.0, scales 1.3.0

```r
install.packages(c("rstan", "loo", "data.table", "MASS", "scoringRules",
                   "ggplot2", "patchwork", "scales"))
```

## Running the analysis

Set the working directory to the repository root and run

```
Rscript analysis.R
```

The script runs the following sections in order:

1. Packages and settings
2. Stan programs (joint frequency-severity model, severity-only model, simulation check model)
3. Functions
4. Data preparation
5. Model compilation and density checks
6. Severity margin selection (lognormal, gamma, Burr XII by PSIS-LOO; penalized ML fit of the GB2)
7. Joint frequency-severity models (M-A, M-B, M-C1, M-C2, M-D1, M-D2, M-B NB, M-D2 NB)
8. Two-part GLM benchmark (M-E)
9. Out-of-sample prediction for 2010 (log score, CRPS, PIT, Gini, conformal intervals, portfolio risk)
10. Pricing, ground-up recovery and dependence
11. Simulation study (24 scenario cells, 960 replications, 4 estimators each)
12. Simulation check against MCMC
13. Figures
14. Session information

All results are written to `outputs/`, which is created on the first run.
Fitted joint models are stored in `outputs/fits/` and simulation cells in
`outputs/simulation/`. Models and cells that already exist are skipped, so an
interrupted run can be resumed.

## Computing time

On a 4-core desktop the joint models take about 1-2 hours each, the severity
margin selection about 3 hours, the simulation study several hours and the
MCMC check about 40 minutes per replicate. The remaining sections take
minutes.

## Reproducibility

The global seed is 20260703. Each section seeds `20260703 + k`, and the
simulation replications use seed `20260703 + 1000 * cell + replicate`.
