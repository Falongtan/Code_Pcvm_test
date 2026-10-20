
# Code_Pcvm_test

## Overview
[Asymptotic Distribution-Free Tests for Ultra-high Dimensional Parametric Regressions via Projected Empirical Processes and p-value Combination](https://arxiv.org/abs/2601.00541v1).
This repository contains the R source code to reproduce the Monte Carlo simulations and real data analyses presented in the manuscript and its Supplementary Material. This paper constructs a Cramer-von Mises type test based on a martingale-transformed, projected residual-marked empirical process, which renders the test asymptotically distribution-free. Furthermore, it proposes a novel hybrid test that aggregates empirical process-based and local smoothing tests using Cauchy combination, which is powerful against both low-frequency and high-frequency alternatives.

## Overview of R Scripts and Correspondence to Manuscript Tables

| Script Name | Description | Corresponding Table/Figure |
| :--- | :--- | :--- |
| `Table1_Pcvm_Pls_hybrid_linear_H11.R`<br>`Table_S1_Pcvm_Pls_hybrid_linear_H12.R`<br>`Table_S2_Pcvm_Pls_hybrid_linear_H13.R` | Main simulation study for high-dimensional **linear regression** models ($H_{11}, H_{12}, H_{13}$). | [Table 1, Table S1, Table S2] |
| `Table2_Pcvm_Pls_hybrid_logit_H21.R`<br>`Table_S9_Pcvm_Pls_hybrid_logit_H22.R` | Main simulation study for high-dimensional **logistic regression** models ($H_{21}, H_{22}$). | [Table 2, Table S9] |
| `Table_S3-S5_Pcvm_Pls_hybrid_logit_H21_bandwidth.R`<br>`Table_S6-S8_Pcvm_Pls_hybrid_logit_H22_bandwidth.R` | Bandwidth sensitivity analysis ($0.5 h_0, 1.0 h_0, 2.0 h_0$) for the martingale transformation in logistic models under $H_{21}$ and $H_{22}$. | [Table S3-S5, Table S6-S8] |
| `Table34_Pcvm_Pls_hybrid_time.R` | Computational runtime comparing our methods against competing resampling/projection tests. | [Table 3, Table 4] |
| `real_data_analysis.R` | Real data applications on the Communities and Crime dataset (linear) and OHSU AML dataset (logistic). The corresponding datasets are provided in the `communities+and+crime/` and `aml_ohsu_2022/` directories, respectively. | [Figure 1, Figure 2] |
| `utils.R` | Shared utility functions across all simulations, including the empirical Brownian motion integral calculation (`pvalue_integ_Brown`) and GCV bandwidth selection (`bandwidth_choice`). | [ ] |
