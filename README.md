
# Code_Pcvm_test

## Overview
[Asymptotic Distribution-Free Tests for Ultra-high Dimensional Parametric Regressions via Projected Empirical Processes and p-value Combination](https://arxiv.org/abs/2601.00541v1)
This repository contains the R source code to reproduce the Monte Carlo simulations and real data analyses presented in the manuscript and its Supplementary Material. This paper constructs a Cramer-von Mises type test based on a martingale-transformed, projected residual-marked empirical process, which renders the test asymptotically distribution-free. Furthermore, it proposes a novel hybrid test that aggregates empirical process-based and local smoothing tests using Cauchy combination, which is powerful against both low-frequency and high-frequency alternatives.

## Overview of R Scripts and Correspondence to Manuscript Tables

| Script Name | Description | Corresponding Table/Figure |
| :--- | :--- | :--- |
| `Pcvm_Pls_hybrid_linear.R; Pcvm_Pls_hybrid_linear_H12.R; Pcvm_Pls_hybrid_linear_H13.R` | Main simulation study for high-dimensional **linear regression** models ($H_{11}, H_{12}, H_{13}$). | [Table 1, Table S1, Table S2] |
| `Pcvm_Pls_hybrid_logit_H21.R; Pcvm_Pls_hybrid_logit_H22.R` | Main simulation study for high-dimensional **logistic regression** models ($H_{21}, H_{22}$). | [Table 2, Table S3] |
| `Pcvm_Pls_hybrid_time.R` | Computational runtime comparing our methods against competing resampling/projection tests. | [Table S10-S11] |
| `Pcvm_Pls_hybrid_logit_h.R` | Bandwidth sensitivity analysis ($0.5 h_0, 1.0 h_0, 2.0 h_0$) for the martingale transformation in logistic models. | [Table S4-S9] |
| `Pcvm_Pls_hybrid_linear_H11_t0.R` | Integration limit ($t_0$) sensitivity analysis (0.99, 0.95, and 0.90 quantiles) for linear models under $H_{11}$. | [Table S13-S15] |
| `Pcvm_Pls_hybrid_logit_H21_t0.R` | Integration limit ($t_0$) sensitivity analysis (0.99, 0.95, and 0.90 quantiles) for logistic models under $H_{21}$. | [Table S16-S18] |
| `Pcvm_Pls_hybrid_logit_H21_alpha.R` | Collinearity diagnosis between estimated $\hat{\beta}$ and residual-based SDR projections ($\hat{\alpha}$) before orthogonalization. | [ ] |
| `real_data_analysis.R` | Real data applications on the Communities and Crime dataset (linear) and OHSU AML dataset (logistic). | [Figure 1, Figure 2] |
