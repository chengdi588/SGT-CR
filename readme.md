# SGT-CR

Title: Inference and diagnostics for censored linear regression model with skewed generalized t distribution

Author: Chengdi Lian, Yaohua Rong, Jinwen Liang, Ruijie Guan and Weihu Cheng\*

In recent years, some interested data can be recorded only if the values fall within an interval range, and the responses are often subject to censoring. Attempting to perform effective statistical analysis with censored, especially heavy-tailed and asymmetric data, can be difficult. In this paper, we develop a novel linear regression model based on the proposed skewed generalized $t$ distribution for censored data. The likelihood-based inference and diagnostic analysis are established using the Expectation Conditional Maximization Either (ECME) algorithm in conjunction with smoothing approximate functions. We derive relavant measures to perform global influence for this novel model and develop local influence analysis based on the conditional expectation of the complete-data log-likelihood function. Some useful perturbation schemes are discussed. We illustrate the finite sample performance and the robustness of the proposed method by simulation studies. The proposed model is compared with other procedures based on a real dataset, and a sensitivity analysis is also conducted.

- "Use SGT.CR.applications.R to fit the SGT-CR model to real data."

- "Prior to usage, ensure that SGT.CR.basic.functions.R is loaded."

- "Use SGT.CR.simu.R for conducting simulations."
