# Code for Online Updating Proportional Hazards Test

This repository contains `R` code for *An Online Updating Approach
for Testing the Proportional Hazards Assumption with Streams of
Survival Data* by Yishu Xue, HaiYing Wang, Jun Yan and Elizabeth D. Schifano.

1. `bczph.cee.R` implements the proposed online updating cumulative
test statistic with the matrices and Schoenfeld residuals evaluated
at the cumulative estimating equation (CEE) estimators of beta.

2. `bczph.cuee.R` implements the proposed online updating cumulative
test statistic with the matrices and Schoenfeld residuals evaluated
at the cumulatively updated estimating equation (CUEE) estimators of beta.

3. `bwzph.cee.R` implements the proposed online updating window  test statistic
with the matrices and Schoenfeld residuals evaluated at the CEE estimators
of beta.

4. `simulate.survival.data.R` contains two functions to simulate survival data,
one for the proportional hazards model with a constant baseline hazard function,
and another for the frailty model.

5. `example.R` provides an example how online updating can be performed with
simulated data whose structure resembles that of SEER data. 

6. `simulated_survdata.csv` contains a total of 1,000,000 simulated observations
with three covariates, `Age`, `Sex` and `Black`. The average censoring rate
is 50.2%.
