### Title: Associations of the LIBRA index with cognitive resilience to genetic susceptibility to dementia 

The R code provided in this repository can be used to:

(i) model trajectories of global cognition using a latent process mixed model for multivariate longitudinal outcomes

(ii) assign CRgen status based on corrected individual slopes estimated in (i)

(iii) examine the association of LIBRA on CRgen status using a logistic regression model

(iv) do parametric bootstrap to account for the uncertainty in the corrected slopes estimated in the first stage of our modeling strategy

(v) represent a forest plot of Odds Ratios obtained in (iii) and the correct 95% confidence interval calculated using the total variance obtained by boostrap in (iv).


Figure. Multivariable-adjusted associations between the baseline LIfestyle for BRAin health (LIBRA) risk score, its individual components, and the odds of cognitive resilience in ApoE-ɛ4 carriers. 

Estimates are derived from two separate logistic regression models considering CRgen status as the outcome and the reversed scale (i.e., higher is brain-healthier) of the continuous LIBRA risk score (Model 1) or the twelve, unweighted, binary, LIBRA components (Model 2) as predictors of interest, adjusted for sex, age at baseline, education, and study center. 95% confidence intervals were obtained based on parametric bootstraps  with 1,000 replicates. For the components’ estimates, we ordered the protective/risk factors according to the strength of the associations. 

![img](ForestPlot.png) 
