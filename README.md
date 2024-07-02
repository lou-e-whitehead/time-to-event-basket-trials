# time-to-event-basket-trials
Bayesian borrowing in basket trials with time to event outcomes

JAGS code together with R functions to implement methods detailed in 'Bayesian borrowing for basket trials with time-to-event outcomes' (preprint available soon).

All code is based upon code from: https://github.com/BasketTrials/Bayesian-analysis-models and https://github.com/longitudinal-basket-trials/Bayesian-analysis-methods

Files contained in this repository can be used to reproduce simulation results reported in the paper entitled: Bayesian borrowing for basket trials with time-to-event outcomes (Whitehead et al, 2024 TBC - preprint available soon)

The files "BHM-W.txt", "EXNEX-W.txt" and "Strat-W.txt" are for implementing the proposed Bayesian models through JAGS and may be used for analysing basket trial data with time-to-event outcomes.

In particular, "BHM-W.txt" is the standard Bayesian hierarchical model which assumes the subgroup-specific parameters to be exchangeable.

"EXNEX-W.txt" is the borrowing method proposed by Neuenschwander et al. (2016), implemented for time-to-event endpoints.

"Strat-W.txt" is the approach of no borrowing or independent Bayesian analysis in each subgroup.

The motivating example reported in Section 3 of the manuscript can be reproduced by calling these models in R with the R2JAGS package. Example R code to implement each of the models via R2JAGS has been provided in the file "Implementation.R".

Contact: l.whitehead2@ncl.ac.uk
