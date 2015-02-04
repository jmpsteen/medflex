# medflex 0.5-0

## Major changes
* Robust sandwich estimator for standard errors as alternative for bootstrap SEs (option 'se = "robust"' in neModel function)
* Advanced options for effect decomposition (neEffdecomp function now allows to specify reference exposure levels and covariate levels via 'xRef' and 'covLev' arguments, respectively)

## Minor changes
* Additional option to weight for multicategorical variables (either via inverse-probability-of-treatment weighting for multivariate exposures or via ratio-of-mediator-probability weighting for multivariate mediators)
* ratio-of-mediator-probability weights can be extracted from the expanded dataset object via weights function
