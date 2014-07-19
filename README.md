psm
===

Code relating to propensity score matching


# Troubleshooting a poorly balanced PS Match #

## Preliminaries ##

### Presence of a Propensity Score Model ###

Prior to matching a propensity score was fit, by some GLiM such as logistic regression, and this PS is at least part of the matching criterion.  Assume that the PS model had model formula `z ~ x1 + <...> + xk`, in which some of the `x`es may themselves be interactions or transforms of other `x`es.  

### Meaning of "Poorly Balanced" ##

For now, this means the following.   Let the match, based wholly or in part on the PS model, be encoded in factor `m`.  The match is assumed to be poorly balanced in the sense that when conditional logistic regression is applied with formula  `z ~ x1 + <...> + xk + strata(m)`, the Rao score test fails to reject the hypothesis that coefficients on `x1`, .., `xk` are all 0.   (Ie, `RItools:xBalance` gives a non-significant p-value.)

## Possible PS Model Issues ##

### Likely suspects
- Leverage points in the PS fit
- substantial changes in PS coeffs if unmatched subjects are dropped from the fit
- ?Isolate offending observations, w/ a custom-developed residual plot? 

(To my knowledge no such residual plot has been proposed, but I have a pretty clear idea about how to get started.
  
1. Locate poorly balanced covariates.  Dig up the fitted PS model.  Restrict consideration to subjects included in match. 
2. For each poorly balanced covariate `xk`, calculate contributions to relevant score equation of the fitted PS model, ie $(z - \hat{\pi}) \cdot x_k$, where "$\cdot$" denotes element-wise product.
3. Write $\bar{z}_m$ for the vector encoding matched set means of $z$.  For each poorly balanced covariate `xk`, calculate $(z-\bar{z}_m) \cdot x_k$.

This might reveal something to us about who's causing the trouble, over and above which variables are causing the trouble. "A residual plot for improving poorly balanced propensity score matches" is a very plausible title for a _J Comp Graph Stats_ paper.)

### Likely Red Herrings; low priority ###

-Wrong link function. (If `x`es were drawn from a continuous dist'n this couldn't possibly cause trouble.  Unclear whether it can with discrete `x`es.)
- PS model needs more interaction terms. (why?)


## Possible PS Matching Issues ##
I.e., fitted score is good enough but haven't matched closely enough on it

### Max of PS Differences Too Large ###

Test matched differences on the fitted PS for distinguishability from 0.  Apply Bonferroni. 

### Mean PS Difference Too Large ###

- Calculate l2 distance on fitted PS, compare its square to null expectation of squared l2 distance on fitted PS
- 

