## Propensity score calipers and the overlap condition

Propensity scores (Rosenbaum and Rubin, 1983) are used widely to
address measured confounding in quasiexperiments. They also arise in
connection with the antecedent question of whether non-equivalent
treatment and control groups are suitable for comparison at all, with or
without covariate adjustments.

"Common support," the assumption that propensity scores are
bounded away from 1, is so named because it means that the support of
the treatment group in covariate spaces is contained within that of the
control group.  This is less simple to check than is often
supposed: even if treatment and control groups' values of the
_true_ propensity score overlap, when arranged in order of _estimated_
propensity scores they may appear not to. The naive
method of discarding those members of the treatment group whose
estimated propensity scores fall above all the controls', and those
members of the control group whose propensity scores fall below all
those estimated within the treatment group, is needlessly wasteful of
sample size.

It is possible to address common support by restricting the range of the
estimated propensity score within which comparisons are permitted, but
this requires careful determination of the tolerance enforced for
matching discrepancies; available heuristics and guidelines attend only
to some of the issues that must be considered. I present a new formula
for determining caliper widths for matching based on propensity scores,
and other constructed index functions. The method is compatible with
conventional means of propensity score estimation, although it relaxes
some of the more tenuous of the conventional assumptions, in particular
permitting the dimension of the parameter to grow with n.
