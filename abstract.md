## Diminishing caliper matching for propensity and other index scores


To estimate intervention effects without the benefit of random
assignment, an often useful beginning is to pair intervention group
members to ostensibly similar counterparts receiving a control
condition.  In practice exact matches are rare, particularly if there
are many measured covariates. Instead, matches may be made within
_calipers_ (Althauser & Rubin, 1970) of a unidimensional index (such
as Rosenbaum & Rubin's [1983] propensity score).

This talk presents a new way to determine calipers.  Compatible with
common index model specifications, its widths diminish as _n_
increases, toward an asymptote of 0.  If the index model is
consistently estimated, then matched contrast-based impact estimates
will be consistent as well, provided matches are made within these
diminishing calipers. This result assumes no hidden bias, an
untestable condition, alongside of additional conditions that can be
enforced.
