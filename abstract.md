## Diminishing caliper matching for propensity and other index scores


To estimate intervention effects without the benefit of random
assignment, an often useful beginning is to pair intervention group
members to ostensibly similar counterparts receiving a control
condition.  In practice exact matches are rare, particularly if there
are many measured covariates. Instead, matches may be made within
_calipers_ (Althauser & Rubin, 1970) of a unidimensional index.
Modern indices arise by modeling specific aspects of the data.  The
most widely used matching indices are propensity scores (Rosenbaum &
Rubin, 1983), followed by risk or prognostic scores (Miettinen, 1976;
Hansen, 2008).

Adjudicating how close is close enough for matching is the murkiest
aspect of the undertaking.  Heuristics in wide use today pre-date the
use of model-based matching indices, fail to adapt to the size of the
model and sample, and lack theoretical support.  In some cases
these heuristics allow pairings of demonstrably dissimilar subjects;
in others they declare wide swaths of the sample to be unmatchable,
needlessly wasting data.

This talk presents a new way to determine calipers.  Compatible with
common index model specifications, its widths diminish as _n_
increases, toward an asymptote of 0.  If the index model is
consistently estimated, then matched contrast-based impact estimates
will be consistent as well, provided matches are made within these
diminishing calipers. This result assumes no hidden bias, an
untestable condition, alongside of additional conditions that can be
enforced. In particular, it restricts growth of the index parameter's
dimension relative to _n_, to a rate intermediate to those required for
ordinary M-estimates to be consistent or root-_n_ consistent (He &
Shao, 2000).