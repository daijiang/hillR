# hillR 0.5.0

- `hill_phylo` set of functions are now much faster after getting rid of `ade4::newick2phylog()`.
- Replaced all `.` in arguments with `_` to standardize styles.
- `hill_phylo_parti_pairwise` now calculate the node/tips by site matrix first to speed things up.
- Added progress bar to all `hill_xxx_parti_pairwise` functions with `.progress` argument, which can be set to `FALSE` to turn it off.

# hillR 0.3.0

- Added a `NEWS.md` file to track changes to the package.
- `hill_func` and `hill_func_parti` no longer standardize distance matrix to have max value of 1 by default.
  + set `stand_dij = TRUE` if want to do so.
- `hill_func_parti` function now reports site **similarity** instead of dissimilarity, to be consistent with the `hill_taxa_parti` functions. Make sure to change your interpretations.
- Added `hill_phylo` functions.
