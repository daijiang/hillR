# hillR 0.4.0

# hillR 0.3.0

- Added a `NEWS.md` file to track changes to the package.
- `hill_func` and `hill_func_parti` no longer standardize distance matrix to have max value of 1 by default.
  + set `stand_dij = TRUE` if want to do so.
- `hill_func_parti` function now reports site **similarity** instead of dissimilarity, to be consistent with the `hill_taxa_parti` functions. Make sure to change your interpretations.
- Added `hill_phylo` functions.
