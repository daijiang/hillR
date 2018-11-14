
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![DOI](http://joss.theoj.org/papers/10.21105/joss.01041/status.svg)](https://doi.org/10.21105/joss.01041)
[![Build
Status](https://travis-ci.org/daijiang/hillR.svg?branch=master)](https://travis-ci.org/daijiang/hillR)
[![Coverage
status](https://codecov.io/gh/daijiang/hillR/branch/master/graph/badge.svg)](https://codecov.io/github/daijiang/hillR?branch=master)

# hillR

This package contains R functions to calculate taxonomic, functional,
and phylogenetic diversity and site similarity through Hill Numbers. The
underlying methods are based on Chao, Chiu and Jost 2014 and Chiu & Chao
2014.

# Installation

To install this package, run the following code:

    devtools::install_github("daijiang/hillR")

# Examples

``` r
dummy = FD::dummy
comm = dummy$abun
traits = dummy$trait
tree = ape::rtree(n = ncol(comm), tip.label = paste0("sp", 1:ncol(comm)))
library(hillR)
```

## Calculate taxonomic, functional, and phylogenetic diversity of each site

``` r
hill_taxa(comm, q = 0) # taxonomic alpha diversity
##  com1  com2  com3  com4  com5  com6  com7  com8  com9 com10 
##     4     3     3     2     3     5     3     4     5     4

hill_func(comm, traits, q = 0) # functional alpha diversity
##           com1      com2      com3      com4      com5      com6      com7
## Q    0.4016663 0.1922618 0.2780442 0.1146261 0.3816159  0.404177 0.2934143
## FDis 0.3481687 0.1670560 0.2375808 0.1146261 0.3211366  0.330233 0.2532751
## D_q  4.0974923 3.6518111 3.2454591 3.0237158 3.0655375  5.233241 3.1470056
## MD_q 1.6458245 0.7021037 0.9023810 0.3465969 1.1698580  2.115156 0.9233765
## FD_q 6.7437533 2.5639502 2.9286406 1.0480104 3.5862436 11.069121 2.9058708
##           com8       com9     com10
## Q    0.3343662  0.4156546 0.3844765
## FDis 0.2877931  0.3421687 0.3503927
## D_q  4.3998540  5.2114653 4.1694097
## MD_q 1.4711625  2.1661695 1.6030400
## FD_q 6.4729004 11.2889174 6.6837303

hill_phylo(comm, tree, q = 0) # phylogenetic alpha diversity
##     com1     com2     com3     com4     com5     com6     com7     com8 
## 6.587145 5.025573 3.290982 4.113949 4.878235 6.188022 5.099523 5.386479 
##     com9    com10 
## 7.069928 5.751388
```

## Calculate taxonomic, functional, and phylogenetic diversity across multiple sites

``` r
hill_taxa_parti(comm, q = 0) # taxonomic diversity across all sites
##   q TD_gamma TD_alpha  TD_beta M_homog local_similarity region_similarity
## 1 0        8      3.6 2.222222    0.45        0.8641975         0.3888889

hill_func_parti(comm, traits, q = 0) # functional diversity across all sites
##   q raoQ_gamma FD_gamma FD_alpha FD_beta local_similarity
## 1 0  0.4529152 29.66099 14.15941 2.09479        0.9889415
##   region_similarity
## 1         0.4720957

hill_phylo_parti(comm, tree, q = 0) # phylogenetic diversity across all sites
##   q PD_gamma PD_alpha PD_beta local_similarity region_similarity
## 1 0 8.714781 5.339122 1.63225          0.92975         0.5696126
```

## Calculate pairwise taxonomic, functional, and phylogenetic diversity

``` r
hill_taxa_parti_pairwise(comm, q = 0, show.warning = F) # pairwise taxonomic diversity
## # A tibble: 45 x 8
##        q site1 site2 TD_gamma TD_alpha TD_beta local_similarity
##    <dbl> <fct> <fct>    <dbl>    <dbl>   <dbl>            <dbl>
##  1     0 com1  com2         6      3.5    1.71            0.286
##  2     0 com1  com3         5      3.5    1.43            0.571
##  3     0 com2  com3         5      3      1.67            0.333
##  4     0 com1  com4         5      3      1.67            0.333
##  5     0 com2  com4         5      2.5    2               0    
##  6     0 com3  com4         4      2.5    1.6             0.400
##  7     0 com1  com5         6      3.5    1.71            0.286
##  8     0 com2  com5         4      3      1.33            0.667
##  9     0 com3  com5         6      3      2               0    
## 10     0 com4  com5         4      2.5    1.6             0.400
## # ... with 35 more rows, and 1 more variable: region_similarity <dbl>

hill_func_parti_pairwise(comm, traits, q = 0, show.warning = F) # pairwise functional diversity
## # A tibble: 45 x 8
##        q site1 site2 FD_gamma FD_alpha FD_beta local_similarity
##    <dbl> <fct> <fct>    <dbl>    <dbl>   <dbl>            <dbl>
##  1     0 com1  com2     15.9     10.3     1.54            0.821
##  2     0 com1  com3     10.7      7.96    1.35            0.883
##  3     0 com2  com3     11.1      7.03    1.58            0.807
##  4     0 com1  com4     11.6      7.90    1.47            0.843
##  5     0 com2  com4     11.7      6.85    1.70            0.765
##  6     0 com3  com4      6.60     4.45    1.48            0.839
##  7     0 com1  com5     17.3     11.3     1.54            0.821
##  8     0 com2  com5      7.86     5.92    1.33            0.891
##  9     0 com3  com5     16.2      9.72    1.66            0.780
## 10     0 com4  com5      8.00     5.32    1.50            0.832
## # ... with 35 more rows, and 1 more variable: region_similarity <dbl>

hill_phylo_parti_pairwise(comm, tree, q = 0, show.warning = F) # pairwise phylogenetic diversity
## # A tibble: 45 x 8
##        q site1 site2 PD_gamma PD_alpha PD_beta local_similarity
##    <dbl> <fct> <fct>    <dbl>    <dbl>   <dbl>            <dbl>
##  1     0 com1  com2      7.63     5.81    1.31            0.685
##  2     0 com1  com3      6.84     4.94    1.39            0.615
##  3     0 com2  com3      6.44     4.16    1.55            0.450
##  4     0 com1  com4      7.54     5.35    1.41            0.590
##  5     0 com2  com4      6.76     4.57    1.48            0.521
##  6     0 com3  com4      5.85     3.70    1.58            0.419
##  7     0 com1  com5      8.33     5.73    1.45            0.546
##  8     0 com2  com5      5.98     4.95    1.21            0.792
##  9     0 com3  com5      7.40     4.08    1.81            0.189
## 10     0 com4  com5      5.66     4.50    1.26            0.742
## # ... with 35 more rows, and 1 more variable: region_similarity <dbl>
```

# Licenses

Licensed under the [MIT license](LICENSE). ([More information
here](http://en.wikipedia.org/wiki/MIT_License).)

# Citation

Please cite this package if you use it. The citation information can be
obtained by running `citation('hillR')` in R.

> Li, (2018). hillR: taxonomic, functional, and phylogenetic diversity
> and similarity through Hill Numbers. Journal of Open Source Software,
> 3(31), 1041. <https://doi.org/10.21105/joss.01041>

``` bibtex
@Article{,
    title = {hillR: taxonomic, functional, and phylogenetic diversity and similarity through Hill Numbers},
    author = {Daijiang Li},
    journal = {Journal of Open Source Software},
    year = {2018},
    volume = {3},
    number = {31},
    pages = {1041},
    url = {https://doi.org/10.21105/joss.01041},
 }
```

# Reference

  - [Chao, Anne, Chun-Huo Chiu, and Lou Jost. “Unifying Species
    Diversity, Phylogenetic Diversity, Functional Diversity, and Related
    Similarity and Differentiation Measures Through Hill Numbers.”
    Annual Review of Ecology, Evolution, and Systematics 45, no. 1
    (2014): 297–324.
    doi:10.1146/annurev-ecolsys-120213-091540.](https://doi.org/10.1146/annurev-ecolsys-120213-091540)
  - [Chiu, Chun-Huo, and Anne Chao. “Distance-Based Functional Diversity
    Measures and Their Decomposition: A Framework Based on Hill
    Numbers.” PLoS ONE 9, no. 7 (July 7, 2014): e100014.
    doi:10.1371/journal.pone.0100014.](https://doi.org/10.1371/journal.pone.0100014)

# Contributing

Contributions are welcome. You can provide comments and feedback or ask
questions by filing an issue on Github
[here](https://github.com/daijiang/hillR/issues) or making pull
requests.
