
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![DOI](http://joss.theoj.org/papers/10.21105/joss.01041/status.svg)](https://doi.org/10.21105/joss.01041)
[![Build
Status](https://travis-ci.org/daijiang/hillR.svg?branch=master)](https://travis-ci.org/daijiang/hillR)
[![Coverage
status](https://codecov.io/gh/daijiang/hillR/branch/master/graph/badge.svg)](https://codecov.io/github/daijiang/hillR?branch=master)
[![CRAN
status](https://www.r-pkg.org/badges/version/hillR)](https://cran.r-project.org/package=hillR)
[![CRAN RStudio mirror
downloads](http://cranlogs.r-pkg.org/badges/hillR)](http://www.r-pkg.org/pkg/hillR)
[![CRAN RStudio mirror
downloads](http://cranlogs.r-pkg.org/badges/grand-total/hillR?color=green)](http://www.r-pkg.org/pkg/hillR)

# hillR

This package contains R functions to calculate taxonomic, functional,
and phylogenetic diversity and site similarity through Hill Numbers. The
underlying methods are based on Chao, Chiu and Jost 2014 and Chiu & Chao
2014.

# Installation

To install this package, run the following code:

``` r
install.packages("hillR")
# or install from Github
devtools::install_github("daijiang/hillR")
```

# Examples

``` r
dummy = FD::dummy
comm = dummy$abun
traits = dummy$trait
set.seed(123)
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
## 5.430079 4.684280 4.461773 2.551395 5.830078 6.088533 4.763594 6.046474 
##     com9    com10 
## 6.262164 5.080340
```

## Calculate taxonomic, functional, and phylogenetic diversity across multiple sites

``` r
hill_taxa_parti(comm, q = 0) # taxonomic diversity across all sites
##   q TD_gamma TD_alpha  TD_beta M_homog local_similarity region_similarity
## 1 0        8      3.6 2.222222    0.45        0.8641975         0.3888889

hill_func_parti(comm, traits, q = 0) # functional diversity across all sites
##   q raoQ_gamma FD_gamma FD_alpha FD_beta local_similarity region_similarity
## 1 0  0.4529152 29.66099 14.15941 2.09479        0.9889415         0.4720957

hill_phylo_parti(comm, tree, q = 0) # phylogenetic diversity across all sites
##   q PD_gamma PD_alpha  PD_beta local_similarity region_similarity
## 1 0 8.292885 5.119871 1.619745        0.9311395          0.574868
```

## Calculate pairwise taxonomic, functional, and phylogenetic diversity

``` r
# pairwise taxonomic diversity
hill_taxa_parti_pairwise(comm, q = 0, show_warning = FALSE, .progress = FALSE) 
## # A tibble: 45 x 8
##        q site1 site2 TD_gamma TD_alpha TD_beta local_similarity region_similari…
##    <dbl> <chr> <chr>    <dbl>    <dbl>   <dbl>            <dbl>            <dbl>
##  1     0 com1  com2         6      3.5    1.71            0.286            0.167
##  2     0 com1  com3         5      3.5    1.43            0.571            0.400
##  3     0 com2  com3         5      3      1.67            0.333            0.200
##  4     0 com1  com4         5      3      1.67            0.333            0.200
##  5     0 com2  com4         5      2.5    2               0                0    
##  6     0 com3  com4         4      2.5    1.6             0.400            0.25 
##  7     0 com1  com5         6      3.5    1.71            0.286            0.167
##  8     0 com2  com5         4      3      1.33            0.667            0.5  
##  9     0 com3  com5         6      3      2               0                0    
## 10     0 com4  com5         4      2.5    1.6             0.400            0.25 
## # … with 35 more rows
# pairwise functional diversity
hill_func_parti_pairwise(comm, traits, q = 0, show_warning = FALSE, .progress = FALSE) 
## # A tibble: 45 x 8
##        q site1 site2 FD_gamma FD_alpha FD_beta local_similarity region_similari…
##    <dbl> <chr> <chr>    <dbl>    <dbl>   <dbl>            <dbl>            <dbl>
##  1     0 com1  com2     15.9     10.3     1.54            0.821            0.534
##  2     0 com1  com3     10.7      7.96    1.35            0.883            0.655
##  3     0 com2  com3     11.1      7.03    1.58            0.807            0.511
##  4     0 com1  com4     11.6      7.90    1.47            0.843            0.573
##  5     0 com2  com4     11.7      6.85    1.70            0.765            0.449
##  6     0 com3  com4      6.60     4.45    1.48            0.839            0.566
##  7     0 com1  com5     17.3     11.3     1.54            0.821            0.534
##  8     0 com2  com5      7.86     5.92    1.33            0.891            0.671
##  9     0 com3  com5     16.2      9.72    1.66            0.780            0.469
## 10     0 com4  com5      8.00     5.32    1.50            0.832            0.554
## # … with 35 more rows
# pairwise phylogenetic diversity
hill_phylo_parti_pairwise(comm, tree, q = 0, show_warning = FALSE, .progress = FALSE) 
## # A tibble: 45 x 8
##        q site1 site2 PD_gamma PD_alpha PD_beta local_similarity region_similari…
##    <dbl> <chr> <chr>    <dbl>    <dbl>   <dbl>            <dbl>            <dbl>
##  1     0 com1  com2      6.79     5.06    1.34           0.657            0.489 
##  2     0 com1  com3      6.14     4.95    1.24           0.759            0.611 
##  3     0 com2  com3      6.75     4.57    1.48           0.523            0.355 
##  4     0 com1  com4      6.38     3.99    1.60           0.400            0.250 
##  5     0 com2  com4      7.13     3.62    1.97           0.0284           0.0144
##  6     0 com3  com4      5.42     3.51    1.54           0.455            0.295 
##  7     0 com1  com5      7.04     5.63    1.25           0.750            0.599 
##  8     0 com2  com5      6.54     5.26    1.24           0.756            0.608 
##  9     0 com3  com5      7.71     5.15    1.50           0.502            0.335 
## 10     0 com4  com5      6.42     4.19    1.53           0.467            0.305 
## # … with 35 more rows
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
