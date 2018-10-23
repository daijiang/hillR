---
title: 'hillR: taxonomic, functional, and phylogenetic diversity and similarity through Hill Numbers'
tags:
  - R
  - Hill numbers
  - biodiversity
  - phylogenetic diversity
  - functional diversity
authors:
  - name: Daijiang Li
    orcid: 0000-0002-0925-3421
    affiliation: 1
affiliations:
 - name: Department of Wildlife Ecology and Conservation, University of Florida, Gainesville, FL 32611
   index: 1
date: 12 September 2018
bibliography: paper.bib
---

# Summary

A unified framework to calculate biodiversity of ecological communities through Hill numbers [@hill1973diversity; @jost2006entropy] was recently proposed by @chao2014unifying. This framework can be applied to three facets of biodiversity: taxonomic, functional, and phylogenetic diversity. These three facets of biodiversity can be expressed as Hill numbers with the same units (effective number of species), facilitating direct comparisons between them [@chao2014unifying]. In addition, this framework can account for species abundance by changing the parameter _q_, which can be any non-negative number (e.g., 0, 0.99, 2). As _q_ increases, the diversity values become more sensitive to common species. When _q_ = 0, species abundance is ignored; _q_ = 1, all species are weighted by their abundance equally (i.e., Shannon's diversity for taxonomic diversity); _q_ = 2, common species get more weight than rare species (i.e., inverse of Simpson diversity for taxonomic diversity). Furthermore, both alpha and beta diversity can be calculated with this framework. Given the above advantages, this framework has been adopted by an increasing number of researchers.

The R package `hillR` implements the framework proposed by @chao2014unifying and makes it easy to calculate taxonomic, functional, and phylogenetic diversity of ecological communities as Hill numbers. For each facet of diversity, `hillR` has three functions. The first set of functions (`hill_taxa()`, `hill_func()`, and `hill_phylo()`) calculates alpha diversity of each site. The second set of functions (`hill_taxa_parti()`, `hill_func_parti()`, and `hill_phylo_parti()`) calculates diversity across all sites. The third set of functions (`hill_taxa_parti_pairwise()`, `hill_func_parti_pairwise()`, and `hill_phylo_parti_pairwise()`) calculates all possible pairwise diversity across all sites. Users can set the argument _q_ to control the weight of species abundance. For usage examples of these functions, see the [README file on Github](https://github.com/daijiang/hillR/blob/master/README.md).

The R package `hillR` is available on [Github](https://github.com/daijiang/hillR), where issues can be opened.

# Reference
