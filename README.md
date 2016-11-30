## Taxonomic, functional, and phylogenetic diversity and similarity through Hill Numbers

This package contains R functions to calculate diversity and similarity measures based on Chao, Chiu and Jost 2014 and Chiu & Chao 2014.

>[Chao, Anne, Chun-Huo Chiu, and Lou Jost. “Unifying Species Diversity, Phylogenetic Diversity, Functional Diversity, and Related Similarity and Differentiation Measures Through Hill Numbers.” Annual Review of Ecology, Evolution, and Systematics 45, no. 1 (2014): 297–324. doi:10.1146/annurev-ecolsys-120213-091540.](http://dx.doi.org/10.1146/annurev-ecolsys-120213-091540)
>

>[Chiu, Chun-Huo, and Anne Chao. “Distance-Based Functional Diversity Measures and Their Decomposition: A Framework Based on Hill Numbers.” PLoS ONE 9, no. 7 (July 7, 2014): e100014. doi:10.1371/journal.pone.0100014.](http://dx.doi.org/10.1371/journal.pone.0100014)


At this moment, only taxonomic and functional diversity and similarity were included.

## Installation
To install this package, just run:

    library(devtools)
    install_github("daijiang/hillR")
    
## Warning
~~This package is just initiated WITHOUT test! Be careful.~~

update: I just compared results with [vegetarian](http://cran.r-project.org/web/packages/vegetarian/index.html) R package and they are the same, for taxonomic diversity. I also have tested functional part and everything looks OK.

## Licenses

Licensed under the [MIT license](LICENSE). ([More information here](http://en.wikipedia.org/wiki/MIT_License).)

## Citation

Please cite this package if you use it.

Daijiang Li. hillR: taxonomic, functional, and phylogenetic diversity and similarity through Hill Numbers. R package.


## Note:
You need `library(dplyr)` to make these functions work. I need to put it as a dependency package later.
