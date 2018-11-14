#' hillR: Diversity Through Hill Numbers
#'
#' The R package `hillR` implements the framework proposed by Chao, et al. 2014 and makes it easy to calculate taxonomic, functional, and phylogenetic diversity of ecological communities as Hill numbers. For each facet of diversity, `hillR` has three functions. The first set of functions (\code{hill_taxa}, \code{hill_func}, and \code{hill_phylo}) calculates alpha diversity of each site. The second set of functions (\code{hill_taxa_parti}, \code{hill_func_parti}, and \code{hill_phylo_parti}) calculates diversity across all sites. The third set of functions (\code{hill_taxa_parti_pairwise}, \code{hill_func_parti_pairwise}, and \code{hill_phylo_parti_pairwise}) calculates all possible pairwise diversity across all sites. Users can set the argument _q_ to control the weight of species abundance.
#'
#' Users may be interested in other similar packages such as \href{ https://CRAN.R-project.org/package=vegetarian}{vegetarian} and \href{ https://CRAN.R-project.org/package=iNEXT}{iNEXT}.
#'
#' @section Taxonomic Hill Numbers:
#' \code{\link{hill_taxa}}, \code{\link{hill_taxa_parti}}, \code{\link{hill_taxa_parti_pairwise}}
#'
#' @section Functional Hill Numbers:
#' \code{\link{hill_func}}, \code{\link{hill_func_parti}}, \code{\link{hill_func_parti_pairwise}}
#'
#' @section Phylogenetic Hill Numbers:
#' \code{\link{hill_phylo}}, \code{\link{hill_phylo_parti}}, \code{\link{hill_phylo_parti_pairwise}}
#'
#' @references  Chao, Anne, Chun-Huo Chiu, and Lou Jost. Unifying Species Diversity, Phylogenetic Diversity, Functional Diversity, and Related Similarity and Differentiation Measures Through Hill Numbers. Annual Review of Ecology, Evolution, and Systematics 45, no. 1 (2014): 297–324. <doi:10.1146/annurev-ecolsys-120213-091540>.
#'
#'
#' @docType package
#' @name hillR
NULL
