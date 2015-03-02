
# Pairwise comparisons for all sites.
#' \code{hill_func_parti_pairwise} to calculate pairwise functional gamma, alpha, and beta diversity for communities, as
#'  well as site similarity. It is based on \code{\link{hill_func_parti}}. If comm has >2 sites, this function will give results for all
#'  pairwise comparisons.
#'
#' @author Daijiang Li
#'
#' @param comm data frame of vegtation data. Sites as rows, species as columns.
#' @param traits data frame of species functional traits data. Species as rows, traits as columns.
#' It can include both continuous and categorical data.
#' @param traits_as_is if FALSE (default) traits data frame will be transformed into a distance
#' matrix using `FD::gowdis(traits)`. Otherwise, will use as is (i.e. traits is a symmetric distance matrix).
#' @param q hill number, q = 0 (default) to get species richness, q = 1 to get shannon entropy, q = 2 will give inverse Simpson.
#' @param rel_then_pool default is TRUE. Abundance of species are first changed to relative abundance within sites,
#'  then pooled into one assemblage. If FALSE, sites are pooled first, then change abundance of species
#'  to relative abundance.
#' @export
#' @return a data frame with results for all pairwise comparisons.
#' @seealso \code{\link{hill_func_parti}}
#' @examples
#' library(FD); data(dummy)
#' hill_func_parti_pairwise(comm = dummy$abun, traits = dummy$trait, q = 0)
#' hill_func_parti_pairwise(comm = dummy$abun, traits = dummy$trait, q = 1)
#' hill_func_parti_pairwise(comm = dummy$abun, traits = dummy$trait, q = 0.9999)
#' hill_func_parti_pairwise(comm = dummy$abun, traits = dummy$trait, q = 2)
#' hill_func_parti_pairwise(comm = dummy$abun, traits = dummy$trait, q = 3)
#'
#'
hill_func_parti_pairwise = function(comm, traits, traits_as_is = FALSE,
                                    q = 0, rel_then_pool = TRUE){
  nsite = nrow(comm)
  site.comp = as.matrix(expand.grid(1:nsite, 1:nsite))
  adply(site.comp, 1, function(x){
    data.frame(site1 = row.names(comm)[x[1]],
               site2 = row.names(comm)[x[2]],
               hill_func_parti(comm = comm[c(x[1], x[2]), ], traits = traits,
                               traits_as_is = traits_as_is, q = q,
                               rel_then_pool = rel_then_pool))
  })
}
