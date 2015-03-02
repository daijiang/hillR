# Pairwise comparisons for all sites.
#' \code{hill_taxa_parti_pairwise} to calculate pairwise taxonomic gamma, alpha, and beta diversity for communities, as
#'  well as site similarity. It is based on \code{\link{hill_taxa_parti}}. If comm has >2 sites, this function will give results for all
#'  pairwise comparisons.
#'
#' @author Daijiang Li
#'
#' @param comm data frame of vegtation data. Sites as rows, species as columns.
#' @param q hill number, q = 0 (default) to get species richness, q = 1 to get shannon entropy, q = 2 will give inverse Simpson.
#' @param rel_then_pool default is TRUE. Abundance of species are first changed to relative abundance within sites,
#'  then pooled into one assemblage. If FALSE, sites are pooled first, then change abundance of species
#'  to relative abundance.
#' @export
#' @return a data frame with results for all pairwise comparisons.
#' @seealso \code{\link{hill_taxa_parti}}
#' @examples
#' library(FD); data(dummy)
#' hill_taxa_parti_pairwise(comm = dummy$abun, q = 0)
#' hill_taxa_parti_pairwise(comm = dummy$abun, q = 1)
#' hill_taxa_parti_pairwise(comm = dummy$abun, q = 0.9999999)
#' hill_taxa_parti_pairwise(comm = dummy$abun, q = 0.9999999, rel_then_pool = F)
#' hill_taxa_parti_pairwise(comm = dummy$abun, q = 1, rel_then_pool = F)
#' hill_taxa_parti_pairwise(comm = dummy$abun, q = 2)
#' hill_taxa_parti_pairwise(comm = dummy$abun, q = 3)
#'
#'
hill_taxa_parti_pairwise = function(comm, q = 0, rel_then_pool = TRUE){
  nsite = nrow(comm)
  site.comp = as.matrix(expand.grid(1:nsite, 1:nsite))
  adply(site.comp, 1, function(x){
    data.frame(site1 = row.names(comm)[x[1]],
               site2 = row.names(comm)[x[2]],
               hill_taxa_parti(comm = comm[c(x[1], x[2]), ], q = q,
                               rel_then_pool = rel_then_pool))
  })
}
