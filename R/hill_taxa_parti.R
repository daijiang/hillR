#' Decompostion of Taxonomic diversity through Hill Numbers
#'
#' Calculate taxonomic gamma, alpha, and beta diversity for communities, as
#' well as site similarity. If comm has 2 sites, this function gives pair comparison.
#' If comm has >2 sites, gamma diversity is the diversity of the pooled assemblage,
#' alpha is the average diversity across all site, beta is across all communities.
#'
#' @author Daijiang Li
#'
#' @param comm data frame of vegtation data. Sites as rows, species as columns.
#' @param q hill number, q = 0 (default) to get species richness, q = 1 to get shannon entropy, q = 2 will give inverse Simpson.
#' @param base default is exp(1), the base of log.
#' @param rel_then_pool default is TRUE. Abundance of species are first changed to relative abundance within sites,
#'  then pooled into one assemblage. If FALSE, sites are pooled first, then change abundance of species
#'  to relative abundance.
#'  @param show.warning whether to print warning, default is TRUE
#' @export
#' @return a data frame with one row, including these columns: q, gamma diversity, alpha diveristy,
#' beta diversity, MacArthur's homogeneity measure, local similarity (species overlap),
#' and region similarity (species overlap).
#' See Chao, Chiu and Jost 2014 Table 2 for more information.
#' @seealso \code{\link{hill_taxa_parti}}
#' @examples
#' dummy = FD::dummy
#' hill_taxa_parti(comm = dummy$abun, q = 0)
#' hill_taxa_parti(comm = dummy$abun, q = 1)
#' hill_taxa_parti(comm = dummy$abun, q = 0.9999999)
#' hill_taxa_parti(comm = dummy$abun, q = 0.9999999, rel_then_pool = FALSE)
#' hill_taxa_parti(comm = dummy$abun, q = 1, rel_then_pool = FALSE)
#' hill_taxa_parti(comm = dummy$abun, q = 2)
#' hill_taxa_parti(comm = dummy$abun, q = 3)
#'
#'
hill_taxa_parti = function(comm, q = 0, base = exp(1),
                           rel_then_pool = TRUE, show.warning = TRUE){
  if (any(comm < 0)) stop("Negative value in comm data")
  if(any(colSums(comm) == 0) & show.warning) warning("Some species in comm data were not observed in any site,\n delete them...")
  comm = comm[, colSums(comm) != 0]
  N = nrow(comm)
  S = ncol(comm)
  comm = as.matrix(comm)
  if(rel_then_pool){
    comm_gamma = colSums(sweep(comm, 1, rowSums(comm, na.rm = TRUE), "/"))/ N
    # relative abun
  } else {
    comm_gamma = colSums(comm)/sum(comm)
  }
  if(!all.equal(sum(comm_gamma), 1)) stop("Accumlative relative abundance should be 1")

  if(rel_then_pool){
    comm_alpha = sweep(comm, 1, rowSums(comm, na.rm = TRUE), "/") # relative abun
  } else {
    comm_alpha = comm
  }

  TD_q_gamma = hill_taxa(comm_gamma, q = q)

  # TD_q_alpha
  if(q == 0){
    TD_q_alpha = sum(hill_taxa(comm_alpha, q = 0))/N # mean of sp richness
  } else {
    if(q == 1){
      TD_q_alpha = exp(-1 * sum((comm_alpha/sum(comm_alpha)) *
                                  log((comm_alpha/sum(comm_alpha)), base), na.rm = T)) / N
    } else{
      TD_q_alpha = (1/N) * (sum((comm_alpha/sum(comm_alpha))^q)^(1/(1-q)))
    }
  }

  TD_q_beta = TD_q_gamma / TD_q_alpha

  if(q == 1){ # why is negative and abs()>1 ??
    if(rel_then_pool){
      local_taxa_overlap = (log(N, base) - log(TD_q_gamma) + log(TD_q_alpha))/log(N, base)
    } else{
      local_taxa_overlap = (log(TD_q_alpha, base) - log(TD_q_gamma, base) -
                              sum((rowSums(comm_alpha)/sum(comm_alpha))*
                                    log(rowSums(comm_alpha)/sum(comm_alpha), base)))/log(N, base)}
  } else {
    local_taxa_overlap = (N^(1*(1-q)) - TD_q_beta^(1-q)) / (N^(1*(1-q)) - 1)
  }

  if(q == 1){
    region_taxa_overlap = local_taxa_overlap
  } else {
    region_taxa_overlap = ((1/TD_q_beta)^(1-q) - (1/N)^(1*(1-q))) / (1 - (1/N)^(1*(1-q)))
  }

  return(data.frame(q = q,
                    TD_gamma = TD_q_gamma,
                    TD_alpha = TD_q_alpha,
                    TD_beta = TD_q_beta,
                    M_homog = 1 / TD_q_beta,
                    local_similarity = local_taxa_overlap,
                    region_similarity = region_taxa_overlap))
}
