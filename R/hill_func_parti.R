# Decompostion of functional diversity through Hill Numbers
#' \code{hill_func_parti} to calculate functional gamma, alpha, and beta diversity for communities, as
#'  well as site similarity.
#'
#' @author Daijiang Li
#'
#' @param comm data frame of vegtation data. Sites as rows, species as columns.
#' @param traits data frame of species functional traits data. Species as rows, traits as columns.
#' It can include both continuous and categorical data.
#' @param traits_as_is if FALSE (default) traits data frame will be transformed into a distance
#' matrix using `FD::gowdis(traits)`. Otherwise, will use as is (i.e. traits is a symmetric distance matrix).
#' @param q hill number, q = 0 (default) to get species richness, q = 1 to get shannon entropy, q = 2 will give inverse Simpson.
#' @param base default is exp(1), the base of log.
#' @param rel_then_pool default is TRUE. Abundance of species are first changed to relative abundance within sites,
#'  then pooled into one assemblage. If FALSE, sites are pooled first, then change abundance of species
#'  to relative abundance.
#' @export
#' @return  a data frame with one row, including these columns: q, RaoQ of pooled assemblage,
#' gamma diversity, alpha diveristy, beta diversity, local species overlap, and region species
#' overlap. See Chiu and Chao 2014 Table 3 for more information.
#' @seealso \code{\link{hill_taxa_parti}}, \code{\link{hill_func}}
#'
#' @examples
#' library(FD); data(dummy)
#' hill_func_parti(comm = dummy$abun, traits = dummy$trait, q = 0)
#' hill_func_parti(comm = dummy$abun, traits = dummy$trait, q = 1)
#' hill_func_parti(comm = dummy$abun, traits = dummy$trait, q = 0.9999)
#' hill_func_parti(comm = dummy$abun, traits = dummy$trait, q = 2)
#' hill_func_parti(comm = dummy$abun, traits = dummy$trait, q = 3)
#'
#'
hill_func_parti = function(comm, traits, traits_as_is = FALSE, q = 0,
                           base = exp(1), checkdata=TRUE,
                           rel_then_pool = TRUE, ...){
  if (checkdata) {
    if (any(comm < 0))
      stop("Negative value in comm data")
    if (is.null(rownames(traits))) {
      stop("\n Traits have no row names\n")
    }
    if (is.null(colnames(comm))) {
      stop("\n Comm data have no col names\n")
    }
  }

  if(any(colSums(comm) == 0)) warning("Some species in comm data were not observed in any site,\n delete them...")
  comm = comm[, colSums(comm) != 0]
  comm = as.matrix(comm)
  N = nrow(comm)
  S = ncol(comm)

  if(any(!colnames(comm) %in% rownames(traits))){
    warning("\n There are species from community data that are not on traits matrix\nDelete these species from comm data...\n")
    comm = comm[, colnames(comm) %in% rownames(traits)]
  }

  if(any(!rownames(traits) %in% colnames(comm))){
    warning("\n There are species from trait data that are not on comm matrix\nDelete these species from trait data...\n")
    traits = traits[rownames(traits) %in% colnames(comm), ]
  }

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

  if(traits_as_is){
    traits.dist = traits
  } else {
    traits.dist = as.matrix(gowdis(x=traits, ...))
  }

  dij = traits.dist^2 # trait distance matrix
  Q_gamma = as.vector(comm_gamma %*% dij %*% matrix(comm_gamma, ncol = 1))

  ## FD_q_gamma
  if(q == 1){
    FD_q_gamma = exp(-1 * sum(dij *
                                (outer(comm_gamma, comm_gamma, FUN = "*")/Q_gamma) *
                                log(outer(comm_gamma, comm_gamma, FUN = "*")/Q_gamma)))
    # Chiu & Chao 2014 p.7, equ 6b
  } else{ #q != 0 or 1
    FD_q_gamma = sum(dij *
                       ((outer(comm_gamma, comm_gamma, FUN = "*")/Q_gamma)^q))^(1/(1-q))
    # Chiu & Chao 2014 p.7, equ 6a
  }

  ## FD_q_alpha
  if(q == 1){
    x = (outer(comm_alpha, comm_alpha, FUN = "*"))/(Q_gamma * (sum(comm_alpha)^2))
    x[x==0] = NA
    xx = x * log(x, base)
    # then * dij
    for(k in 1:N){
      for(j in 1:S){
        xx[,,k,j] = sweep(xx[,,k,j], 2, dij[j,], "*")
      }
    }
    FD_q_alpha = exp(-1*sum(xx, na.rm=T))/(N^2)
    # Chiu & Chao 2014 p.8, equ 7b
  } else{ #q != 0 or 1
    if( q == 0){ #Chiu & Chao 2014 p.8, (2) when q = 0, ...
      FAD_pair = matrix(0, N, N)
      for(k in 1:N){
        for(m in 1:N){
          s1 = names(comm_alpha[k,][comm_alpha[k,]>0])
          s2 = names(comm_alpha[m,][comm_alpha[m,]>0])
          FAD_pair[k,m] = sum(dij[unique(c(s1, s2)), unique(c(s1, s2))])
        }
      }
      FD_q_alpha = sum(FAD_pair)/(N^2)
    } else{
      x = (outer(comm_alpha, comm_alpha, FUN = "*"))/(Q_gamma * (sum(comm_alpha)^2))
      x = x^q
      # then * dij
      for(k in 1:N){
        for(j in 1:S){
          x[,,k,j] = sweep(x[,,k,j], 2, dij[j,], "*")
        }
      }
      FD_q_alpha = (1/N^2) * (sum(x, na.rm = T)^(1/(1-q)))
      # Chiu & Chao 2014 p.8, equ 7a
    }}

  FD_q_beta = FD_q_gamma / FD_q_alpha

  if(q == 1){
    local_dist_overlap = 1 - ((log(FD_q_gamma) - log(FD_q_alpha))/(2*log(N)))
  } else {
    local_dist_overlap = (N^(2*(1-q)) - FD_q_beta^(1-q))/
      (N^(2*(1-q)) - 1)
  }

  if(q == 1){
    region_dist_overlap = 1 - ((log(FD_q_gamma) - log(FD_q_alpha))/(2*log(N)))
  } else {
    region_dist_overlap = ((1/FD_q_beta)^(1-q) - (1/N)^(2*(1-q)))/
      (1 - (1/N)^(2*(1-q)))
  }

  return(data.frame(q = q,
                    raoQ_gamma = Q_gamma,
                    FD_gamma = FD_q_gamma,
                    FD_alpha = FD_q_alpha,
                    FD_beta = FD_q_beta,
                    local_dist_overlap = local_dist_overlap,
                    region_dist_overlap = region_dist_overlap))
}
