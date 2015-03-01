# Functional diversity through Hill Numbers
#' \code{hill_func} to calculate functional diversity for each site.
#'
#' @author Daijiang Li
#'
#' @param comm data frame of vegtation data. Sites as rows, species as columns.
#' @param traits data frame of species functional traits data. Species as rows, traits as columns.
#' It can include both continuous and categorical data. It will be transformed into a distance
#' matrix using `FD::gowdis(traits)`.
#' @param q hill number, q (default is 0) to control weights of species abundance.
#' @param base default is exp(1), the base of log.
#' @param checkdata default is TRUE.
#' @export
#' @return a matrix, with these information for each site: Q (Rao's Q); D_q (functional hill number,
#'  the effective number of equally abundant and functionally equally distince species);
#'  MD_q (mean functional diversity per species, the effective sum of pairwise distances between
#'  a fixed species and all other species); FD_q (total functional diversity, the effective total functional
#'  distance between species of the assemblage). See Chiu and Chao 2014 page 4 for more information.
#'
#' @examples
#' library(FD); data(dummy)
#' hill_func(comm = dummy$abun, traits = dummy$trait, q = 0)
#' hill_func(comm = dummy$abun, traits = dummy$trait, q = 1)
#' hill_func(comm = dummy$abun, traits = dummy$trait, q = 0.9999)
#' hill_func(comm = dummy$abun, traits = dummy$trait, q = 2)
#' hill_func(comm = dummy$abun, traits = dummy$trait, q = 3)
#'
#'
hill_func = function(comm, traits, q = 0, base = exp(1), checkdata=TRUE,...){
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

  if(any(!colnames(comm) %in% rownames(traits))){
    warning("\n There are species from community data that are not on traits matrix\nDelete these species from comm data...\n")
    comm = comm[, colnames(comm) %in% rownames(traits)]
  }

  if(any(!rownames(traits) %in% colnames(comm))){
    warning("\n There are species from trait data that are not on comm matrix\nDelete these species from trait data...\n")
    traits = traits[rownames(traits) %in% colnames(comm), ]
  }

  comm = as.matrix(comm)
  N = nrow(comm)
  S = ncol(comm)
  comm = sweep(comm, 1, rowSums(comm, na.rm = TRUE), "/") # relative abun
  traits.dist = as.matrix(gowdis(x=traits, ...))
  dij = traits.dist^2 # trait distance matrix
  #   inter = comm %*% dij # \sum_i,j_S(p_i * dij)
  #   Q = rowSums(sweep(comm,1,inter,"*", check.margin = F))/2 # \sum_j_S\sum_i,j_S(p_i * dij)
  Q = vector("numeric", length = N)
  names(Q) = dimnames(comm)[[1]]
  for(k in 1:N){
    Q[k] = (comm[k,]) %*% dij %*% matrix(comm[k,], ncol = 1) # /2
    # Q[k] = sum(dij * outer(comm[k,], comm[k,], "*"))/2
  }

  ## D_q
  FD_q = MD_q = D_q = vector("numeric", length = N)
  names(D_q) = dimnames(comm)[[1]]
  names(MD_q) = dimnames(comm)[[1]]
  names(FD_q) = dimnames(comm)[[1]]

  if(q == 0){
    for(k in 1:N){
      df2 = comm[k,][comm[k,] > 0]
      dis2 = dij[names(df2), names(df2)]
      D_q[k] = sum(dis2/Q[k])^0.5
      MD_q[k] = D_q[k] * Q[k]
      FD_q[k] = (D_q[k])^2 * Q[k]
    }
  } else{
    if(q == 1){
      for(k in 1:N){
        df2 = comm[k,][comm[k,] > 0]
        dis2 = dij[names(df2), names(df2)]
        D_q[k] = exp(-0.5 * sum(dis2/Q[k] *
                                  outer(df2, df2, FUN = "*") *
                                  log(outer(df2, df2, FUN = "*"),base)))
        # exp(-0.5 * (dij/Q) * pi*pj * log(pi*pj)
        MD_q[k] = D_q[k] * Q[k]
        FD_q[k] = (D_q[k])^2 * Q[k]
      }
    } else{ #q != 0 or 1
      for(k in 1:N){
        df2 = comm[k,][comm[k,] > 0]
        dis2 = dij[names(df2), names(df2)]
        D_q[k] = sum(dis2/Q[k] * (outer(df2, df2, FUN = "*")^q))^(1/(2*(1-q)))
        MD_q[k] = D_q[k] * Q[k]
        FD_q[k] = (D_q[k])^2 * Q[k]
      }
    }
  }
  rbind(Q, D_q, MD_q, FD_q)
}
