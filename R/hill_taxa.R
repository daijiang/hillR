#' @importFrom stats na.omit dist
NULL

#' Taxonomic diversity through Hill Numbers
#'
#' Calculate taxonomic diversity for each site (alpha diversity).
#'
#' @author Daijiang Li
#'
#' @param comm data frame of vegtation data. Sites as rows, species as columns.
#' @param q hill number, q = 0 (default) to get species richness,
#' q = 1 to get shannon entropy, q = 2 will give inverse Simpson.
#' @param MARGIN default is 1, if sites are columns, set MARGIN to 2.
#' @param base default is exp(1), the base of log.
#' @export
#' @return a named vector, diversity values for each site in the comm.
#' @rdname hill_taxa
#' @examples
#' dummy = FD::dummy
#' hill_taxa(comm = dummy$abun, q = 0)
#' # same as: vegan::specnumber(dummy$abun)
#' hill_taxa(comm = dummy$abun, q = 1)
#' # same as: exp(vegan::diversity(x = dummy$abun, index = "shannon"))
#' hill_taxa(comm = dummy$abun, q = 2)
#' # same as: vegan::diversity(x = dummy$abun, index = "invsimpson")
#' hill_taxa(comm = dummy$abun, q = 0.999)
#'
#'
hill_taxa = function (comm, q = 0, MARGIN = 1, base = exp(1))
{
  comm <- drop(as.matrix(comm))
  if (length(dim(comm)) > 1) { # get relative abundance
    total <- apply(comm, MARGIN, sum)
    comm <- sweep(comm, MARGIN, total, "/")
  }
  else {
    comm <- comm/sum(comm)
  }

  if(q == 0){# richness
    if (length(dim(comm)) > 1) {
      hill <- apply(comm > 0, MARGIN, sum, na.rm = TRUE)
    } else {
      hill <- sum(comm > 0, na.rm = TRUE)
    }
  } else {
    if(q == 1){ # shannon
      comm <- -comm * log(comm, base)
      if (length(dim(comm)) > 1) {
        hill <- exp(apply(comm, MARGIN, sum, na.rm = TRUE))
      } else {
        hill <- exp(sum(comm, na.rm = TRUE))
      }
    } else { # q != 0,1, simpson, etc.
      comm <- comm^q # p_i^q
      if (length(dim(comm)) > 1) {
        hill <- (apply(comm, MARGIN, sum, na.rm = TRUE))^(1/(1-q))
      } else {
        hill <- (sum(comm, na.rm = TRUE))^(1/(1-q))
      }
    }
  }
  hill
}
