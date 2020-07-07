#' @importFrom stats na.omit dist
NULL

#' Taxonomic diversity through Hill Numbers
#'
#' Calculate taxonomic diversity for each site (alpha diversity).
#'
#' @param comm A data frame of vegetation data. Sites as rows, species as columns.
#' @param q Hill number, \code{q} = 0 (default) to get species richness,
#'   \code{q} = 1 to get shannon entropy, \code{q} = 2 will give inverse Simpson.
#' @param MARGIN default is 1, if sites are columns, set \code{MARGIN} to 2.
#' @param base default is \code{exp(1)}, the base of log.
#' @export
#' @return A named vector, diversity values for each site in the comm.
#' @rdname hill_taxa
#' @references Chao, Anne, Chun-Huo Chiu, and Lou Jost. Unifying Species Diversity, Phylogenetic Diversity, Functional Diversity, and Related Similarity and Differentiation Measures Through Hill Numbers. Annual Review of Ecology, Evolution, and Systematics 45, no. 1 (2014): 297–324. <doi:10.1146/annurev-ecolsys-120213-091540>.
#'
#' Jost, Lou. Entropy and diversity. Oikos 113, no. 2 (2006): 363-375. <doi:10.1111/j.2006.0030-1299.14714.x>.
#' @examples
#' dummy = FD::dummy
#' hill_taxa(comm = dummy$abun, q = 0)
#' # same as: vegan::specnumber(dummy$abun)
#' hill_taxa(comm = dummy$abun, q = 1)
#' # same as: exp(vegan::diversity(x = dummy$abun, index = 'shannon'))
#' hill_taxa(comm = dummy$abun, q = 2)
#' # same as: vegan::diversity(x = dummy$abun, index = 'invsimpson')
#' hill_taxa(comm = dummy$abun, q = 0.999)
#'
hill_taxa <- function(comm, q = 0, MARGIN = 1, base = exp(1)) {
    comm <- drop(as.matrix(comm))
    if (length(dim(comm)) > 1) {
        # get relative abundance
        total <- apply(comm, MARGIN, sum)
        comm <- sweep(comm, MARGIN, total, "/")
    } else {
        comm <- comm/sum(comm)
    }

    if (q == 0) {
        # richness
        if (length(dim(comm)) > 1) {
            hill <- apply(comm > 0, MARGIN, sum, na.rm = TRUE)
        } else {
            hill <- sum(comm > 0, na.rm = TRUE)
        }
    } else {
        if (q == 1) {
            # shannon
            comm <- -comm * log(comm, base)
            if (length(dim(comm)) > 1) {
                hill <- exp(apply(comm, MARGIN, sum, na.rm = TRUE))
            } else {
                hill <- exp(sum(comm, na.rm = TRUE))
            }
        } else {
            # q != 0,1, simpson, etc.
            comm <- comm^q  # p_i^q
            if (length(dim(comm)) > 1) {
                hill <- (apply(comm, MARGIN, sum, na.rm = TRUE))^(1/(1 - q))
            } else {
                hill <- (sum(comm, na.rm = TRUE))^(1/(1 - q))
            }
        }
    }
    hill
}
