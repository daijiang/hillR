#' Decompostion of Taxonomic diversity through Hill Numbers
#'
#' Calculate taxonomic gamma, alpha, and beta diversity across all communities, as
#' well as site similarity. If comm has 2 sites, this function gives pair comparison.
#' If comm has >2 sites, gamma diversity is the diversity of the pooled assemblage,
#' alpha is the average diversity across all site, beta is across all communities.
#'
#' @inheritParams hill_taxa
#' @inheritParams hill_func
#' @param rel_then_pool default is \code{TRUE.} Abundance of species are first changed to relative abundance within sites,
#' then pooled into one assemblage. If \code{FALSE}, sites are pooled first, then change abundance of species
#' to relative abundance.
#' @param show_warning whether to print warning, default is \code{TRUE}.
#' @export
#' @return A data frame with one row (across all sites), including these columns: q, gamma diversity, alpha diverisity,
#' beta diversity, MacArthur's homogeneity measure, local similarity (species overlap, similar to Sorensen),
#' and region similarity (species overlap, similar to Jaccard).
#' See Chao, Chiu and Jost 2014 Table 2 for more information.
#' @references Chao, Anne, Chun-Huo Chiu, and Lou Jost. Unifying Species Diversity, Phylogenetic Diversity, Functional Diversity, and Related Similarity and Differentiation Measures Through Hill Numbers. Annual Review of Ecology, Evolution, and Systematics 45, no. 1 (2014): 297–324. <doi:10.1146/annurev-ecolsys-120213-091540>.
#'
#' Jost, Lou. Entropy and diversity. Oikos 113, no. 2 (2006): 363-375. <doi:10.1111/j.2006.0030-1299.14714.x>.
#' @seealso \code{\link{hill_taxa}}
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
hill_taxa_parti <- function(comm, q = 0, base = exp(1), rel_then_pool = TRUE,
                            show_warning = TRUE, check_data = TRUE) {
    if(check_data){
        if (any(comm < 0))
            stop("Negative value in comm data")
        if (any(colSums(comm) == 0) & show_warning)
            warning("Some species in comm data were not observed in any site,\n delete them...")
    }

    comm <- comm[, colSums(comm) != 0, drop = FALSE]
    N <- nrow(comm)
    S <- ncol(comm)
    comm <- as.matrix(comm)
    if (rel_then_pool) {
        comm_gamma <- colSums(sweep(comm, 1, rowSums(comm, na.rm = TRUE), "/"))/N
        # relative abun
    } else {
        comm_gamma <- colSums(comm)/sum(comm)
    }
    if (!all.equal(sum(comm_gamma), 1))
        stop("Accumlative relative abundance should be 1")

    if (rel_then_pool) {
        comm_alpha <- sweep(comm, 1, rowSums(comm, na.rm = TRUE), "/")  # relative abun
    } else {
        comm_alpha <- comm
    }

    TD_q_gamma <- hill_taxa(comm_gamma, q = q)

    # TD_q_alpha
    if (q == 0) {
        TD_q_alpha <- sum(hill_taxa(comm_alpha, q = 0))/N  # mean of sp richness
    } else {
        if (q == 1) {
            TD_q_alpha <- exp(-1 * sum((comm_alpha/sum(comm_alpha)) * log((comm_alpha/sum(comm_alpha)),
                base), na.rm = T))/N
        } else {
            TD_q_alpha <- (1/N) * (sum((comm_alpha/sum(comm_alpha))^q)^(1/(1 - q)))
        }
    }

    TD_q_beta <- TD_q_gamma/TD_q_alpha

    if (q == 1) {
        # why is negative and abs()>1 ??
        if (rel_then_pool) {
            local_taxa_overlap <- (log(N, base) - log(TD_q_gamma) + log(TD_q_alpha))/log(N,
                base)
        } else {
            local_taxa_overlap <- (log(TD_q_alpha, base) - log(TD_q_gamma, base) - sum((rowSums(comm_alpha)/sum(comm_alpha)) *
                log(rowSums(comm_alpha)/sum(comm_alpha), base)))/log(N, base)
        }
    } else {
        local_taxa_overlap <- (N^(1 * (1 - q)) - TD_q_beta^(1 - q))/(N^(1 * (1 - q)) -
            1)
    }

    if (q == 1) {
        region_taxa_overlap <- local_taxa_overlap
    } else {
        region_taxa_overlap <- ((1/TD_q_beta)^(1 - q) - (1/N)^(1 * (1 - q)))/(1 - (1/N)^(1 *
            (1 - q)))
    }

    return(data.frame(q = q, TD_gamma = TD_q_gamma, TD_alpha = TD_q_alpha, TD_beta = TD_q_beta,
        M_homog = 1/TD_q_beta, local_similarity = local_taxa_overlap, region_similarity = region_taxa_overlap))
}
