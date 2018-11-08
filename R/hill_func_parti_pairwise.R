#' Pairwise comparisons for all sites.
#'
#' Calculate pairwise functional gamma, alpha, and beta diversity for communities, as
#'  well as site similarity. It is based on \code{\link{hill_func_parti}}.
#'  If comm has >2 sites, this function will give results for all pairwise comparisons.
#'
#' @inheritParams hill_func_parti
#' @inheritParams hill_taxa_parti_pairwise
#' @param ... additional arguments for \code{hill_func_parti}.
#' @export
#' @references Chao, Anne, Chun-Huo Chiu, and Lou Jost. Unifying Species Diversity, Phylogenetic Diversity, Functional Diversity, and Related Similarity and Differentiation Measures Through Hill Numbers. Annual Review of Ecology, Evolution, and Systematics 45, no. 1 (2014): 297–324. <doi:10.1146/annurev-ecolsys-120213-091540>.
#'
#' Chiu, Chun-Huo, and Anne Chao. Distance-Based Functional Diversity Measures and Their Decomposition: A Framework Based on Hill Numbers. PLoS ONE 9, no. 7 (July 7, 2014): e100014. <doi:10.1371/journal.pone.0100014>.
#' @return a data frame with results for all pairwise comparisons.
#' @seealso \code{\link{hill_func_parti}}
#' @examples
#' \dontrun{
#' dummy = FD::dummy
#' hill_func_parti_pairwise(comm = dummy$abun, traits = dummy$trait, q = 0)
#' hill_func_parti_pairwise(comm = dummy$abun, traits = dummy$trait, q = 0,
#'                          output = 'matrix')
#' hill_func_parti_pairwise(comm = dummy$abun, traits = dummy$trait, q = 0,
#'                          output = 'matrix', pairs = 'full')
#' hill_func_parti_pairwise(comm = dummy$abun, traits = dummy$trait, q = 1)
#' hill_func_parti_pairwise(comm = dummy$abun, traits = dummy$trait, q = 0.9999)
#' hill_func_parti_pairwise(comm = dummy$abun, traits = dummy$trait, q = 2)
#' hill_func_parti_pairwise(comm = dummy$abun, traits = dummy$trait, q = 3)
#' }
hill_func_parti_pairwise <- function(comm, traits, traits_as_is = FALSE, q = 0, rel_then_pool = TRUE,
    output = c("data.frame", "matrix"), pairs = c("unique", "full"), ...) {
    output <- match.arg(output)
    pairs <- match.arg(pairs)
    nsite <- nrow(comm)
    temp <- matrix(1, nsite, nsite)
    dimnames(temp) <- list(row.names(comm), row.names(comm))
    gamma_pair <- alpha_pair <- beta_pair <- local_simi <- region_simi <- temp
    for (i in 1:nsite) {
        for (j in i:nsite) {
            o <- hill_func_parti(comm = comm[c(i, j), ], traits = traits, traits_as_is = traits_as_is,
                q = q, rel_then_pool = rel_then_pool, ...)
            gamma_pair[i, j] <- o$FD_gamma
            gamma_pair[j, i] <- o$FD_gamma
            alpha_pair[i, j] <- o$FD_alpha
            alpha_pair[j, i] <- o$FD_alpha
            beta_pair[i, j] <- o$FD_beta
            beta_pair[j, i] <- o$FD_beta
            local_simi[i, j] <- o$local_similarity
            local_simi[j, i] <- o$local_similarity
            region_simi[i, j] <- o$region_similarity
            region_simi[j, i] <- o$region_similarity
        }
    }

    if (pairs == "full") {
        if (output == "matrix") {
            out <- list(q = q, FD_gamma = gamma_pair, FD_alpha = alpha_pair, FD_beta = beta_pair,
                local_similarity = local_simi, region_similarity = region_simi)
        }

        if (output == "data.frame") {
            site.comp <- as.matrix(expand.grid(row.names(comm), row.names(comm)))
            out <- plyr::adply(site.comp, 1, function(x) {
                data.frame(q = q, site1 = x[1], site2 = x[2], FD_gamma = gamma_pair[x[1],
                  x[2]], FD_alpha = alpha_pair[x[1], x[2]], FD_beta = beta_pair[x[1],
                  x[2]], local_similarity = local_simi[x[1], x[2]], region_similarity = region_simi[x[1],
                  x[2]])
            })[, -1]
            out <- tibble::as.tibble(out)
        }
        out
    }

    if (pairs == "unique") {
        gamma_pair[lower.tri(gamma_pair, diag = TRUE)] <- NA
        alpha_pair[lower.tri(alpha_pair, diag = TRUE)] <- NA
        beta_pair[lower.tri(beta_pair, diag = TRUE)] <- NA
        local_simi[lower.tri(local_simi, diag = TRUE)] <- NA
        region_simi[lower.tri(region_simi, diag = TRUE)] <- NA

        if (output == "matrix") {
            out <- list(q = q, FD_gamma = gamma_pair, FD_alpha = alpha_pair, FD_beta = beta_pair,
                local_similarity = local_simi, region_similarity = region_simi)
        }

        if (output == "data.frame") {
            site.comp <- as.matrix(expand.grid(row.names(comm), row.names(comm)))
            out <- plyr::adply(site.comp, 1, function(x) {
                data.frame(q = q, site1 = x[1], site2 = x[2], FD_gamma = gamma_pair[x[1],
                  x[2]], FD_alpha = alpha_pair[x[1], x[2]], FD_beta = beta_pair[x[1],
                  x[2]], local_similarity = local_simi[x[1], x[2]], region_similarity = region_simi[x[1],
                  x[2]])
            })
            out <- na.omit(out)[, -1]
            row.names(out) <- NULL
            out <- tibble::as.tibble(out)
        }
    }
    out
}
