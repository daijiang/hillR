#' Pairwise comparisons for all sites
#'
#' Calculate pairwise taxonomic gamma, alpha, and beta diversity for communities, as
#' well as site similarity. It is based on \code{\link{hill_taxa_parti}}.
#' If comm has >2 sites, this function will give results for all pairwise comparisons.
#'
#' @inheritParams hill_taxa
#' @inheritParams hill_taxa_parti
#' @param output output type: data.frame (default) or matrix. If matrix, then this function will return a list of matrices.
#' @param pairs full or unique (default). Do you want to compare all possible pairs (i.e. n^2) or just unique pairs (i.e. \code{choose(n, 2))}?
#' @param .progress Whether to show progress bar. Default is `TRUE`.
#' @param ... other arguments in \code{hill_taxa_parti()}.
#' @export
#' @return A data frame with results for all pairwise comparisons.
#' @references Chao, Anne, Chun-Huo Chiu, and Lou Jost. Unifying Species Diversity, Phylogenetic Diversity, Functional Diversity, and Related Similarity and Differentiation Measures Through Hill Numbers. Annual Review of Ecology, Evolution, and Systematics 45, no. 1 (2014): 297–324. <doi:10.1146/annurev-ecolsys-120213-091540>.
#'
#' Jost, Lou. Entropy and diversity. Oikos 113, no. 2 (2006): 363-375. <doi:10.1111/j.2006.0030-1299.14714.x>.
#' @seealso \code{\link{hill_taxa_parti}}
#' @examples
#' \dontrun{
#' dummy = FD::dummy
#' hill_taxa_parti_pairwise(comm = dummy$abun, q = 0)
#' hill_taxa_parti_pairwise(comm = dummy$abun, q = 0, output = 'matrix')
#' hill_taxa_parti_pairwise(comm = dummy$abun, q = 1)
#' hill_taxa_parti_pairwise(comm = dummy$abun, q = 0.9999999)
#' hill_taxa_parti_pairwise(comm = dummy$abun, q = 0.9999999, rel_then_pool = FALSE)
#' hill_taxa_parti_pairwise(comm = dummy$abun, q = 1, rel_then_pool = FALSE)
#' hill_taxa_parti_pairwise(comm = dummy$abun, q = 2)
#' hill_taxa_parti_pairwise(comm = dummy$abun, q = 3)
#' }
hill_taxa_parti_pairwise <- function(comm, q = 0, rel_then_pool = TRUE,
                                     output = c("data.frame", "matrix"),
                                     pairs = c("unique", "full"),
                                     .progress = TRUE,
                                     show_warning = TRUE, ...) {
    if (any(comm < 0))
        stop("Negative value in comm data")
    if (any(colSums(comm) == 0) & show_warning)
        warning("Some species in comm data were not observed in any site,\n delete them...")

    output <- match.arg(output)
    pairs <- match.arg(pairs)
    nsite <- nrow(comm)
    temp <- matrix(1, nsite, nsite)
    dimnames(temp) <- list(row.names(comm), row.names(comm))
    gamma_pair <- alpha_pair <- beta_pair <- local_simi <- region_simi <- temp
    if(.progress)
        progbar = utils::txtProgressBar(min = 0, max = nsite - 1, initial = 0, style = 3)
    for (i in 1:(nsite - 1)) {
        if(.progress) utils::setTxtProgressBar(progbar, i)
        for (j in (i + 1):nsite) {
            o <- hill_taxa_parti(comm[c(i, j), ], q = q, check_data = FALSE, ...)
            gamma_pair[i, j] <- o$TD_gamma
            gamma_pair[j, i] <- o$TD_gamma
            alpha_pair[i, j] <- o$TD_alpha
            alpha_pair[j, i] <- o$TD_alpha
            beta_pair[i, j] <- o$TD_beta
            beta_pair[j, i] <- o$TD_beta
            local_simi[i, j] <- o$local_similarity
            local_simi[j, i] <- o$local_similarity
            region_simi[i, j] <- o$region_similarity
            region_simi[j, i] <- o$region_similarity
        }
    }
    if(.progress) close(progbar)

    if (pairs == "full") {
        if (output == "matrix") {
            out <- list(q = q, TD_gamma = gamma_pair, TD_alpha = alpha_pair, TD_beta = beta_pair,
                local_similarity = local_simi, region_similarity = region_simi)
        }

        if (output == "data.frame") {
            site.comp <- as.matrix(expand.grid(row.names(comm), row.names(comm)))
            out <- plyr::adply(site.comp, 1, function(x) {
                data.frame(q = q, site1 = x[1], site2 = x[2],
                           TD_gamma = gamma_pair[x[1], x[2]],
                           TD_alpha = alpha_pair[x[1], x[2]],
                           TD_beta = beta_pair[x[1], x[2]],
                           local_similarity = local_simi[x[1], x[2]],
                           region_similarity = region_simi[x[1], x[2]])
            })[, -1]  # get rid of X1 column
            out <- tibble::as_tibble(out)
        }
    }

    if (pairs == "unique") {
        gamma_pair[lower.tri(gamma_pair, diag = TRUE)] <- NA
        alpha_pair[lower.tri(alpha_pair, diag = TRUE)] <- NA
        beta_pair[lower.tri(beta_pair, diag = TRUE)] <- NA
        local_simi[lower.tri(local_simi, diag = TRUE)] <- NA
        region_simi[lower.tri(region_simi, diag = TRUE)] <- NA

        if (output == "matrix") {
            out <- list(q = q, TD_gamma = gamma_pair, TD_alpha = alpha_pair, TD_beta = beta_pair,
                local_similarity = local_simi, region_similarity = region_simi)
        }

        if (output == "data.frame") {
            site.comp <- as.matrix(expand.grid(row.names(comm), row.names(comm)))
            out <- plyr::adply(site.comp, 1, function(x) {
                data.frame(q = q, site1 = x[1], site2 = x[2],
                           TD_gamma = gamma_pair[x[1], x[2]],
                           TD_alpha = alpha_pair[x[1], x[2]],
                           TD_beta = beta_pair[x[1], x[2]],
                           local_similarity = local_simi[x[1], x[2]],
                           region_similarity = region_simi[x[1], x[2]])
            })
            out <- na.omit(out)[, -1]
            row.names(out) <- NULL
            out <- tibble::as_tibble(out)
        }
    }
    out
}

