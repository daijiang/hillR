#' Pairwise phylogenetic diversity through Hill numbers
#'
#' Calculate pairwise phylogenetic diversity.
#'
#' @inheritParams hill_phylo
#' @inheritParams hill_taxa_parti_pairwise
#' @param ... additional arguments for \code{hill_func_parti}.
#' @return A data frame or a matrix with results for all pairwise comparisons.
#' @seealso \code{\link{hill_phylo_parti}}
#' @references Chao, Anne, Chun-Huo Chiu, and Lou Jost. Unifying Species Diversity, Phylogenetic Diversity, Functional Diversity, and Related Similarity and Differentiation Measures Through Hill Numbers. Annual Review of Ecology, Evolution, and Systematics 45, no. 1 (2014): 297–324. <doi:10.1146/annurev-ecolsys-120213-091540>.
#' @export
#' @examples
#' \dontrun{
#' comm = dummy = FD::dummy$abun
#' tree = ape::rtree(n = ncol(comm), tip.label = paste0('sp', 1:8))
#' hill_phylo_parti_pairwise(comm, tree, q = 0, show_warning = FALSE)
#' hill_phylo_parti_pairwise(comm, tree, q = 0.999, show_warning = FALSE)
#' hill_phylo_parti_pairwise(comm, tree, q = 1, show_warning = FALSE)
#' hill_phylo_parti_pairwise(comm, tree, q = 2, show_warning = FALSE)
#' }
hill_phylo_parti_pairwise <- function(comm, tree, q = 0, output = c("data.frame", "matrix"),
    pairs = c("unique", "full"), rel_then_pool = TRUE, .progress = TRUE,
    show_warning = TRUE, ...) {
    if (any(comm < 0)) stop("Negative value in comm data")
    if (class(tree) != "phylo")
        stop("tree must be an object with phylo as class")
    # clean phylogeny and community data
    sp_drop <- setdiff(tree$tip.label, colnames(comm))
    if (length(sp_drop)) {
        if (show_warning)
            warning("Some species in the phylogeny but not in comm, \n remove them from the phylogeny...",
                    immediate. = TRUE)
        tree <- ape::drop.tip(tree, sp_drop)
    }
    if (length(setdiff(colnames(comm), tree$tip.label))) {
        if (show_warning)
            warning("Some species in the comm but not in the phylogeny, \n remove them from the comm",
                    immediate. = TRUE)
    }
    comm <- comm[, tree$tip.label] # only select species in the tree, in that order
    comm <- as.matrix(comm)
    if (rel_then_pool) {
        comm <- sweep(comm, 1, rowSums(comm, na.rm = TRUE), "/")  # relative abun
    }
    pabund <- dat_prep_phylo(comm, tree) # pre-calculate node/tip by site matrix

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
            o <- hill_phylo_parti(comm = comm[c(i, j), ], tree, q = q,
                                  phy_abund = pabund, check_data = FALSE)
            gamma_pair[i, j] <- o$PD_gamma
            gamma_pair[j, i] <- o$PD_gamma
            alpha_pair[i, j] <- o$PD_alpha
            alpha_pair[j, i] <- o$PD_alpha
            beta_pair[i, j] <- o$PD_beta
            beta_pair[j, i] <- o$PD_beta
            local_simi[i, j] <- o$local_similarity
            local_simi[j, i] <- o$local_similarity
            region_simi[i, j] <- o$region_similarity
            region_simi[j, i] <- o$region_similarity
        }
    }
    if(.progress) close(progbar)

    if (pairs == "full") {
        if (output == "matrix") {
            out <- list(q = q, PD_gamma = gamma_pair, PD_alpha = alpha_pair, PD_beta = beta_pair,
                local_similarity = local_simi, region_similarity = region_simi)
        }

        if (output == "data.frame") {
            site.comp <- as.matrix(expand.grid(row.names(comm), row.names(comm)))
            out <- plyr::adply(site.comp, 1, function(x) {
                data.frame(q = q, site1 = x[1], site2 = x[2],
                           PD_gamma = gamma_pair[x[1], x[2]],
                           PD_alpha = alpha_pair[x[1], x[2]],
                           PD_beta = beta_pair[x[1], x[2]],
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
            out <- list(q = q, PD_gamma = gamma_pair, PD_alpha = alpha_pair, PD_beta = beta_pair,
                local_similarity = local_simi, region_similarity = region_simi)
        }

        if (output == "data.frame") {
            site.comp <- as.matrix(expand.grid(row.names(comm), row.names(comm)))
            out <- plyr::adply(site.comp, 1, function(x) {
                data.frame(q = q, site1 = x[1], site2 = x[2],
                           PD_gamma = gamma_pair[x[1], x[2]],
                           PD_alpha = alpha_pair[x[1], x[2]],
                           PD_beta = beta_pair[x[1], x[2]],
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
