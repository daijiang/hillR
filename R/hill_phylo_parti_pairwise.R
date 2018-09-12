#' Pairwise phylogenetic diversity through Hill numbers
#'
#' Calculate pairwise phylogenetic diversity.
#'
#' @param comm data frame of vegtation data. Sites as rows, species as columns.
#' @param tree a phylogeny with class 'phylo'.
#' @param q hill number, any non-negative value.
#' @param output output type: data.frame (default) or matrix. If matrix, then this function will return a list of matrices.
#' @param pairs full or unique (default). Do you want to compare all possible pairs (i.e. n^2) or just unique pairs (i.e. choose(n, 2))?
#' @param ... additional arguments for \code{hill_func_parti}.
#' @return a data frame or a matrix with results for all pairwise comparisons.
#' @seealso \code{\link{hill_phylo_parti}}
#' @export
#' @examples
#' comm = dummy = FD::dummy$abun
#' tree = ape::rtree(n = ncol(comm), tip.label = paste0('sp', 1:8))
#' hill_phylo_parti_pairwise(comm, tree, q = 0, show.warning = FALSE)
#' hill_phylo_parti_pairwise(comm, tree, q = 0.999, show.warning = FALSE)
#' hill_phylo_parti_pairwise(comm, tree, q = 1, show.warning = FALSE)
#' hill_phylo_parti_pairwise(comm, tree, q = 2, show.warning = FALSE)
#'
hill_phylo_parti_pairwise <- function(comm, tree, q = 0, output = c("data.frame", "matrix"), 
    pairs = c("unique", "full"), ...) {
    output <- match.arg(output)
    pairs <- match.arg(pairs)
    nsite <- nrow(comm)
    temp <- matrix(1, nsite, nsite)
    dimnames(temp) <- list(row.names(comm), row.names(comm))
    gamma_pair <- alpha_pair <- beta_pair <- local_simi <- region_simi <- temp
    for (i in 1:nsite) {
        for (j in i:nsite) {
            o <- hill_phylo_parti(comm[c(i, j), ], tree, q = q)
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
    
    if (pairs == "full") {
        if (output == "matrix") {
            out <- list(q = q, PD_gamma = gamma_pair, PD_alpha = alpha_pair, PD_beta = beta_pair, 
                local_similarity = local_simi, region_similarity = region_simi)
        }
        
        if (output == "data.frame") {
            site.comp <- as.matrix(expand.grid(row.names(comm), row.names(comm)))
            out <- plyr::adply(site.comp, 1, function(x) {
                data.frame(q = q, site1 = x[1], site2 = x[2], PD_gamma = gamma_pair[x[1], 
                  x[2]], PD_alpha = alpha_pair[x[1], x[2]], PD_beta = beta_pair[x[1], 
                  x[2]], local_similarity = local_simi[x[1], x[2]], region_similarity = region_simi[x[1], 
                  x[2]])
            })[, -1]  # get rid of X1 column
            out <- tibble::as.tibble(out)
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
                data.frame(q = q, site1 = x[1], site2 = x[2], PD_gamma = gamma_pair[x[1], 
                  x[2]], PD_alpha = alpha_pair[x[1], x[2]], PD_beta = beta_pair[x[1], 
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
