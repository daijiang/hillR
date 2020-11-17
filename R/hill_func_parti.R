#' Decompostion of functional diversity through Hill Numbers
#'
#' Calculate functional gamma, alpha, and beta diversity for all communities, as
#'  well as site similarity. These values are based on ALL communities.
#'
#' @inheritParams hill_taxa
#' @inheritParams hill_func
#' @inheritParams hill_taxa_parti
#' @export
#' @references Chao, Anne, Chun-Huo Chiu, and Lou Jost. Unifying Species Diversity, Phylogenetic Diversity, Functional Diversity, and Related Similarity and Differentiation Measures Through Hill Numbers. Annual Review of Ecology, Evolution, and Systematics 45, no. 1 (2014): 297–324. <doi:10.1146/annurev-ecolsys-120213-091540>.
#'
#' Chiu, Chun-Huo, and Anne Chao. Distance-Based Functional Diversity Measures and Their Decomposition: A Framework Based on Hill Numbers. PLoS ONE 9, no. 7 (July 7, 2014): e100014. <doi:10.1371/journal.pone.0100014>.
#' @return  a data frame with one row (across all sites), including these columns: q, RaoQ of pooled assemblage,
#' gamma diversity, alpha diverisity, beta diversity, local species overlap (similar to Sorensen), and region species
#' overlap (similar to Jaccard). See Chiu and Chao 2014 Table 3 for more information.
#' @seealso \code{\link{hill_taxa_parti}}, \code{\link{hill_func}}
#'
#' @examples
#' dummy = FD::dummy
#' hill_func_parti(comm = dummy$abun, traits = dummy$trait, q = 0)
#' hill_func_parti(comm = dummy$abun, traits = dummy$trait, q = 1)
#' hill_func_parti(comm = dummy$abun, traits = dummy$trait, q = 0.9999)
#' hill_func_parti(comm = dummy$abun, traits = dummy$trait, q = 2)
#' hill_func_parti(comm = dummy$abun, traits = dummy$trait, q = 3)
#'
hill_func_parti <- function(comm, traits, traits_as_is = FALSE, q = 0, base = exp(1),
    check_data = TRUE, rel_then_pool = TRUE, ord = c("podani", "metric"), stand_dij = FALSE,
    show_warning = TRUE) {
    if (check_data) {
        if (any(comm < 0))
            stop("Negative value in comm data")
        if (is.null(rownames(traits))) {
            stop("\n Traits have no row names\n")
        }
        if (is.null(colnames(comm))) {
            stop("\n Comm data have no col names\n")
        }
        if (any(colSums(comm) == 0) & show_warning)
            warning("Some species in comm data were not observed in any site,\n
                                      delete them...")
    }

    comm <- comm[, colSums(comm) != 0]

    if (any(!colnames(comm) %in% rownames(traits))) {
        warning("\n There are species from community data that are not on traits matrix\n
            Delete these species from comm data...\n")
        comm <- comm[, colnames(comm) %in% rownames(traits)]
    }

    if (traits_as_is) {
        if (any(!rownames(traits) %in% colnames(comm))) {
            if (show_warning)
                warning("\n There are species from trait data that are not in comm matrix\n
              Delete these species from trait data...\n")
            traits <- traits[rownames(traits) %in% colnames(comm), colnames(traits) %in%
                colnames(comm)]
        }
        dij <- as.matrix(traits)
    } else {
        # traits is not a distance matrix
        traits <- traits[colnames(comm), , drop = FALSE]

        if (ncol(traits) == 1) {
            # only 1 trait
            if (any(is.na(traits)) & show_warning) {
                if (show_warning)
                  warning("Warning: Species with missing trait values have been excluded.",
                    "\n")
                traits <- na.omit(traits)
                comm <- comm[, colnames(comm) %in% rownames(traits)]
            }
            if (is.numeric(traits[, 1])) {
                # 1 numeric trait
                dij <- dist(traits)
            }
            if (is.factor(traits[, 1]) | is.character(traits[, 1])) {
                # 1 categorical trait
                if (is.ordered(traits[, 1])) {
                  traits2 <- data.frame(rank(traits[, 1]))
                  rownames(traits2) <- rownames(traits)
                  names(traits2) <- names(traits)
                  dij <- dist(traits2)
                } else {
                  traits[, 1] <- as.factor(traits[, 1])
                  x.f <- as.factor(traits[, 1])
                  x.dummy <- diag(nlevels(x.f))[x.f, ]
                  x.dummy.df <- data.frame(x.dummy, row.names = rownames(traits))
                  dij <- ade4::dist.binary(x.dummy.df, method = 2)
                }
            }
        } else {
            # more than 1 trait:
            for (i in 1:ncol(traits)) {
                if (is.factor(traits[, i]) & nlevels(traits[, i]) == 2) {
                  traits[, i] <- as.numeric(traits[, i]) - 1  # so to be 0, 1
                }
            }
            if (all(sapply(traits, is.numeric)) & all(!is.na(traits))) {
                dij <- dist(scale(traits, center = TRUE, scale = TRUE))
            } else {
                ord <- match.arg(ord)
                dij <- FD::gowdis(x = traits, asym.bin = NULL, ord = ord)
            }
            # dij = gowdis(x=traits, ...)
        }
    }

    comm <- as.matrix(comm)
    N <- nrow(comm)
    S <- ncol(comm)

    dij <- as.matrix(dij)
    if(any(!is.finite(dij))){
        warning("Some species pairs have distance of NA or NaN, set it to zero (this may be incorrect!)")
        dij[!is.finite(dij)] <- 0
    }

    if (stand_dij)
        dij <- dij/max(dij)

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

    Q_gamma <- as.vector(comm_gamma %*% dij %*% matrix(comm_gamma, ncol = 1))

    ## FD_q_gamma
    if (q == 1) {
        if (Q_gamma == 0) {
            FD_q_gamma <- 0
        } else {
            FD_q_gamma <- exp(-1 * sum(dij * (outer(comm_gamma, comm_gamma, FUN = "*")/Q_gamma) *
                log(outer(comm_gamma, comm_gamma, FUN = "*")/Q_gamma)))
        }
        # Chiu & Chao 2014 p.7, equ 6b q != 0 or 1
    } else {
        if (Q_gamma == 0) {
            FD_q_gamma <- 0
        } else {
            FD_q_gamma <- sum(dij * ((outer(comm_gamma, comm_gamma, FUN = "*")/Q_gamma)^q))^(1/(1 -
                q))
        }
        # Chiu & Chao 2014 p.7, equ 6a
    }

    ## FD_q_alpha if q_gamma is 0, no need to calc alpha
    if (Q_gamma == 0) {
        FD_q_alpha <- 1e-05
    } else {
        if (q == 1) {
            x <- (outer(comm_alpha, comm_alpha, FUN = "*"))/(Q_gamma * (sum(comm_alpha)^2))
            x[x == 0] <- NA
            xx <- x * log(x, base)
            # then * dij
            for (k in 1:N) {
                for (j in 1:S) {
                  xx[, , k, j] <- sweep(xx[, , k, j], 2, dij[j, ], "*")
                }
            }
            FD_q_alpha <- exp(-1 * sum(xx, na.rm = T))/(N^2)
            # Chiu & Chao 2014 p.8, equ 7b q != 0 or 1 Chiu & Chao 2014 p.8, (2) when q = 0, ...
        } else {
            if (q == 0) {
                FAD_pair <- matrix(0, N, N)
                for (k in 1:N) {
                  for (m in 1:N) {
                    s1 <- names(comm_alpha[k, ][comm_alpha[k, ] > 0])
                    s2 <- names(comm_alpha[m, ][comm_alpha[m, ] > 0])
                    FAD_pair[k, m] <- sum(dij[unique(c(s1, s2)), unique(c(s1, s2))])
                  }
                }
                FD_q_alpha <- sum(FAD_pair)/(N^2)
            } else {
                x <- (outer(comm_alpha, comm_alpha, FUN = "*"))/(Q_gamma * (sum(comm_alpha)^2))
                x <- x^q
                # then * dij
                for (k in 1:N) {
                  for (j in 1:S) {
                    x[, , k, j] <- sweep(x[, , k, j], 2, dij[j, ], "*")
                  }
                }
                FD_q_alpha <- (1/N^2) * (sum(x, na.rm = T)^(1/(1 - q)))
                # Chiu & Chao 2014 p.8, equ 7a
            }
        }
    }

    FD_q_beta <- FD_q_gamma/FD_q_alpha

    if (q == 1) {
        local_dist_overlap <- 1 - ((log(FD_q_gamma) - log(FD_q_alpha))/(2 * log(N)))
    } else {
        local_dist_overlap <- (N^(2 * (1 - q)) - FD_q_beta^(1 - q))/(N^(2 * (1 - q)) -
            1)
    }

    if (q == 1) {
        region_dist_overlap <- 1 - ((log(FD_q_gamma) - log(FD_q_alpha))/(2 * log(N)))
    } else {
        region_dist_overlap <- ((1/FD_q_beta)^(1 - q) - (1/N)^(2 * (1 - q)))/(1 - (1/N)^(2 *
            (1 - q)))
    }

    return(data.frame(q = q, raoQ_gamma = Q_gamma, FD_gamma = FD_q_gamma, FD_alpha = FD_q_alpha,
        FD_beta = FD_q_beta, local_similarity = local_dist_overlap, region_similarity = region_dist_overlap))
}
