#' Functional diversity through Hill Numbers
#'
#' Calculate functional diversity for each site (alpha diversity).
#'
#' @inheritParams hill_taxa
#' @param traits A data frame of species functional traits data. Species as rows, traits as columns.
#' It can include both continuous and categorical data. It will be transformed into a distance
#' matrix using `FD::gowdis(traits)`. If all traits are numeric, then it will use Euclidean distance.
#' @param traits_as_is if \code{FALSE} (default) traits data frame will be transformed into a distance
#' matrix. Otherwise, will use as is (i.e. traits is a symmetric distance matrix).
#' @param check_data whether to check data first? Default is \code{TRUE}.
#' @param div_by_sp as FD calculated in this way will be highly correlated with taxonomic diversity,
#' one potential simple way to correct this is to divide the results by the number of species.
#' However, a more common way to deal with correlations is to use null models and calculate standardized effect sizes.
#' Therefore, I set the default to be \code{FALSE}.
#' @param ord ord in \code{FD::gowdis}.
#' @param fdis whether to calculated FDis, default is \code{TRUE}
#' @param stand_dij whether to standardize distance matrix to have max value of 1? Default is \code{FALSE}.
#' @export
#' @references Chao, Anne, Chun-Huo Chiu, and Lou Jost. Unifying Species Diversity, Phylogenetic Diversity, Functional Diversity, and Related Similarity and Differentiation Measures Through Hill Numbers. Annual Review of Ecology, Evolution, and Systematics 45, no. 1 (2014): 297–324. <doi:10.1146/annurev-ecolsys-120213-091540>.
#'
#' Chiu, Chun-Huo, and Anne Chao. Distance-Based Functional Diversity Measures and Their Decomposition: A Framework Based on Hill Numbers. PLoS ONE 9, no. 7 (July 7, 2014): e100014. <doi:10.1371/journal.pone.0100014>.
#' @return A matrix, with these information for each site: Q (Rao's Q); D_q (functional hill number,
#'  the effective number of equally abundant and functionally equally distinct species);
#'  MD_q (mean functional diversity per species, the effective sum of pairwise distances between
#'  a fixed species and all other species); FD_q (total functional diversity, the effective total functional
#'  distance between species of the assemblage). See Chiu and Chao 2014 page 4 for more information.
#'
#' @examples
#' dummy = FD::dummy
#' hill_func(comm = dummy$abun, traits = dummy$trait, q = 0)
#' hill_func(comm = dummy$abun, traits = dummy$trait, q = 1)
#' hill_func(comm = dummy$abun, traits = dummy$trait, q = 0.9999)
#' hill_func(comm = dummy$abun, traits = dummy$trait, q = 2)
#' hill_func(comm = dummy$abun, traits = dummy$trait, q = 3)
#'
hill_func <- function(comm, traits, traits_as_is = FALSE, q = 0, base = exp(1), check_data = TRUE,
                      div_by_sp = FALSE, ord = c("podani", "metric"), fdis = TRUE, stand_dij = FALSE) {
    if (check_data) {
        if (any(comm < 0))
            stop("Negative value in comm data")
        if(traits_as_is){
            if(is.null(attributes(traits)$Labels)) stop("\n Traits distance matrix has no labels\n")
        } else {
            if (is.null(rownames(traits))) stop("\n Traits have no row names\n")
        }
        if (is.null(colnames(comm))) {
            stop("\n Comm data have no col names\n")
        }
    }

    if(traits_as_is){
        trait_sp = attributes(traits)$Labels
    } else {
        trait_sp = rownames(traits)
    }
    if (any(!colnames(comm) %in% trait_sp)) {
        warning("\n There are species from community data that are not on traits matrix\n
                Delete these species from comm data...\n")
        comm <- comm[, trait_sp]
    }

    # all(rownames(traits) == names(comm))

    if (traits_as_is) {
        # traits is already a distance matrix
        dij <- traits
        if(!inherits(dij, "dist")) stop("`traits` is not a distance object yet `trait_as_is` is TRUE\n")
    } else {
        # traits is not a distance matrix
        traits <- traits[trait_sp, , drop = FALSE]

        if (ncol(traits) == 1) {
            # only 1 trait
            if (any(is.na(traits))) {
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

        # if (!is.euclid(dij)) { if (corr == 'lingoes') { dij2 <- lingoes(dij)
        # warning('Species x species distance matrix was not Euclidean. Lingoes correction was
        # applied.','\n') } if (corr == 'cailliez') { dij2 <- cailliez(dij) warning('Species
        # x species distance matrix was not Euclidean. Cailliez correction was
        # applied.','\n') } if (corr == 'sqrt') { dij2 <- sqrt(dij) # check if sqrt
        # correction actually worked if(!is.euclid(dij2) ) stop('Species x species distance
        # matrix was still is not Euclidean after 'sqrt' correction. Use another correction
        # method.','\n') if (is.euclid(dij2) ) warning('Species x species distance matrix was
        # not Euclidean. 'sqrt' correction was applied.','\n') } if (corr == 'none') { dij2
        # <- quasieuclid(dij) warning('Species x species distance was not Euclidean, but no
        # correction was applied. Only the PCoA axes with positive eigenvalues were
        # kept.','\n') } dij = dij2 }
    }

    if (fdis) {
        # calculate fdis
        FDis <- FD::fdisp(d = dij, a = as.matrix(comm))$FDis
    }

    comm <- as.matrix(comm)
    N <- nrow(comm)
    S <- ncol(comm)
    SR <- rowSums(comm > 0)  # species richness of each site
    comm <- sweep(comm, 1, rowSums(comm, na.rm = TRUE), "/")  # relative abun

    dij <- as.matrix(dij)
    if(all(is.na(dij[upper.tri(dij)])))
        stop("All pairwise distance is NA, do all species have the same trait values?")
    if (stand_dij)
        dij <- dij/max(dij)

    # inter = comm %*% dij # \sum_i,j_S(p_i * dij) Q = rowSums(sweep(comm,1,inter,'*',
    # check.margin = F))/2 # \sum_j_S\sum_i,j_S(p_i * dij)
    Q <- vector("numeric", length = N)
    names(Q) <- dimnames(comm)[[1]]
    for (k in 1:N) {
        Q[k] <- (comm[k, ]) %*% dij %*% matrix(comm[k, ], ncol = 1)  # /2
        # Q[k] = sum(dij * outer(comm[k,], comm[k,], '*'))/2
    }

    ## D_q
    FD_q <- MD_q <- D_q <- vector("numeric", length = N)
    names(D_q) <- dimnames(comm)[[1]]
    names(MD_q) <- dimnames(comm)[[1]]
    names(FD_q) <- dimnames(comm)[[1]]

    if (q == 0) {
        for (k in 1:N) {
            df2 <- comm[k, ][comm[k, ] > 0]
            dis2 <- dij[names(df2), names(df2)]
            if (Q[k] == 0) {
                D_q[k] <- 0
            } else {
                D_q[k] <- sum(dis2/Q[k])^0.5
            }
            MD_q[k] <- D_q[k] * Q[k]
            FD_q[k] <- (D_q[k])^2 * Q[k]
        }
    } else {
        if (q == 1) {
            for (k in 1:N) {
                df2 <- comm[k, ][comm[k, ] > 0]
                dis2 <- dij[names(df2), names(df2)]
                if (Q[k] == 0) {
                    D_q[k] <- 0
                } else {
                    D_q[k] <- exp(-0.5 * sum(dis2/Q[k] * outer(df2, df2, FUN = "*") * log(outer(df2,
                                                                                                df2, FUN = "*"), base)))
                    # exp(-0.5 * (dij/Q) * pi*pj * log(pi*pj)
                }
                MD_q[k] <- D_q[k] * Q[k]
                FD_q[k] <- (D_q[k])^2 * Q[k]
            }
        } else {
            # q != 0 or 1
            for (k in 1:N) {
                df2 <- comm[k, ][comm[k, ] > 0]
                dis2 <- dij[names(df2), names(df2)]
                if (Q[k] == 0) {
                    D_q[k] <- 0
                } else {
                    din <- sum(dis2/Q[k] * (outer(df2, df2, FUN = "*")^q))
                    if (din == 0) {
                        D_q[k] <- 0
                    } else {
                        D_q[k] <- din^(1/(2 * (1 - q)))
                    }
                }
                MD_q[k] <- D_q[k] * Q[k]
                FD_q[k] <- (D_q[k])^2 * Q[k]
            }
        }
    }

    if (fdis) {
        if (div_by_sp == TRUE) {
            return(rbind(Q, FDis, D_q/SR, MD_q/SR, FD_q/choose(SR, 2)))
        } else {
            return(rbind(Q, FDis, D_q, MD_q, FD_q))
        }
    } else {
        if (div_by_sp == TRUE) {
            return(rbind(Q, D_q/SR, MD_q/SR, FD_q/choose(SR, 2)))
        } else {
            return(rbind(Q, D_q, MD_q, FD_q))
        }
    }
}
