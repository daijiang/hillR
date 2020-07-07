#' Phylogenetic diversity of multiple sites
#'
#' Calculate overall phylogenetic diversity and site similarity across multiple sites.
#'
#' @inheritParams hill_phylo
#' @param phy_abund A matrix of phylogeny node and tips by community matrix derived
#' from `dat_prep_phylo()`. Can be specified to speed up `hill_phylo_parti_pairwise()`.
#' @param check_data Whether to check the community data and phylogeny. Default is `TRUE`.
#' Can be set to `FALSE` to speed up `hill_phylo_parti_pairwise()`.
#' @export
#' @author Chiu & Chao, Daijiang Li
#' @references Chao, Anne, Chun-Huo Chiu, and Lou Jost. Unifying Species Diversity, Phylogenetic Diversity, Functional Diversity, and Related Similarity and Differentiation Measures Through Hill Numbers. Annual Review of Ecology, Evolution, and Systematics 45, no. 1 (2014): 297–324. <doi:10.1146/annurev-ecolsys-120213-091540>.
#' @return A data frame with one row (across all sites) and six columns: q, gamma diversity, alpha diversity,
#' beta diversity, local similarity (similar to Sorensen), and region similarity (similar to Jaccard).
#' @examples
#' comm = dummy = FD::dummy$abun
#' tree = ape::rtree(n = ncol(comm), tip.label = paste0('sp', 1:8))
#' hill_phylo_parti(comm, tree, q = 0)
#' hill_phylo_parti(comm, tree, q = 0.999)
#' hill_phylo_parti(comm, tree, q = 1)
#' hill_phylo_parti(comm, tree, q = 2)
#'
hill_phylo_parti <- function(comm, tree, q = 0, base = exp(1), rel_then_pool = TRUE,
                             show_warning = TRUE, phy_abund = NULL, check_data = TRUE) {
    if(check_data){
        if (any(comm < 0))
            stop("Negative value in comm data")
        # if(any(colSums(comm) == 0) & show_warning) warning('Some species in comm data were
        # not observed in any site,\n delete them...') comm = comm[, colSums(comm) != 0] #

        comm_sp <- intersect(colnames(comm), tree$tip.label)

        if (class(tree) != "phylo")
            stop("tree must be an object with phylo as class")
        if (length(setdiff(tree$tip.label, comm_sp))) {
            if (show_warning)
                warning("Some species in the phylogeny but not in comm, \n remove them from the phylogeny...")
            tree <- ape::keep.tip(tree, comm_sp)
        }

        if (length(setdiff(colnames(comm), comm_sp))) {
            if (show_warning)
                warning("Some species in the comm but not in the phylogeny, \n remove them from the comm")
            comm <- comm[, comm_sp]
        }

        comm <- comm[, tree$tip.label]  # resort sp
        comm <- as.matrix(comm)

        if (rel_then_pool) {
            comm <- sweep(comm, 1, rowSums(comm, na.rm = TRUE), "/")  # relative abun
        }
    }

    if(is.null(phy_abund)){
        pabun <- dat_prep_phylo(comm, tree)
    } else{ # already calculated
        pabun <- phy_abund[, rownames(comm)]
    }

    plength <- tree$edge.length

    N <- ncol(pabun)
    Tabun <- rowSums(pabun)
    stopifnot(length(Tabun) == length(plength))
    gT <- sum(Tabun * plength)
    gI <- which(Tabun > 0)

    if (q == 1) {
        gPD <- exp(-sum(plength[gI] * (Tabun[gI]/gT) * log(Tabun[gI]/gT, base)))
        L <- matrix(rep(plength, N), ncol = N)
        aI <- which(pabun > 0)
        aPD <- exp(-sum(L[aI] * (pabun[aI]/gT) * log(pabun[aI]/gT, base)))/N
        bPD <- gPD/aPD
        phyloCqN <- 1 - log(bPD, base)/log(N, base)
        phyloUqN <- phyloCqN
    } else {
        gPD <- sum(plength[gI] * (Tabun[gI]/gT)^q)^(1/(1 - q))
        L <- matrix(rep(plength, N), ncol = N)
        aI <- which(pabun > 0)
        aPD <- sum(L[aI] * (pabun[aI]/gT)^q)^(1/(1 - q))/N
        bPD <- gPD/aPD
        phyloCqN <- 1 - (bPD^(1 - q) - 1)/(N^(1 - q) - 1)
        phyloUqN <- 1 - (bPD^(q - 1) - 1)/(N^(q - 1) - 1)
    }

    return(data.frame(q = q, PD_gamma = gPD, PD_alpha = aPD, PD_beta = bPD, local_similarity = phyloCqN,
                      region_similarity = phyloUqN))
}
