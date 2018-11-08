#' Phylogenetic diversity of multiple sites
#'
#' Calculate overall phylogenetic diversity and site similarity across multiple sites.
#'
#' @inheritParams hill_phylo
#' @export
#' @author Chiu & Chao, Daijiang Li
#' @references Chao, Anne, Chun-Huo Chiu, and Lou Jost. Unifying Species Diversity, Phylogenetic Diversity, Functional Diversity, and Related Similarity and Differentiation Measures Through Hill Numbers. Annual Review of Ecology, Evolution, and Systematics 45, no. 1 (2014): 297–324. <doi:10.1146/annurev-ecolsys-120213-091540>.
#' @return A data frame with one row (across all sites) and six columns: q, gamma diversity, alpha diveristy,
#' beta diversity, local similarity, and region similarity.
#' @examples
#' comm = dummy = FD::dummy$abun
#' tree = ape::rtree(n = ncol(comm), tip.label = paste0('sp', 1:8))
#' hill_phylo_parti(comm, tree, q = 0)
#' hill_phylo_parti(comm, tree, q = 0.999)
#' hill_phylo_parti(comm, tree, q = 1)
#' hill_phylo_parti(comm, tree, q = 2)
#'
hill_phylo_parti <- function(comm, tree, q = 0, base = exp(1), rel_then_pool = TRUE, show.warning = TRUE) {
    if (any(comm < 0))
        stop("Negative value in comm data")
    # if(any(colSums(comm) == 0) & show.warning) warning('Some species in comm data were
    # not observed in any site,\n delete them...') comm = comm[, colSums(comm) != 0] #
    # when only 2 sp, ade4::newick2phylog has trouble, so keep zeros

    comm_sp <- intersect(colnames(comm), tree$tip.label)

    if (class(tree) != "phylo")
        stop("tree must be an object with phylo as class")
    if (length(setdiff(tree$tip.label, comm_sp))) {
        if (show.warning)
            warning("Some species in the phylogeny but not in comm, \n remove them from the phylogeny...")
        tree <- ape::drop.tip(tree, tree$tip.label[!tree$tip.label %in% comm_sp])
    }

    if (length(setdiff(colnames(comm), comm_sp))) {
        if (show.warning)
            warning("Some species in the comm but not in the phylogeny, \n remove them from the comm")
        comm <- comm[, comm_sp]
    }

    comm <- comm[, tree$tip.label]  # resort sp
    comm <- as.matrix(comm)

    if (rel_then_pool) {
        comm <- sweep(comm, 1, rowSums(comm, na.rm = TRUE), "/")  # relative abun
    }

    dat <- dat_prep_phylo(comm, tree)
    pabun <- dat$pcomm
    plength <- dat$pLength

    N <- ncol(pabun)
    Tabun <- rowSums(pabun)
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
