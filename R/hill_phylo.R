#' code for hill_phylo and hill_phylo_parti are mostly from Chiu & Chao.
#'
#' @param comm data frame of vegtation data. Sites as rows, species as columns.
#' @param tree a phylogeny with class 'phylo'.
#' @author Chiu & Chao
#' @noRd
dat_prep_phylo <- function(comm, tree) {
    if (class(tree) == "phylo")
        tree <- ape::write.tree(tree)
    phyloData <- ade4::newick2phylog(tree)
    comm <- as.matrix(comm[, names(phyloData$leaves)])  # resort sp
    nodenames <- c(names(phyloData$leaves), names(phyloData$nodes))

    M <- matrix(0, nrow = ncol(comm), ncol = length(nodenames), dimnames = list(names(phyloData$leaves),
        nodenames))

    for (i in 1:nrow(M)) {
        M[i, ][unlist(phyloData$paths[i])] <- 1
    }

    phylo_comm <- comm %*% M
    phyloLength <- c(phyloData$leaves, phyloData$nodes)
    treeH <- sum(phyloLength * phylo_comm[1, ]/sum(comm[1, ]))

    return(list(pcomm = t(phylo_comm), pLength = phyloLength, treeH = treeH))
}

#' Phylogenetic diversity through Hill Numbers
#'
#' Calculate alpha phylogenetic diversity based on Hill numbers
#'
#' @inheritParams hill_taxa
#' @param tree a phylogeny with class 'phylo'.
#' @inheritParams hill_taxa_parti
#' @author Chiu & Chao
#' @return A vector of hill number based phylogenetic diversity for all sites.
#' @export
#' @references Chao, Anne, Chun-Huo Chiu, and Lou Jost. Unifying Species Diversity, Phylogenetic Diversity, Functional Diversity, and Related Similarity and Differentiation Measures Through Hill Numbers. Annual Review of Ecology, Evolution, and Systematics 45, no. 1 (2014): 297–324. <doi:10.1146/annurev-ecolsys-120213-091540>.
#' @examples
#' comm = dummy = FD::dummy$abun
#' tree = ape::rtree(n = ncol(comm), tip.label = paste0('sp', 1:8))
#' hill_phylo(comm, tree, q = 0)
#' hill_phylo(comm, tree, q = 0.999)
#' hill_phylo(comm, tree, q = 1)
#' hill_phylo(comm, tree, q = 2)
#'
hill_phylo <- function(comm, tree, q = 0, base = exp(1), rel_then_pool = TRUE, show.warning = TRUE) {
    if (any(comm < 0))
        stop("Negative value in comm data")
    # if(any(colSums(comm) == 0) & show.warning) warning('Some species in comm data were
    # not observed in any site,\n delete them...') comm = comm[, colSums(comm) != 0]

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
    PD <- numeric(N)
    names(PD) <- row.names(comm)

    if (q == 1) {
        for (i in 1:N) {
            TT <- sum(pabun[, i] * plength)
            I <- which(pabun[, i] > 0)
            PD[i] <- exp(-sum(plength[I] * (pabun[, i][I]/TT) * log(pabun[, i][I]/TT,
                base)))
        }
    } else {
        for (i in 1:N) {
            TT <- sum(pabun[, i] * plength)
            I <- which(pabun[, i] > 0)
            PD[i] <- sum(plength[I] * (pabun[, i][I]/TT)^q)^(1/(1 - q))
        }
    }

    PD
}
