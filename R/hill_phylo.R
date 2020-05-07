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
        tree <- ape::keep.tip(tree, comm_sp)
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

    N <- nrow(comm)
    PD <- numeric(N)
    names(PD) <- row.names(comm)

    br_L <- tree$edge.length
    node_tips <- lapply(tree$edge[, 2], function(nd) geiger::tips(tree, nd))

    for (i in 1:N) {
        # for each node, which of its descends are in this site? how abundant in total?
        a_i <- sapply(node_tips, function(tips_per_node) sum(comm[i, ][tips_per_node]))
        TT <- sum(br_L * a_i) # weight by branch length
        # remove nodes with zeros
        I <- which(a_i > 0)
        br_i = br_L[I]
        a_i = a_i[I]
        if(q == 1){
            PD[i] <- exp(-sum(br_i * (a_i/TT) * log(a_i/TT, base)))
        } else {
            PD[i] <- sum(br_i * (a_i/TT)^q)^(1/(1 - q))
        }
    }

    PD
}
