#' code for hill_phylo and hill_phylo_parti are mostly from Chiu & Chao.
#'
#' @param comm data frame of vegetation data. Sites as rows, species as columns.
#' @param tree a phylogeny with class 'phylo'.
#' @noRd
dat_prep_phylo <- function(comm, tree) {
    if(is.null(tree$edge.length))
        stop("tree must have branch length")
    node_tips <- lapply(tree$edge[, 2], function(nd) geiger::tips(tree, nd))
    n_node = ape::Nnode(tree, internal = FALSE)
    xn = tree$edge[, 2] # put tip and node names there
    xn[xn < min(tree$edge[,1])] = tree$tip.label

    M2 <- matrix(0, nrow = ncol(comm), ncol = n_node - 1, # no root needed
                 dimnames = list(tree$tip.label, xn))
    for (i in 1:ncol(M2)) {
        M2[, i][unlist(node_tips[[i]])] = 1
    }

    t(comm %*% M2)
}

#' Phylogenetic diversity through Hill Numbers
#'
#' Calculate alpha phylogenetic diversity based on Hill numbers
#'
#' @inheritParams hill_taxa
#' @param tree A phylogeny with class 'phylo'.
#' @param return_dt Whether to return the Phylogenetic Hill numbers Dt, default is `FALSE`.
#' @inheritParams hill_taxa_parti
#' @author Chiu & Chao & Daijiang Li
#' @return A vector of hill number based phylogenetic diversity (`PD(T)`, effective total branch length) for all sites. If a site has less than 2 species in the final dataset, NA will be returned.
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
hill_phylo <- function(comm, tree, q = 0, base = exp(1), rel_then_pool = TRUE,
                       show_warning = TRUE, return_dt = FALSE) {
    if (any(comm < 0))
        stop("Negative value in comm data")
    # if(any(colSums(comm) == 0) & show_warning) warning('Some species in comm data were
    # not observed in any site,\n delete them...') comm = comm[, colSums(comm) != 0]

  # remove species with no observations
   if (any(colSums(comm) == 0) & show_warning)
     warning("Some species in comm data were not observed in any site,\n delete them...")
   if (any(rowSums(comm) == 0) & show_warning)
     warning("Some sites in comm data do not have any species,\n delete them...")

   comm <- comm[rowSums(comm) != 0, colSums(comm) != 0, drop = FALSE]

    comm_sp <- intersect(colnames(comm), tree$tip.label)

    if (!inherits(tree, "phylo"))
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

    if(ncol(comm) < 2) stop("The final community data have less than 2 species, which does not make much sense for phylogenetic diversity calculation")

    if(q == 0) comm[comm > 0] = 1

    if (rel_then_pool) {
        comm <- sweep(comm, 1, rowSums(comm, na.rm = TRUE), "/")  # relative abun
    }

    pabun <- dat_prep_phylo(comm, tree)
    plength <- tree$edge.length
    N <- ncol(pabun)
    PD <- D_t <- numeric(N)
    names(PD) <- row.names(comm)
    names(D_t) <- row.names(comm)
    for (i in 1:N) {
      # In cases that a community does not have species from the old clades,
      # when analyze that community separately, the old clades would be trimmed from the tree,
      # leading to an internal node to be the root of that sub-tree. When analyzed separately,
      # the root of that sub-tree, which is an internal node here, will not be calculated.
      # However, in the case of analyzing all communities together, this internal node will
      # be included, with an accumulative relative abundance of 1.
      # To make sure we get the same results whether analyzing all communities together,
      # or analyzing them individually, we set this internal node's accumulative abundance to 0,
      # i.e., not included it in the calculation.
        pabun[, i][which(abs(pabun[, i] - 1) < 0.000001)] = 0

        TT <- sum(pabun[, i] * plength)
        I <- which(pabun[, i] > 0)
        if(length(I) == 0) {
          PD[i] <- NA # only one species
          D_t[i] <- NA
          next()
        }
        D_t[i] <- TT
        if(q == 1){
            PD[i] <- exp(-sum(plength[I] * (pabun[, i][I]/TT) * log(pabun[, i][I]/TT, base)))
        } else {
            PD[i] <- sum(plength[I] * (pabun[, i][I]/TT)^q)^(1/(1 - q))
        }
    }

    if(return_dt){
        D_t = PD/D_t
        return(rbind(D_t, PD))
    } else {
        return(PD)
    }
}

