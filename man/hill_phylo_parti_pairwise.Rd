% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hill_phylo_parti_pairwise.R
\name{hill_phylo_parti_pairwise}
\alias{hill_phylo_parti_pairwise}
\title{Pairwise phylogenetic diversity through Hill numbers}
\usage{
hill_phylo_parti_pairwise(
  comm,
  tree,
  q = 0,
  output = c("data.frame", "matrix"),
  pairs = c("unique", "full"),
  rel_then_pool = TRUE,
  .progress = TRUE,
  show_warning = TRUE,
  ...
)
}
\arguments{
\item{comm}{A data frame of vegetation data. Sites as rows, species as columns.}

\item{tree}{A phylogeny with class 'phylo'.}

\item{q}{Hill number, \code{q} = 0 (default) to get species richness,
\code{q} = 1 to get shannon entropy, \code{q} = 2 will give inverse Simpson.}

\item{output}{output type: data.frame (default) or matrix. If matrix, then this function will return a list of matrices.}

\item{pairs}{full or unique (default). Do you want to compare all possible pairs (i.e. n^2) or just unique pairs (i.e. \code{choose(n, 2))}?}

\item{rel_then_pool}{default is \code{TRUE.} Abundance of species are first changed to relative abundance within sites,
then pooled into one assemblage. If \code{FALSE}, sites are pooled first, then change abundance of species
to relative abundance.}

\item{.progress}{Whether to show progress bar. Default is `TRUE`.}

\item{show_warning}{whether to print warning, default is \code{TRUE}.}

\item{...}{additional arguments for \code{hill_func_parti}.}
}
\value{
A data frame or a matrix with results for all pairwise comparisons.
}
\description{
Calculate pairwise phylogenetic diversity.
}
\examples{
\dontrun{
comm = dummy = FD::dummy$abun
tree = ape::rtree(n = ncol(comm), tip.label = paste0('sp', 1:8))
hill_phylo_parti_pairwise(comm, tree, q = 0, show_warning = FALSE)
hill_phylo_parti_pairwise(comm, tree, q = 0.999, show_warning = FALSE)
hill_phylo_parti_pairwise(comm, tree, q = 1, show_warning = FALSE)
hill_phylo_parti_pairwise(comm, tree, q = 2, show_warning = FALSE)
}
}
\references{
Chao, Anne, Chun-Huo Chiu, and Lou Jost. Unifying Species Diversity, Phylogenetic Diversity, Functional Diversity, and Related Similarity and Differentiation Measures Through Hill Numbers. Annual Review of Ecology, Evolution, and Systematics 45, no. 1 (2014): 297–324. <doi:10.1146/annurev-ecolsys-120213-091540>.
}
\seealso{
\code{\link{hill_phylo_parti}}
}
