% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hill_func_parti.R
\name{hill_func_parti}
\alias{hill_func_parti}
\title{Decompostion of functional diversity through Hill Numbers}
\usage{
hill_func_parti(
  comm,
  traits,
  traits_as_is = FALSE,
  q = 0,
  base = exp(1),
  check_data = TRUE,
  rel_then_pool = TRUE,
  ord = c("podani", "metric"),
  stand_dij = FALSE,
  show_warning = TRUE
)
}
\arguments{
\item{comm}{A data frame of vegetation data. Sites as rows, species as columns.}

\item{traits}{A data frame of species functional traits data. Species as rows, traits as columns.
It can include both continuous and categorical data. It will be transformed into a distance
matrix using `FD::gowdis(traits)`. If all traits are numeric, then it will use Euclidean distance.}

\item{traits_as_is}{if \code{FALSE} (default) traits data frame will be transformed into a distance
matrix. Otherwise, will use as is (i.e. traits is a symmetric distance matrix).}

\item{q}{Hill number, \code{q} = 0 (default) to get species richness,
\code{q} = 1 to get shannon entropy, \code{q} = 2 will give inverse Simpson.}

\item{base}{default is \code{exp(1)}, the base of log.}

\item{check_data}{whether to check data first? Default is \code{TRUE}.}

\item{rel_then_pool}{default is \code{TRUE.} Abundance of species are first changed to relative abundance within sites,
then pooled into one assemblage. If \code{FALSE}, sites are pooled first, then change abundance of species
to relative abundance.}

\item{ord}{ord in \code{FD::gowdis}.}

\item{stand_dij}{whether to standardize distance matrix to have max value of 1? Default is \code{FALSE}.}

\item{show_warning}{whether to print warning, default is \code{TRUE}.}
}
\value{
a data frame with one row (across all sites), including these columns: q, RaoQ of pooled assemblage,
gamma diversity, alpha diversity, beta diversity, local species overlap (similar to Sorensen), and region species
overlap (similar to Jaccard). See Chiu and Chao 2014 Table 3 for more information.
}
\description{
Calculate functional gamma, alpha, and beta diversity for all communities, as
 well as site similarity. These values are based on ALL communities.
}
\examples{
dummy = FD::dummy
hill_func_parti(comm = dummy$abun, traits = dummy$trait, q = 0)
hill_func_parti(comm = dummy$abun, traits = dummy$trait, q = 1)
hill_func_parti(comm = dummy$abun, traits = dummy$trait, q = 0.9999)
hill_func_parti(comm = dummy$abun, traits = dummy$trait, q = 2)
hill_func_parti(comm = dummy$abun, traits = dummy$trait, q = 3)

}
\references{
Chao, Anne, Chun-Huo Chiu, and Lou Jost. Unifying Species Diversity, Phylogenetic Diversity, Functional Diversity, and Related Similarity and Differentiation Measures Through Hill Numbers. Annual Review of Ecology, Evolution, and Systematics 45, no. 1 (2014): 297–324. <doi:10.1146/annurev-ecolsys-120213-091540>.

Chiu, Chun-Huo, and Anne Chao. Distance-Based Functional Diversity Measures and Their Decomposition: A Framework Based on Hill Numbers. PLoS ONE 9, no. 7 (July 7, 2014): e100014. <doi:10.1371/journal.pone.0100014>.
}
\seealso{
\code{\link{hill_taxa_parti}}, \code{\link{hill_func}}
}
