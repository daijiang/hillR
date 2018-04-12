#' Pairwise comparisons for all sites
#'
#' Calculate pairwise taxonomic gamma, alpha, and beta diversity for communities, as
#' well as site similarity. It is based on \code{\link{hill_taxa_parti}}.
#' If comm has >2 sites, this function will give results for all pairwise comparisons.
#'
#' @author Daijiang Li
#'
#' @param comm data frame of vegtation data. Sites as rows, species as columns.
#' @param q hill number, q = 0 (default) to get species richness,
#' q = 1 to get shannon entropy, q = 2 will give inverse Simpson.
#' @param rel_then_pool default is TRUE. Abundance of species are first changed to relative abundance within sites,
#'  then pooled into one assemblage. If FALSE, sites are pooled first, then change abundance of species
#'  to relative abundance.
#' @param output output type: data.frame (default) or matrix. If matrix, then this function will return a list of matrices.
#' @param pairs full or unique (default). Do you want to compare all possible pairs (i.e. n^2) or just unique pairs (i.e. choose(n, 2))?
#' @export
#' @return a data frame with results for all pairwise comparisons.
#' @seealso \code{\link{hill_taxa_parti}}
#' @examples
#' dummy = FD::dummy
#' hill_taxa_parti_pairwise(comm = dummy$abun, q = 0)
#' hill_taxa_parti_pairwise(comm = dummy$abun, q = 0, output = "matrix")
#' hill_taxa_parti_pairwise(comm = dummy$abun, q = 1)
#' hill_taxa_parti_pairwise(comm = dummy$abun, q = 0.9999999)
#' hill_taxa_parti_pairwise(comm = dummy$abun, q = 0.9999999, rel_then_pool = F)
#' hill_taxa_parti_pairwise(comm = dummy$abun, q = 1, rel_then_pool = F)
#' hill_taxa_parti_pairwise(comm = dummy$abun, q = 2)
#' hill_taxa_parti_pairwise(comm = dummy$abun, q = 3)
#'
#'
hill_taxa_parti_pairwise = function(comm, q = 0, rel_then_pool = TRUE,
                                    output = c("data.frame", "matrix"),
                                    pairs = c( "unique", "full")){
  output <- match.arg(output)
  pairs <- match.arg(pairs)
  nsite = nrow(comm)
  temp = matrix(1, nsite, nsite)
  dimnames(temp) = list(row.names(comm), row.names(comm))
  gamma_pair = alpha_pair = beta_pair = local_simi = region_simi = temp
  for(i in 1:nsite){
    for(j in i:nsite){
      o = hill_taxa_parti(comm[c(i,j), ], q = q)
      gamma_pair[i,j] = o$TD_gamma; gamma_pair[j,i] = o$TD_gamma
      alpha_pair[i,j] = o$TD_alpha; alpha_pair[j,i] = o$TD_alpha
      beta_pair[i,j] = o$TD_beta; beta_pair[j,i] = o$TD_beta
      local_simi[i,j] = o$local_taxa_overlap; local_simi[j,i] = o$local_taxa_overlap
      region_simi[i,j] = o$region_taxa_overlap; region_simi[j,i] = o$region_taxa_overlap
    }
  }

  if(pairs == "full"){
    if(output == "matrix"){
      out = list(q = q, TD_gamma = gamma_pair, TD_alpha = alpha_pair, TD_beta = beta_pair,
                 local_taxa_overlap = local_simi, region_taxa_overlap = region_simi)
    }

    if(output == "data.frame"){
      site.comp = as.matrix(expand.grid(row.names(comm), row.names(comm)))
      out = plyr::adply(site.comp, 1, function(x){
        data.frame(q = q,
                   site1 = x[1],
                   site2 = x[2],
                   TD_gamma = gamma_pair[x[1], x[2]],
                   TD_alpha = alpha_pair[x[1], x[2]],
                   TD_beta = beta_pair[x[1], x[2]],
                   local_taxa_overlap = local_simi[x[1], x[2]],
                   region_taxa_overlap = region_simi[x[1], x[2]])
      })[, -1] # get rid of X1 column
    }
  }

  if(pairs == "unique"){
    gamma_pair[lower.tri(gamma_pair, diag = TRUE)] = NA
    alpha_pair[lower.tri(alpha_pair, diag = TRUE)] = NA
    beta_pair[lower.tri(beta_pair, diag = TRUE)] = NA
    local_simi[lower.tri(local_simi, diag = TRUE)] = NA
    region_simi[lower.tri(region_simi, diag = TRUE)] = NA

    if(output == "matrix"){
      out = list(q = q, TD_gamma = gamma_pair, TD_alpha = alpha_pair, TD_beta = beta_pair,
                 local_taxa_overlap = local_simi, region_taxa_overlap = region_simi)
    }

    if(output == "data.frame"){
      site.comp = as.matrix(expand.grid(row.names(comm), row.names(comm)))
      out = plyr::adply(site.comp, 1, function(x){
        data.frame(q = q,
                   site1 = x[1],
                   site2 = x[2],
                   TD_gamma = gamma_pair[x[1], x[2]],
                   TD_alpha = alpha_pair[x[1], x[2]],
                   TD_beta = beta_pair[x[1], x[2]],
                   local_taxa_overlap = local_simi[x[1], x[2]],
                   region_taxa_overlap = region_simi[x[1], x[2]])
      })
      out = na.omit(out)[, -1]
      row.names(out) = NULL
    }
  }
out
}

