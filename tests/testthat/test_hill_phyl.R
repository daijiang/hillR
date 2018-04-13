context("testing phylogenetic diveristy")
dummy = FD::dummy
n = ncol(dummy$abun)
tree = ape::rtree(n, tip.label = paste0("sp", 1:n))

test_that("PD of each site", {
  a = hillR::hill_phylo(comm = dummy$abun, tree, q = 0)
})

test_that("phylogenetic similarity should be between 0 and 1; q = 0", {
  a = hillR::hill_phylo_parti(comm = dummy$abun, tree, q = 0)
  expect_lte(a$local_similarity, 1)
  expect_gte(a$local_similarity, 0)
  expect_lte(a$region_similarity, 1)
  expect_gte(a$region_similarity, 0)
})

test_that("phylogenetic similarity should be between 0 and 1; q = 1", {
  a = hillR::hill_phylo_parti(comm = dummy$abun, tree, q = 1)
  expect_lte(a$local_similarity, 1)
  expect_gte(a$local_similarity, 0)
  expect_lte(a$region_similarity, 1)
  expect_gte(a$region_similarity, 0)
})

test_that("phylogenetic similarity should be between 0 and 1; q = 2", {
  a = hillR::hill_phylo_parti(comm = dummy$abun, tree, q = 2)
  expect_lte(a$local_similarity, 1)
  expect_gte(a$local_similarity, 0)
  expect_lte(a$region_similarity, 1)
  expect_gte(a$region_similarity, 0)
})

test_that("pairwise similarity", {
  a = hillR::hill_phylo_parti_pairwise(comm = dummy$abun, tree, q = 0)
  expect_equal(nrow(a), choose(nrow(dummy$abun), 2))

  a2 = hillR::hill_phylo_parti_pairwise(comm = dummy$abun, tree, q = 1)
  expect_equal(nrow(a2), choose(nrow(dummy$abun), 2))
})
