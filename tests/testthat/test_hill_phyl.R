context("testing phylogenetic diversity")
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
  a = hillR::hill_phylo_parti_pairwise(comm = dummy$abun, tree, q = 0, .progress = FALSE)
  expect_equal(nrow(a), choose(nrow(dummy$abun), 2))

  a2 = hillR::hill_phylo_parti_pairwise(comm = dummy$abun, tree, q = 1, .progress = FALSE)
  expect_equal(nrow(a2), choose(nrow(dummy$abun), 2))
})

test_that("when N = 2, hill_phylo_parti equals Sorensen", {
  toy.comm = matrix(nrow = 2, ncol = 6)
  rownames(toy.comm) = c("A","B")
  colnames(toy.comm) = c("sp1","sp2","sp3","sp4","sp5","sp6")
  toy.comm[1,] = c(1,1,1,0,0,0)
  toy.comm[2,] = c(0,0,1,1,1,1)

  toy.tree = ape::read.tree(text="(((sp1:1,sp2:1):5,(sp3:3,sp4:3):3):2,(sp5:7,sp6:7):1);")

  skip_if_not_installed("betapart")

  toy.phylobetamulti = betapart::phylo.beta.multi(toy.comm, toy.tree, index.family = "sor")
  toy.phylobetamulti2 = betapart::phylo.beta.multi(toy.comm, toy.tree, index.family = "jac")
  toy.hill = hillR::hill_phylo_parti(toy.comm, toy.tree)
  # local_similarity equals with phylo Sorensen
  expect_equal(1 - toy.phylobetamulti$phylo.beta.SOR, toy.hill$local_similarity)
  # regional_similarity equals with phylo Jaccard (1 - UniFrac)
  expect_equal(1 - toy.phylobetamulti2$phylo.beta.JAC, toy.hill$region_similarity)
})
