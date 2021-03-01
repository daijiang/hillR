context("comparing taxa diversity")
dummy = FD::dummy

test_that("vegetarian package vs hillR, taxa alpha diversity, q = 0", {

  skip_if_not_installed("vegetarian")

  a = apply(dummy$abun, 1, vegetarian::d, lev = "alpha", q = 0)
  b = hillR::hill_taxa(dummy$abun, q = 0)
  expect_equal(a, b)
})

test_that("vegetarian package vs hillR, taxa alpha diversity, q = 0.9999", {

  skip_if_not_installed("vegetarian")

  a = apply(dummy$abun, 1, vegetarian::d, lev = "alpha", q = 0.9999)
  b = hillR::hill_taxa(dummy$abun, q = 0.9999)
  expect_equal(a, b)
})

test_that("vegetarian package vs hillR, taxa alpha diversity, q = 1", {

  skip_if_not_installed("vegetarian")

  a = apply(dummy$abun, 1, vegetarian::d, lev = "alpha", q = 1)
  b = hillR::hill_taxa(dummy$abun, q = 1)
  expect_equal(a, b)
})

test_that("vegetarian package vs hillR, taxa alpha diversity, q = 2", {

  skip_if_not_installed("vegetarian")

  a = apply(dummy$abun, 1, vegetarian::d, lev = "alpha", q = 2)
  b = hillR::hill_taxa(dummy$abun, q = 2)
  expect_equal(a, b)
})

test_that("vegetarian package vs hillR, taxa beta diversity, q = 0", {

  skip_if_not_installed("vegetarian")

  a = vegetarian::d(dummy$abun, lev = "beta", q = 0)
  b = hillR::hill_taxa_parti(dummy$abun, q = 0)
  expect_equal(a, b$TD_beta)

  a2 = vegetarian::M.homog(dummy$abun, q = 0)
  expect_equal(a2, b$M_homog)

  a3 = vegetarian::similarity(dummy$abun, q = 0)
  expect_equal(a3, b$local_similarity)
})

test_that("vegetarian package vs hillR, taxa beta diversity, q = 0.9999", {

  skip_if_not_installed("vegetarian")

  a = vegetarian::d(dummy$abun, lev = "beta", q = 0.9999)
  b = hillR::hill_taxa_parti(dummy$abun, q = 0.9999)
  expect_equal(a, b$TD_beta)

  a2 = vegetarian::M.homog(dummy$abun, q = 0.9999)
  expect_equal(a2, b$M_homog)

  a3 = vegetarian::similarity(dummy$abun, q = 0.9999)
  expect_equal(a3, b$local_similarity)
})

test_that("vegetarian package vs hillR, taxa beta diversity, q = 1", {

  skip_if_not_installed("vegetarian")

  a = vegetarian::d(dummy$abun, lev = "beta", q = 1)
  b = hillR::hill_taxa_parti(dummy$abun, q = 1)
  expect_equal(a, b$TD_beta)

  a2 = vegetarian::M.homog(dummy$abun, q = 1)
  expect_equal(a2, b$M_homog)

  a3 = vegetarian::similarity(dummy$abun, q = 1)
  expect_equal(a3, b$local_similarity)
})

test_that("vegetarian package vs hillR, taxa beta diversity, q = 2", {

  skip_if_not_installed("vegetarian")

  a = vegetarian::d(dummy$abun, lev = "beta", q = 2)
  b = hillR::hill_taxa_parti(dummy$abun, q = 2)
  expect_equal(a, b$TD_beta)

  a2 = vegetarian::M.homog(dummy$abun, q = 2)
  expect_equal(a2, b$M_homog)

  a3 = vegetarian::similarity(dummy$abun, q = 2)
  expect_equal(a3, b$local_similarity)
})

test_that("vegetarian package vs hillR, taxa beta diversity, q = 2, two sites", {

  skip_if_not_installed("vegetarian")

  a = vegetarian::d(dummy$abun[1:2,], lev = "beta", q = 2)
  b = hillR::hill_taxa_parti(dummy$abun[1:2,], q = 2)
  expect_equal(a, b$TD_beta)

  a2 = vegetarian::M.homog(dummy$abun[1:2,], q = 2)
  expect_equal(a2, b$M_homog)

  a3 = vegetarian::similarity(dummy$abun[1:2,], q = 2)
  expect_equal(a3, b$local_similarity)
})

test_that("taxanomic similarity should be between 0 and 1; q = 0", {
  a = hillR::hill_taxa_parti(comm = dummy$abun, q = 0)
  expect_lte(a$local_similarity, 1)
  expect_gte(a$local_similarity, 0)
  expect_lte(a$region_similarity, 1)
  expect_gte(a$region_similarity, 0)
})

test_that("taxanomic similarity should be between 0 and 1; q = 1", {
  a = hillR::hill_taxa_parti(comm = dummy$abun, q = 1)
  expect_lte(a$local_similarity, 1)
  expect_gte(a$local_similarity, 0)
  expect_lte(a$region_similarity, 1)
  expect_gte(a$region_similarity, 0)
})

test_that("taxanomic similarity should be between 0 and 1; q = 2", {
  a = hillR::hill_taxa_parti(comm = dummy$abun, q = 2)
  expect_lte(a$local_similarity, 1)
  expect_gte(a$local_similarity, 0)
  expect_lte(a$region_similarity, 1)
  expect_gte(a$region_similarity, 0)
})

test_that("pairwise similarity", {
  a = hillR::hill_taxa_parti_pairwise(comm = dummy$abun, q = 0, .progress = FALSE)
  expect_equal(nrow(a), choose(nrow(dummy$abun), 2))

  a2 = hillR::hill_taxa_parti_pairwise(comm = dummy$abun, q = 1, .progress = FALSE)
  expect_equal(nrow(a2), choose(nrow(dummy$abun), 2))
})

test_that("when N = 2, hill_taxa_parti equals Sorensen", {

  skip_if_not_installed("betapart")

  toy.comm = matrix(nrow = 2, ncol = 6)
  rownames(toy.comm) = c("A","B")
  colnames(toy.comm) = c("sp1","sp2","sp3","sp4","sp5","sp6")
  toy.comm[1,] = c(1,1,1,0,0,0)
  toy.comm[2,] = c(0,0,1,1,1,1)

  toy.betamulti = betapart::beta.multi(toy.comm, index.family = "sor")
  toy.betamulti2 = betapart::beta.multi(toy.comm, index.family = "jac")
  toy.hill = hillR::hill_taxa_parti(toy.comm)
  # local_similarity equals with Sorensen
  expect_equal(1 - toy.betamulti$beta.SOR, toy.hill$local_similarity)
  # regional_similarity equals with Jaccard
  expect_equal(1 - toy.betamulti2$beta.JAC, toy.hill$region_similarity)
})
