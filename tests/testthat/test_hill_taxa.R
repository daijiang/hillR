context("comparing taxa diveristy")
dummy = FD::dummy

test_that("vegetariam package vs hillR, taxa alpha diversity, q = 0", {
  a = apply(dummy$abun, 1, vegetarian::d, lev = "alpha", q = 0)
  b = hillR::hill_taxa(dummy$abun, q = 0)
  expect_equal(a, b)
})

test_that("vegetariam package vs hillR, taxa alpha diversity, q = 0.9999", {
  a = apply(dummy$abun, 1, vegetarian::d, lev = "alpha", q = 0.9999)
  b = hillR::hill_taxa(dummy$abun, q = 0.9999)
  expect_equal(a, b)
})

test_that("vegetariam package vs hillR, taxa alpha diversity, q = 1", {
  a = apply(dummy$abun, 1, vegetarian::d, lev = "alpha", q = 1)
  b = hillR::hill_taxa(dummy$abun, q = 1)
  expect_equal(a, b)
})

test_that("vegetariam package vs hillR, taxa alpha diversity, q = 2", {
  a = apply(dummy$abun, 1, vegetarian::d, lev = "alpha", q = 2)
  b = hillR::hill_taxa(dummy$abun, q = 2)
  expect_equal(a, b)
})

test_that("vegetariam package vs hillR, taxa beta diversity, q = 0", {
  a = vegetarian::d(dummy$abun, lev = "beta", q = 0)
  b = hillR::hill_taxa_parti(dummy$abun, q = 0)
  expect_equal(a, b$TD_beta)

  a2 = vegetarian::M.homog(dummy$abun, q = 0)
  expect_equal(a2, b$M_homog)

  a3 = vegetarian::similarity(dummy$abun, q = 0)
  expect_equal(a3, b$local_similarity)
})

test_that("vegetariam package vs hillR, taxa beta diversity, q = 0.9999", {
  a = vegetarian::d(dummy$abun, lev = "beta", q = 0.9999)
  b = hillR::hill_taxa_parti(dummy$abun, q = 0.9999)
  expect_equal(a, b$TD_beta)

  a2 = vegetarian::M.homog(dummy$abun, q = 0.9999)
  expect_equal(a2, b$M_homog)

  a3 = vegetarian::similarity(dummy$abun, q = 0.9999)
  expect_equal(a3, b$local_similarity)
})

test_that("vegetariam package vs hillR, taxa beta diversity, q = 1", {
  a = vegetarian::d(dummy$abun, lev = "beta", q = 1)
  b = hillR::hill_taxa_parti(dummy$abun, q = 1)
  expect_equal(a, b$TD_beta)

  a2 = vegetarian::M.homog(dummy$abun, q = 1)
  expect_equal(a2, b$M_homog)

  a3 = vegetarian::similarity(dummy$abun, q = 1)
  expect_equal(a3, b$local_similarity)
})

test_that("vegetariam package vs hillR, taxa beta diversity, q = 2", {
  a = vegetarian::d(dummy$abun, lev = "beta", q = 2)
  b = hillR::hill_taxa_parti(dummy$abun, q = 2)
  expect_equal(a, b$TD_beta)

  a2 = vegetarian::M.homog(dummy$abun, q = 2)
  expect_equal(a2, b$M_homog)

  a3 = vegetarian::similarity(dummy$abun, q = 2)
  expect_equal(a3, b$local_similarity)
})

test_that("vegetariam package vs hillR, taxa beta diversity, q = 2, two sites", {
  a = vegetarian::d(dummy$abun[1:2,], lev = "beta", q = 2)
  b = hillR::hill_taxa_parti(dummy$abun[1:2,], q = 2)
  expect_equal(a, b$TD_beta)

  a2 = vegetarian::M.homog(dummy$abun[1:2,], q = 2)
  expect_equal(a2, b$M_homog)

  a3 = vegetarian::similarity(dummy$abun[1:2,], q = 2)
  expect_equal(a3, b$local_similarity)
})
