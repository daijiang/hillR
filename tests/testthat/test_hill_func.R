context("testing functional diveristy")
dummy = FD::dummy

test_that("FDis, FD vs hillR", {
  a = hillR::hill_func(comm = dummy$abun, traits = dummy$trait, q = 0)
  b = FD::fdisp(FD::gowdis(dummy$trait), dummy$abun)$FDis
  expect_equal(a["FDis", ], b)
})

test_that("functional similarity should be between 0 and 1; q = 0", {
  a = hillR::hill_func_parti(comm = dummy$abun, traits = dummy$trait, q = 0)
  expect_lte(a$local_similarity, 1)
  expect_gte(a$local_similarity, 0)
  expect_lte(a$region_similarity, 1)
  expect_gte(a$region_similarity, 0)
})

test_that("functional similarity should be between 0 and 1; q = 1", {
  a = hillR::hill_func_parti(comm = dummy$abun, traits = dummy$trait, q = 1)
  expect_lte(a$local_similarity, 1)
  expect_gte(a$local_similarity, 0)
  expect_lte(a$region_similarity, 1)
  expect_gte(a$region_similarity, 0)
})

test_that("functional similarity should be between 0 and 1; q = 2", {
  a = hillR::hill_func_parti(comm = dummy$abun, traits = dummy$trait, q = 2)
  expect_lte(a$local_similarity, 1)
  expect_gte(a$local_similarity, 0)
  expect_lte(a$region_similarity, 1)
  expect_gte(a$region_similarity, 0)
})

test_that("pairwise similarity", {
  a = hillR::hill_func_parti_pairwise(comm = dummy$abun, traits = dummy$trait, q = 0)
  expect_equal(nrow(a), choose(nrow(dummy$abun), 2))

  a2 = hillR::hill_func_parti_pairwise(comm = dummy$abun, traits = dummy$trait, q = 1)
  expect_equal(nrow(a2), choose(nrow(dummy$abun), 2))
})
