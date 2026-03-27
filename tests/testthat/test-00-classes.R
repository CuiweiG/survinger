test_that("package loads without error", {
  expect_true(requireNamespace("survinger", quietly = TRUE))
})
