test_that("runAE checks sampleID", {
  x <- matrix(
    c(10, 5, 8, 2),
    nrow = 2,
    dimnames = list(c("gene1", "gene2"), c("s1", "s2"))
  )

  expect_error(runAE(x, "s3", filterByFPKM = FALSE), "sampleID")
})
