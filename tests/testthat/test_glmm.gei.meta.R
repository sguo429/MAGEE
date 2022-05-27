context("single-variant GEI and joint test meta-analysis")

test_that("single-variant GEI and Joint test meta-analysis", {
  infile1 <- system.file("extdata", "meta1.txt", package = "MAGEE")
  infile2 <- system.file("extdata", "meta2.txt", package = "MAGEE")
  infile3 <- system.file("extdata", "meta3.txt", package = "MAGEE")
  infile4 <- system.file("extdata", "meta4.txt", package = "MAGEE")
  infile5 <- system.file("extdata", "meta5.txt", package = "MAGEE")
  outfile <- tempfile()
  glmm.gei.meta(files = c(infile1, infile2, infile3, infile4, infile5), outfile = outfile, interaction="sex")
  out <- read.table(outfile, header = TRUE, as.is = TRUE)
  outPval<-c(7.398859e-01, 1.057351e-03, 4.710037e-01, 9.201342e-01, 6.607810e-01,7.933586e-02, 4.072482e-01, 1.997179e-01, 8.607651e-01, 2.508775e-03,9.637767e-01, 2.540691e-01, 9.556030e-01, 6.061708e-01, 6.718866e-01,8.917068e-01, 5.476881e-01, 6.643518e-01, 4.905400e-01, 9.769166e-01,2.244303e-01, 9.112266e-02, 3.718101e-05, 8.142838e-01, 3.259739e-01,4.175679e-01, 5.000037e-01, 3.301131e-01, 6.235216e-01, 8.682540e-01,2.557019e-01, 1.277121e-01, 7.233825e-06, 3.025667e-01, 6.202095e-01,4.754977e-01, 9.279986e-02, 9.792363e-01)
  expect_equal(signif(out$P_Value_Marginal, 5), signif(outPval, 5))
  unlink(outfile)
})

