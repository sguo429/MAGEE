context("Mixed model Association tests for GEne-Environment interactions (MAGEE)")

test_that("400 gaussian", {
  
  gdsfile  <- system.file("extdata", "geno.gds",  package = "MAGEE")
  bgenfile <- system.file("extdata", "geno.bgen", package = "MAGEE")
  samplefile <- system.file("extdata", "geno.sample", package = "MAGEE")
  group.file <- system.file("extdata", "SetID.withweights.txt", package = "MAGEE")
  data(example)
  suppressWarnings(RNGversion("3.5.0"))
  set.seed(123)
  pheno <- rbind(example$pheno, example$pheno[1:100, ])
  pheno$id <- 1:500
  pheno$disease[sample(1:500,20)] <- NA
  pheno$age[sample(1:500,20)] <- NA
  pheno$sex[sample(1:500,20)] <- NA
  pheno <- pheno[sample(1:500,450), ]
  pheno <- pheno[pheno$id <= 400, ]
  kins <- example$GRM
  
  ### single thread with kins
  obj1 <- glmmkin(trait ~ age + sex, data = pheno, kins = kins, id = "id", family = gaussian(link = "identity"))
  out1 <- MAGEE(null.obj=obj1, interaction="sex", geno.file=gdsfile, group.file=group.file, tests = c("JV", "JF", "JD"), use.minor.allele = T, auto.flip = T)
  out1_bgen <- MAGEE(null.obj=obj1, interaction="sex", geno.file=bgenfile, bgen.samplefile = samplefile, group.file=group.file, tests = c("JV", "JF", "JD"), use.minor.allele = T, auto.flip = T)
  
  expect_equal(signif(range(out1$IV.pval)), signif(c(0.0442781, 0.9028424)))
  expect_equal(signif(range(out1$IF.pval)), signif(c(0.03147219, 0.97079474)))
  expect_equal(signif(range(out1$JV.pval)), signif(c(0.06884383, 0.82438351)))
  expect_equal(signif(range(out1$JF.pval)), signif(c(0.07428697, 0.90053289)))
  expect_equal(signif(range(out1$JD.pval)), signif(c(0.07417719, 0.85422359)))
  expect_equal(out1, out1_bgen)
  
  ### single thread without kins
  obj2 <- glmmkin(trait ~ age + sex, data = pheno, kins = NULL, id = "id", family = gaussian(link = "identity"))
  out2 <- MAGEE(null.obj=obj2, interaction="sex", geno.file=gdsfile, group.file=group.file, tests = c("JV", "JF", "JD"), use.minor.allele = T, auto.flip = T)
  out2_bgen <- MAGEE(null.obj=obj2, interaction="sex", geno.file=bgenfile, bgen.samplefile = samplefile, group.file=group.file, tests = c("JV", "JF", "JD"), use.minor.allele = T, auto.flip = T)
  expect_equal(signif(range(out2$IV.pval)), signif(c(0.04641134, 0.89956347)))
  expect_equal(signif(range(out2$IF.pval)), signif(c(0.03844619, 0.94288432)))
  #  expect_equal(signif(range(out2$JV.pval)), signif(c(0.05562695, 0.63240910)))
  #  expect_equal(signif(range(out2$JF.pval)), signif(c(0.08496645, 0.79924790)))
  expect_equal(signif(range(out2$JD.pval)), signif(c(0.08749729, 0.81941241)))
  expect_equal(out2, out2_bgen)
  
  ### multi-thread
  skip_on_cran()
  out1.tmp <- MAGEE(null.obj=obj1, interaction="sex", geno.file=gdsfile, group.file=group.file, tests = c("JV", "JF", "JD"), use.minor.allele = T, auto.flip = T, ncores = 2)
  expect_equal(out1, out1.tmp)
  out1_bgen.tmp <- MAGEE(null.obj=obj1, interaction="sex", geno.file=bgenfile, bgen.samplefile = samplefile, group.file=group.file, tests = c("JV", "JF", "JD"), use.minor.allele = T, auto.flip = T, ncores = 2)
  expect_equal(out1, out1_bgen.tmp)
  out2.tmp <- MAGEE(null.obj=obj2, interaction="sex", geno.file=gdsfile, group.file=group.file, tests = c("JV", "JF", "JD"), use.minor.allele = T, auto.flip = T, ncores = 2)
  expect_equal(out2, out2.tmp)
  out2_bgen.tmp <- MAGEE(null.obj=obj2, interaction="sex", geno.file=bgenfile, bgen.samplefile = samplefile, group.file=group.file, tests = c("JV", "JF", "JD"), use.minor.allele = T, auto.flip = T, ncores = 2)
  expect_equal(out2, out2_bgen.tmp)
  
  ### re-order id
  idx <- sample(nrow(pheno))
  pheno <- pheno[idx, ]
  obj1 <- glmmkin(trait ~ age + sex, data = pheno, kins = kins, id = "id", family = gaussian(link = "identity"))
  tmpout <- MAGEE(null.obj=obj1, interaction="sex", geno.file=gdsfile,  group.file=group.file, tests = c("JV", "JF", "JD"), use.minor.allele = T, auto.flip = T)
  expect_equal(out1, tmpout)
  tmpout <- MAGEE(null.obj=obj1, interaction="sex", geno.file=bgenfile, bgen.samplefile = samplefile, group.file=group.file, tests = c("JV", "JF", "JD"), use.minor.allele = T, auto.flip = T)
  expect_equal(out1_bgen, tmpout)
  obj2 <- glmmkin(trait ~ age + sex, data = pheno, kins = NULL, id = "id", family = gaussian(link = "identity"))
  tmpout <- MAGEE(null.obj=obj2, interaction="sex", geno.file=gdsfile, group.file=group.file, tests = c("JV", "JF", "JD"), use.minor.allele = T, auto.flip = T)
  expect_equal(out2, tmpout)
  tmpout <- MAGEE(null.obj=obj2, interaction="sex", geno.file=bgenfile, bgen.samplefile = samplefile, group.file=group.file, tests = c("JV", "JF", "JD"), use.minor.allele = T, auto.flip = T)
  expect_equal(out2_bgen, tmpout)
  
  ### re-order id
  idx <- sample(nrow(kins))
  kins <- kins[idx, idx]
  obj1 <- glmmkin(trait ~ age + sex, data = pheno, kins = kins, id = "id", family = gaussian(link = "identity"))
  tmpout <- MAGEE(null.obj=obj1, interaction="sex", geno.file=gdsfile, group.file=group.file, tests = c("JV", "JF", "JD"), use.minor.allele = T, auto.flip = T)
  expect_equal(out1, tmpout)
  tmpout <- MAGEE(null.obj=obj1, interaction="sex", geno.file=bgenfile, bgen.samplefile = samplefile, group.file=group.file, tests = c("JV", "JF", "JD"), use.minor.allele = T, auto.flip = T)
  expect_equal(out1_bgen, tmpout)
  obj2 <- glmmkin(trait ~ age + sex, data = pheno, kins = NULL, id = "id", family = gaussian(link = "identity"))
  tmpout <- MAGEE(null.obj=obj2, interaction="sex", geno.file=gdsfile, group.file=group.file, tests = c("JV", "JF", "JD"), use.minor.allele = T, auto.flip = T)
  expect_equal(out2, tmpout)
  tmpout <- MAGEE(null.obj=obj2, interaction="sex", geno.file=bgenfile, bgen.samplefile = samplefile, group.file=group.file, tests = c("JV", "JF", "JD"), use.minor.allele = T, auto.flip = T)
  expect_equal(out2_bgen, tmpout)
})

test_that("400 binomial", {
  skip_on_cran()
  
  gdsfile <- system.file("extdata", "geno.gds", package = "MAGEE")
  group.file <- system.file("extdata", "SetID.withweights.txt", package = "MAGEE")
  data(example)
  suppressWarnings(RNGversion("3.5.0"))
  set.seed(123)
  pheno <- rbind(example$pheno, example$pheno[1:100, ])
  pheno$id <- 1:500
  pheno$disease[sample(1:500,20)] <- NA
  pheno$age[sample(1:500,20)] <- NA
  pheno$sex[sample(1:500,20)] <- NA
  pheno <- pheno[sample(1:500,450), ]
  pheno <- pheno[pheno$id <= 400, ]
  kins <- example$GRM
  
  ### single thread with kins
  obj1 <- glmmkin(disease ~ age + sex, data = pheno, kins = kins, id = "id", family = binomial(link = "logit"))
  out1 <- MAGEE(null.obj=obj1, interaction="sex", geno.file=gdsfile, group.file=group.file, tests = c("JV", "JF", "JD"), use.minor.allele = T, auto.flip = T)
  out1_bgen <- MAGEE(null.obj=obj1, interaction="sex", geno.file=bgenfile, bgen.samplefile = samplefile, group.file=group.file, tests = c("JV", "JF", "JD"), use.minor.allele = T, auto.flip = T)
  expect_equal(signif(range(out1$IV.pval)), signif(c(0.0124653, 0.9469504)))
  expect_equal(signif(range(out1$IF.pval)), signif(c(0.009173469, 0.827455722)))
  expect_equal(signif(range(out1$JV.pval)), signif(c(0.01985053, 0.89813390)))
  expect_equal(signif(range(out1$JF.pval)), signif(c(0.03580556, 0.87867728)))
  expect_equal(signif(range(out1$JD.pval)), signif(c(0.03194091, 0.89022491)))
  expect_equal(out1, out1_bgen)
  
  ### single thread without kins
  obj2 <- glmmkin(disease ~ age + sex, data = pheno, kins = NULL, id = "id", family = binomial(link = "logit"))
  out2 <- MAGEE(null.obj=obj2, interaction="sex", geno.file=gdsfile, group.file=group.file, tests = c("JV", "JF", "JD"), use.minor.allele = T, auto.flip = T)
  out2_bgen <- MAGEE(null.obj=obj2, interaction="sex", geno.file=bgenfile, bgen.samplefile = samplefile, group.file=group.file, tests = c("JV", "JF", "JD"), use.minor.allele = T, auto.flip = T)
  expect_equal(signif(range(out2$IV.pval)), signif(c(0.01105658, 0.94118125)))
  expect_equal(signif(range(out2$IF.pval)), signif(c(0.007766316, 0.809971207)))
  expect_equal(signif(range(out2$JV.pval)), signif(c(0.02021408, 0.90127556)))
  expect_equal(signif(range(out2$JF.pval)), signif(c(0.02663287, 0.81920025)))
  expect_equal(signif(range(out2$JD.pval)), signif(c(0.02450435, 0.83817777)))
  expect_equal(out2, out2_bgen)
  
  ### multi-thread
  skip_on_cran()
  out1.tmp <- MAGEE(null.obj=obj1, interaction="sex", geno.file=gdsfile, group.file=group.file, tests = c("JV", "JF", "JD"), use.minor.allele = T, auto.flip = T, ncores = 2)
  expect_equal(out1, out1.tmp)
  out1_bgen.tmp <- MAGEE(null.obj=obj1, interaction="sex", geno.file=bgenfile, bgen.samplefile = samplefile, group.file=group.file, tests = c("JV", "JF", "JD"), use.minor.allele = T, auto.flip = T, ncores = 2)
  expect_equal(out1, out1_bgen.tmp)
  out2.tmp <- MAGEE(null.obj=obj2, interaction="sex", geno.file=gdsfile, group.file=group.file, tests = c("JV", "JF", "JD"), use.minor.allele = T, auto.flip = T, ncores = 2)
  expect_equal(out2, out2.tmp)
  out2_bgen.tmp <- MAGEE(null.obj=obj2, interaction="sex", geno.file=bgenfile, bgen.samplefile = samplefile, group.file=group.file, tests = c("JV", "JF", "JD"), use.minor.allele = T, auto.flip = T, ncores = 2)
  expect_equal(out2, out2_bgen.tmp)
  
  ### re-order id
  idx <- sample(nrow(pheno))
  pheno <- pheno[idx, ]
  obj1 <- glmmkin(disease ~ age + sex, data = pheno, kins = kins, id = "id", family = binomial(link = "logit"))
  tmpout <- MAGEE(null.obj=obj1, interaction="sex", geno.file=gdsfile, group.file=group.file, tests = c("JV", "JF", "JD"), use.minor.allele = T, auto.flip = T)
  expect_equal(out1, tmpout)
  tmpout <- MAGEE(null.obj=obj1, interaction="sex", geno.file=bgenfile, bgen.samplefile = samplefile, group.file=group.file, tests = c("JV", "JF", "JD"), use.minor.allele = T, auto.flip = T)
  expect_equal(out1_bgen, tmpout)
  obj2 <- glmmkin(disease ~ age + sex, data = pheno, kins = NULL, id = "id", family = binomial(link = "logit"))
  tmpout <- MAGEE(null.obj=obj2, interaction="sex", geno.file=gdsfile, group.file=group.file, tests = c("JV", "JF", "JD"), use.minor.allele = T, auto.flip = T)
  expect_equal(out2, tmpout)
  tmpout <- MAGEE(null.obj=obj2, interaction="sex", geno.file=bgenfile, bgen.samplefile = samplefile, group.file=group.file, tests = c("JV", "JF", "JD"), use.minor.allele = T, auto.flip = T)
  expect_equal(out2_bgen, tmpout)
  
  ### re-order id
  idx <- sample(nrow(kins))
  kins <- kins[idx, idx]
  obj1 <- glmmkin(disease ~ age + sex, data = pheno, kins = kins, id = "id", family = binomial(link = "logit"))
  tmpout <- MAGEE(null.obj=obj1, interaction="sex", geno.file=gdsfile, group.file=group.file, tests = c("JV", "JF", "JD"), use.minor.allele = T, auto.flip = T)
  expect_equal(out1, tmpout)
  tmpout <- MAGEE(null.obj=obj1, interaction="sex", geno.file=bgenfile, bgen.samplefile = samplefile, group.file=group.file, tests = c("JV", "JF", "JD"), use.minor.allele = T, auto.flip = T)
  expect_equal(out1_bgen, tmpout)
  obj2 <- glmmkin(disease ~ age + sex, data = pheno, kins = NULL, id = "id", family = binomial(link = "logit"))
  tmpout <- MAGEE(null.obj=obj2, interaction="sex", geno.file=gdsfile, group.file=group.file, tests = c("JV", "JF", "JD"), use.minor.allele = T, auto.flip = T)
  expect_equal(out2, tmpout)
  tmpout <- MAGEE(null.obj=obj2, interaction="sex", geno.file=bgenfile, bgen.samplefile = samplefile, group.file=group.file, tests = c("JV", "JF", "JD"), use.minor.allele = T, auto.flip = T)
  expect_equal(out2_bgen, tmpout)
})

### multi-phenotype MAGEE
test_that("multiple phenotypes gaussian", {
  skip_on_cran()
  
  gdsfile <- system.file("extdata", "geno.gds", package = "MAGEE")
  bgenfile <- system.file("extdata", "geno.bgen", package = "MAGEE")
  samplefile <- system.file("extdata", "geno.sample", package = "MAGEE")
  group.file <- system.file("extdata", "SetID.withweights.txt", package = "MAGEE")
  data(example)
  suppressWarnings(RNGversion("3.5.0"))
  set.seed(103)
  kins <- example$GRM
  tau1 <- matrix(c(3,0.5,0,0.5,2.5,-0.1,0,-0.1,3),3,3)
  tau2 <- matrix(c(2.5,0.8,0.2,0.8,4.8,-0.1,0.2,-0.1,2.8),3,3)
  kins.chol <- chol(tau1 %x% kins + tau2 %x% diag(400))
  tmp <- as.vector(crossprod(kins.chol, rnorm(1200)))
  x1 <- rnorm(400)
  x2 <- rbinom(400,1,0.5)
  pheno <- data.frame(id = 1:400, x1 = x1, x2 = x2, y1 = 0.5*x1+0.8*x2+tmp[1:400], y2 = x1-0.3*x2+tmp[401:800], y3 = x2+tmp[801:1200])
  ### single thread with kins
  obj1 <- glmmkin(cbind(y1,y2,y3)~x1+x2, data = pheno, kins = kins, id = "id", family = gaussian(link = "identity"))
  out1 <- MAGEE(null.obj=obj1, interaction="x1", geno.file=gdsfile, group.file=group.file, tests = c("JV", "JF", "JD"), use.minor.allele = T, auto.flip = T)
  out1_bgen <- MAGEE(null.obj=obj1, interaction="x1", geno.file=bgenfile, bgen.samplefile = samplefile, group.file=group.file, tests = c("JV", "JF", "JD"), use.minor.allele = T, auto.flip = T)
  
  expect_equal(signif(range(out1$IV.pval)), signif(c(0.1060731, 0.9129642)))
  expect_equal(signif(range(out1$IF.pval)), signif(c( 0.1785101, 0.9041499)))
  expect_equal(signif(range(out1$JV.pval)), signif(c(0.261298, 0.940976)))
  expect_equal(signif(range(out1$JF.pval)), signif(c(0.2591923, 0.9611525)))
  expect_equal(signif(range(out1$JD.pval)), signif(c(0.2641673, 0.9661435)))
  expect_equal(out1, out1_bgen)
  
  ### single thread without kins
  obj2 <- glmmkin(cbind(y1,y2,y3)~x1+x2, data = pheno, kins = NULL, id = "id", family = gaussian(link = "identity"))
  out2 <- MAGEE(null.obj=obj2, interaction="x1", geno.file=gdsfile, group.file=group.file, tests = c("JV", "JF", "JD"), use.minor.allele = T, auto.flip = T)
  out2_bgen <- MAGEE(null.obj=obj2, interaction="x1", geno.file=bgenfile, bgen.samplefile = samplefile, group.file=group.file, tests = c("JV", "JF", "JD"), use.minor.allele = T, auto.flip = T)
  
  expect_equal(signif(range(out2$IV.pval)), signif(c(0.1830878, 0.9157741)))
  expect_equal(signif(range(out2$IF.pval)), signif(c(0.2745413, 0.8634707)))
  expect_equal(signif(range(out2$JV.pval)), signif(c(0.003123721, 0.803773857)))
  expect_equal(signif(range(out2$JF.pval)), signif(c( 0.004704365, 0.737300921)))
  expect_equal(signif(range(out2$JD.pval)), signif(c(0.002961206, 0.761153895)))
  expect_equal(out2, out2_bgen)
  
  ### multi-thread
  skip_on_cran()
  out1.tmp <- MAGEE(null.obj=obj1, interaction="x1", geno.file=gdsfile, group.file=group.file, tests = c("JV", "JF", "JD"), use.minor.allele = T, auto.flip = T, ncores = 2)
  expect_equal(out1, out1.tmp)
  out1_bgen.tmp <- MAGEE(null.obj=obj1, interaction="x1", geno.file=bgenfile, bgen.samplefile = samplefile, group.file=group.file, tests = c("JV", "JF", "JD"), use.minor.allele = T, auto.flip = T, ncores = 2)
  expect_equal(out1, out1_bgen.tmp)
  out2.tmp <- MAGEE(null.obj=obj2, interaction="x1", geno.file=gdsfile, group.file=group.file, tests = c("JV", "JF", "JD"), use.minor.allele = T, auto.flip = T, ncores = 2)
  expect_equal(out2, out2.tmp)
  out2_bgen.tmp <- MAGEE(null.obj=obj2, interaction="x1", geno.file=bgenfile, bgen.samplefile = samplefile, group.file=group.file, tests = c("JV", "JF", "JD"), use.minor.allele = T, auto.flip = T, ncores = 2)
  expect_equal(out2, out2_bgen.tmp)
  
  ### re-order id
  idx <- sample(nrow(pheno))
  pheno <- pheno[idx, ]
  obj1 <- glmmkin(cbind(y1,y2,y3)~x1+x2, data = pheno, kins = kins, id = "id", family = gaussian(link = "identity"))
  tmpout <- MAGEE(null.obj=obj1, interaction="x1", geno.file=gdsfile, group.file=group.file, tests = c("JV", "JF", "JD"), use.minor.allele = T, auto.flip = T)
  expect_equal(out1, tmpout)
  obj2 <- glmmkin(cbind(y1,y2,y3)~x1+x2, data = pheno, kins = NULL, id = "id", family = gaussian(link = "identity"))
  tmpout <- MAGEE(null.obj=obj2, interaction="x1", geno.file=gdsfile, group.file=group.file, tests = c("JV", "JF", "JD"), use.minor.allele = T, auto.flip = T)
  expect_equal(out2, tmpout)
  
  ### re-order id
  idx <- sample(nrow(kins))
  kins <- kins[idx, idx]
  obj1 <- glmmkin(cbind(y1,y2,y3)~x1+x2, data = pheno, kins = kins, id = "id", family = gaussian(link = "identity"))
  tmpout <- MAGEE(null.obj=obj1, interaction="x1", geno.file=gdsfile, group.file=group.file, tests = c("JV", "JF", "JD"), use.minor.allele = T, auto.flip = T)
  expect_equal(out1, tmpout)
  obj2 <- glmmkin(cbind(y1,y2,y3)~x1+x2, data = pheno, kins = NULL, id = "id", family = gaussian(link = "identity"))
  tmpout <- MAGEE(null.obj=obj2, interaction="x1", geno.file=gdsfile, group.file=group.file, tests = c("JV", "JF", "JD"), use.minor.allele = T, auto.flip = T)
  expect_equal(out2, tmpout)
})

### lonngitudinal MAGEE
test_that("longitudinal random time trend gaussian", {
  skip_on_cran()
  
  gdsfile  <- system.file("extdata", "geno.gds",  package = "MAGEE")
  bgenfile <- system.file("extdata", "geno.bgen", package = "MAGEE")
  samplefile <- system.file("extdata", "geno.sample", package = "MAGEE")
  group.file <- system.file("extdata", "SetID.withweights.txt", package = "MAGEE")
  data(example)
  suppressWarnings(RNGversion("3.5.0"))
  set.seed(123)
  pheno <- example$pheno2
  kins <- example$GRM
  ### single thread with kins
  obj1 <- glmmkin(y.trend ~ sex + time, data = pheno, kins = kins, id = "id",random.slope = "time", family = gaussian(link = "identity"))
  out1 <- MAGEE(null.obj=obj1, interaction="sex", geno.file=gdsfile, group.file=group.file, tests = c("JV", "JF", "JD"), time= "time", use.minor.allele = T, auto.flip = T)
  out1_bgen <- MAGEE(null.obj=obj1, interaction="sex", geno.file=bgenfile, bgen.samplefile = samplefile, group.file=group.file, tests = c("JV", "JF", "JD"), time= "time", use.minor.allele = T, auto.flip = T)
  
  expect_equal(signif(range(out1$IV.pval)), signif(c(0.00772688, 0.96278208)))
  expect_equal(signif(range(out1$IF.pval)), signif(c(0.02210924, 0.99445043)))
  expect_equal(signif(range(out1$JV.pval)), signif(c(0.02190347, 0.97984894)))
  expect_equal(signif(range(out1$JF.pval)), signif(c(0.08805881, 0.99536133)))
  expect_equal(signif(range(out1$JD.pval)), signif(c(0.07763854, 0.99368482)))
  expect_equal(out1, out1_bgen)
  
  ### single thread without kins
  obj2 <- glmmkin(y.trend ~ sex + time, data = pheno, kins = NULL, id = "id",random.slope = "time", family = gaussian(link = "identity"))
  out2 <- MAGEE(null.obj=obj2, interaction="sex", geno.file=gdsfile, group.file=group.file, tests = c("JV", "JF", "JD"), time= "time", use.minor.allele = T, auto.flip = T)
  out2_bgen <- MAGEE(null.obj=obj2, interaction="sex", geno.file=bgenfile, bgen.samplefile = samplefile, group.file=group.file, tests = c("JV", "JF", "JD"), time= "time", use.minor.allele = T, auto.flip = T)
  
  expect_equal(signif(range(out2$IV.pval)), signif(c(0.009518904, 0.972481750)))
  expect_equal(signif(range(out2$IF.pval)), signif(c(0.02777976, 0.99456861)))
  expect_equal(signif(range(out2$JV.pval)), signif(c(0.01523866, 0.99165479)))
  expect_equal(signif(range(out2$JF.pval)), signif(c(0.07038413, 0.99273235)))
  expect_equal(signif(range(out2$JD.pval)), signif(c(0.06933245, 0.99219566)))
  expect_equal(out2, out2_bgen)
  
  ### re-order id
  idx <- sample(nrow(pheno))
  pheno <- pheno[idx, ]
  obj1 <- glmmkin(y.trend ~ sex + time, data = pheno, kins = kins, id = "id",random.slope = "time", family = gaussian(link = "identity"))
  tmpout <- MAGEE(null.obj=obj1, interaction="sex", geno.file=gdsfile, group.file=group.file, tests = c("JV", "JF", "JD"), time= "time", use.minor.allele = T, auto.flip = T)
  expect_equal(out1, tmpout)
  obj2 <- glmmkin(y.trend ~ sex + time, data = pheno, kins = NULL, id = "id",random.slope = "time", family = gaussian(link = "identity"))
  tmpout <- MAGEE(null.obj=obj2, interaction="sex", geno.file=gdsfile, group.file=group.file, tests = c("JV", "JF", "JD"), time= "time", use.minor.allele = T, auto.flip = T)
  expect_equal(out2, tmpout)
  
  ### re-order id
  idx <- sample(nrow(kins))
  kins <- kins[idx, idx]
  obj1 <- glmmkin(y.trend ~ sex + time, data = pheno, kins = kins, id = "id",random.slope = "time", family = gaussian(link = "identity"))
  tmpout <- MAGEE(null.obj=obj1, interaction="sex", geno.file=gdsfile, group.file=group.file, tests = c("JV", "JF", "JD"), time= "time", use.minor.allele = T, auto.flip = T)
  expect_equal(out1, tmpout)
  obj2 <- glmmkin(y.trend ~ sex + time, data = pheno, kins = NULL, id = "id",random.slope = "time", family = gaussian(link = "identity"))
  tmpout <- MAGEE(null.obj=obj2, interaction="sex", geno.file=gdsfile, group.file=group.file, tests = c("JV", "JF", "JD"), time= "time", use.minor.allele = T, auto.flip = T)
  expect_equal(out2, tmpout)
})
