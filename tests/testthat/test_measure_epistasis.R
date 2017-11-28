library(martini)
data(examplegwas)

subnet <- igraph::delete_edges(examplegwas$net, sample(igraph::E(examplegwas$net), 39564))
epistasis <- measure_epistasis(examplegwas$gwas, subnet)
edges <- igraph::as_data_frame(subnet, what = "edges")
colnames(edges) <- c("snp1", "snp2")
X <- as(examplegwas$gwas$genotypes, "matrix")
colnames(X) <- as.character(examplegwas$gwas$map[,2])

test_that("we calculate epistasis among the expected edges", {
  expect_equal(dim(epistasis), c(10,3))
  expect_equal(dim(merge(epistasis, edges)), c(10,3))
  expect_equal(colnames(epistasis), c("snp1", "snp2", "score"))
})

test_that("epistasis is calculated correctly", {
  expect_equal(subset(epistasis, snp1 == edges$snp1[1] & snp2 == edges$snp2[1])$score, 
               as.numeric( chisq.test(X[,edges$snp1[1]], X[,edges$snp2[1]])$statistic ))
  expect_equal(subset(epistasis, snp1 == edges$snp1[2] & snp2 == edges$snp2[2])$score, 
               as.numeric( chisq.test(X[,edges$snp1[2]], X[,edges$snp2[2]])$statistic ))
  expect_equal(subset(epistasis, snp1 == edges$snp1[3] & snp2 == edges$snp2[3])$score, 
               as.numeric( chisq.test(X[,edges$snp1[3]], X[,edges$snp2[3]])$statistic ))
  expect_equal(subset(epistasis, snp1 == edges$snp1[4] & snp2 == edges$snp2[4])$score, 
               as.numeric( chisq.test(X[,edges$snp1[4]], X[,edges$snp2[4]])$statistic ))
  expect_equal(subset(epistasis, snp1 == edges$snp1[5] & snp2 == edges$snp2[5])$score, 
               as.numeric( chisq.test(X[,edges$snp1[5]], X[,edges$snp2[5]])$statistic ))
  expect_equal(subset(epistasis, snp1 == edges$snp1[6] & snp2 == edges$snp2[6])$score, 
               as.numeric( chisq.test(X[,edges$snp1[6]], X[,edges$snp2[6]])$statistic ))
  expect_equal(subset(epistasis, snp1 == edges$snp1[7] & snp2 == edges$snp2[7])$score, 
               as.numeric( chisq.test(X[,edges$snp1[7]], X[,edges$snp2[7]])$statistic ))
})