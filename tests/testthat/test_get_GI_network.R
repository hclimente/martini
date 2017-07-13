library(martini)
data(examplegwas)

test_that("we interconnect the right genes", {
  expect_true(
    are.connected(
      get_GI_network(gwas, 
                     data.frame(snp = c("rs111", "rs120", "rs161", "rs162"),
                                      gene = c("A", "A", "B", "B")),
                     data.frame(gene1 = "A", gene2 = "B", stringsAsFactors = FALSE)),
      "rs111", "rs161"))
  expect_false(
    are.connected(
      get_GI_network(gwas, 
                     data.frame(snp = c("rs111", "rs120", "rs161", "rs162"),
                                gene = c("A", "A", "C", "B")),
                     data.frame(gene1 = "A", gene2 = "B", stringsAsFactors = FALSE)),
      "rs111", "rs161"))
})