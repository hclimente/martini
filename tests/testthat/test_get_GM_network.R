library(martini)
data(examplegwas)

test_that("we interconnect snps from a gene", {
  expect_true(are.connected(get_GM_network(gwas, data.frame(snp = c("rs111", "rs120"), gene = "1")), "rs111", "rs120"))
  expect_true(are.connected(get_GM_network(gwas, 
                                           data.frame(snp = c("rs111", "rs120", "rs164"), 
                                                      gene = "1")),
                            "rs111", "rs164"))
  expect_false(are.connected(get_GM_network(gwas, 
                                            data.frame(snp = c("rs111", "rs120", "rs164"), 
                                                       gene = c("1","1","2"))),
                             "rs111", "rs164"))
})

test_that("crash if snp2gene is insufficient to create a GM network", {
  expect_error(get_GM_network(gwas, data.frame(snp = "rs111", gene = "A")),
               "the data frame should contain at least two columns")
})