load("examplegwas.rda")

gwas <- examplegwas[['gwas']]
ld <- Matrix::sparseMatrix(i = c(1,1,2), j = c(1,3,3), x = c(0, 1, 0.5), symmetric = T)
colnames(ld) <- c("rs1101", "rs1102", "rs1103")
rownames(ld) <- colnames(ld)

snps <- data.frame(snp1 = c("rs1101","rs1101","rs1102"),
                   snp2 = c("rs1102","rs1103","rs1103"))
net <- igraph::graph_from_data_frame(snps, directed = FALSE)
net <- igraph::set_edge_attr(net, "weight", value = 1)

# subtraction tests
ldnet <- ldweight_edges(net, ld, method = "subtraction")
ldnetDf <- igraph::as_data_frame(ldnet)

test_that("snps in perfect ld are disconnected", {
  expect_false(igraph::are_adjacent(ldnet, "rs1101", "rs1103"))
  expect_equal(length(igraph::E(ldnet)), 2)
})

test_that("snps not in perfect ld connected", {
  expect_true(all(igraph::E(ldnet) %in% igraph::E(net)))
  expect_equal(subset(ldnetDf, from == "rs1101" & to == "rs1102")$weight, 1)
  expect_equal(subset(ldnetDf, from == "rs1102" & to == "rs1103")$weight, 0.5)
})

test_that("negative and NAs are not possible", {
  
  ld <- Matrix::sparseMatrix(i = c(1,1,2), j = c(1,3,3), x = c(0, 2, 0.5), symmetric = TRUE)
  colnames(ld) <- c("rs1101", "rs1102", "rs1103")
  rownames(ld) <- colnames(ld)
  expect_error(ldweight_edges(net, ld, method = "subtraction"),
               "Edge weights cannot be negative.")
  
  ld <- Matrix::sparseMatrix(i = c(1,1,2), j = c(1,3,3), x = c(NA, NA, NA), symmetric = TRUE)
  colnames(ld) <- c("rs1101", "rs1102", "rs1103")
  rownames(ld) <- colnames(ld)
  expect_error(ldweight_edges(net, ld, method = "subtraction"),
               "NA values as edge weights.")
  
})

# inverse tests
ldnet <- ldweight_edges(net, ld, method = "inverse")
ldnetDf <- igraph::as_data_frame(ldnet)

test_that("edges have the apropriate values", {
  expect_equal(length(igraph::E(ldnet)), 3)
  expect_true(all(igraph::E(net) %in% igraph::E(ldnet)))
  expect_equal(subset(ldnetDf, from == "rs1101" & to == "rs1102")$weight, 1)
  expect_equal(subset(ldnetDf, from == "rs1101" & to == "rs1103")$weight, 0.5)
  expect_equal(subset(ldnetDf, from == "rs1102" & to == "rs1103")$weight, 2/3)
})

test_that("negative and NAs are not possible", {
  
  ld <- Matrix::sparseMatrix(i = c(1,1,2), j = c(1,3,3), x = c(0, -10, 0.5), symmetric = T)
  colnames(ld) <- c("rs1101", "rs1102", "rs1103")
  rownames(ld) <- colnames(ld)
  expect_error(ldweight_edges(net, ld, method = "inverse"),
               "Edge weights cannot be negative.")
  
  ld <- Matrix::sparseMatrix(i = c(1,1,2), j = c(1,3,3), x = c(NA, NA, NA), symmetric = T)
  colnames(ld) <- c("rs1101", "rs1102", "rs1103")
  rownames(ld) <- colnames(ld)
  expect_error(ldweight_edges(net, ld, method = "inverse"),
               "NA values as edge weights.")
  
})

# sigmoid tests
ldnet <- ldweight_edges(net, ld, method = "sigmoid")
ldnetDf <- igraph::as_data_frame(ldnet)

test_that("edges have the apropriate values", {
  expect_equal(length(igraph::E(ldnet)), 3)
  expect_true(all(igraph::E(net) %in% igraph::E(ldnet)))
  expect_equal(subset(ldnetDf, from == "rs1101" & to == "rs1102")$weight, 1)
  expect_equal(subset(ldnetDf, from == "rs1101" & to == "rs1103")$weight, 1 / (1 + exp(10*(1 - 0.5))))
  expect_equal(subset(ldnetDf, from == "rs1102" & to == "rs1103")$weight, 1 / (1 + exp(10*(0.5 - 0.5))))
})