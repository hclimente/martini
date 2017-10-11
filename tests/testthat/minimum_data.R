gwas <- list()
gwas$map <- read.table(text = "
                       chr snp.names cm gpos allele.1 allele.2
                       1 rs1 0 10 A G
                       1 rs2 0 20 A G
                       1 rs3 0 30 A G
                       2 rs4 0 15 A G
                       2 rs5 0 25 A G
                       2 rs6 0 35 A G
                       ", header = TRUE, stringsAsFactors = FALSE)

snpMapping <- read.table(text = "
                       snp gene
                       rs1 A
                       rs2 A
                       rs5 B
                       rs6 B
                       ", header = TRUE, stringsAsFactors = FALSE)
ppi <- read.table(text = "
                  gene1 gene2
                  A B
                  ", header = TRUE, stringsAsFactors = FALSE)

gs <- get_GS_network(gwas)
gm <- get_GM_network(gwas, snpMapping = snpMapping)
gi <- get_GI_network(gwas, snpMapping = snpMapping, ppi = ppi)

result <- read.table(text = "
                       chr snp.names cm gpos allele.1 allele.2 selected C module
                       1 rs1 0 10 A G TRUE 100 1
                       1 rs2 0 20 A G FALSE 10 NA
                       1 rs3 0 30 A G FALSE 10 NA
                       2 rs4 0 15 A G FALSE 10 NA
                       2 rs5 0 25 A G FALSE 10 NA
                       2 rs6 0 35 A G TRUE 100 1
                       ", header = TRUE, stringsAsFactors = FALSE)