minimum_gwas <- list()
minimum_gwas$map <- read.table(text = "
                       chr snp.name cm gpos allele.1 allele.2
                       1 rs1 0 10 A G
                       1 rs2 0 20 A G
                       1 rs3 0 30 A G
                       2 rs4 0 15 A G
                       2 rs5 0 25 A G
                       2 rs6 0 35 A G
                       2 rs7 0 45 A G
                       2 rs8 0 55 A G
                       2 rs9 0 65 A G
                       2 rs10 0 75 A G
                       2 rs11 0 85 A G
                       ", header = TRUE, stringsAsFactors = FALSE)

snpMapping <- read.table(text = "
                       snp gene
                       rs1 A
                       rs2 A
                       rs5 B
                       rs6 B
                       rs8 C
                       rs9 C
                       rs10 D
                       rs11 C
                       rs11 D
                       ", header = TRUE, stringsAsFactors = FALSE)
ppi <- read.table(text = "
                  gene1 gene2
                  A B
                  A C
                  B D
                  ", header = TRUE, stringsAsFactors = FALSE)

mini_gs <- get_GS_network(minimum_gwas)
mini_gm <- get_GM_network(minimum_gwas, snpMapping = snpMapping)
mini_gi <- get_GI_network(minimum_gwas, snpMapping = snpMapping, ppi = ppi)

result <- read.table(text = "
                       chr snp.name cm gpos allele.1 allele.2 selected C module
                       1 rs1 0 10 A G TRUE 100 1
                       1 rs2 0 20 A G FALSE 10 NA
                       1 rs3 0 30 A G FALSE 10 NA
                       2 rs4 0 15 A G FALSE 10 NA
                       2 rs5 0 25 A G FALSE 10 NA
                       2 rs6 0 35 A G TRUE 100 1
                       ", header = TRUE, stringsAsFactors = FALSE)
