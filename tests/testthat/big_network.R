gwas <- list()
gwas$map <- read.table(text = "
                       chr snp.names cm gpos allele.1 allele.2
                       1 1A1 0 10 A G
                       1 1A2 0 20 A G
                       1 1A3 0 30 A G
                       1 1A4 0 40 A G
                       1 1A5 0 50 A G
                       1 1A6 0 60 A G
                       1 1-1 0 65 A G
                       1 1B1 0 70 A G
                       1 1B2 0 80 A G
                       1 1B3 0 90 A G
                       1 1B4 0 100 A G
                       1 1B5 0 110 A G
                       1 1B6 0 115 A G
                       1 1-2 0 120 A G
                       2 2C1 0 35 A G
                       2 2C2 0 45 A G
                       2 2C3 0 55 A G
                       2 2C4 0 65 A G
                       2 2C5 0 75 A G
                       2 2C6 0 80 A G
                       2 2-1 0 85 A G
                       2 2-2 0 95 A G
                       2 2D1 0 105 A G
                       2 2D2 0 115 A G
                       2 2D3 0 85 A G
                       ", header = TRUE, stringsAsFactors = FALSE)

snpMapping <- data.frame(snp = gwas$map$snp.names,
                         gene = substr(gwas$map$snp.names, 2, 2))
snpMapping <- subset(snpMapping, gene != "-")
ppi <- read.table(text = "
                  gene1 gene2
                  A B
                  A C
                  B D
                  ", header = TRUE, stringsAsFactors = FALSE)

gs <- get_GS_network(gwas)
gm <- get_GM_network(gwas, snpMapping = snpMapping)
gi <- get_GI_network(gwas, snpMapping = snpMapping, ppi = ppi)