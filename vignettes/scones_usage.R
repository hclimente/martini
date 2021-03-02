## ----load martini-------------------------------------------------------------
library(martini)

## ----print fam----------------------------------------------------------------
head(minigwas$fam)

## ----print map----------------------------------------------------------------
head(minigwas$map)

## ----print genotypes----------------------------------------------------------
minigwas$genotypes

## ----create GS network--------------------------------------------------------
gs <- get_GS_network(minigwas)
par(mar=c(0,0,0,0)+.1)
plot(gs)

## ----create GM network--------------------------------------------------------
gm <- get_GM_network(minigwas, snpMapping = minisnpMapping)
par(mar=c(0,0,0,0)+.1)
plot(gm)

## ----create GI network--------------------------------------------------------
gi <- get_GI_network(minigwas, snpMapping = minisnpMapping, ppi = minippi)
par(mar=c(0,0,0,0)+.1)
plot(gi)

## ----find cones gridsearch----------------------------------------------------
cones <- scones.cv(minigwas, gi)

## -----------------------------------------------------------------------------
head(cones)

