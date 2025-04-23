## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, fig.align = 'center', message = FALSE)

## ----loadExampleData----------------------------------------------------------
library(CytoSimplex)
data("rnaRaw")
print(paste0("Class of `rnaRaw`: ", class(rnaRaw), ", dimension of `rnaRaw`: ", nrow(rnaRaw), " genes x ", ncol(rnaRaw), " cells"))

data("rnaCluster")
print(table(rnaCluster))

## ----selectTopGenes-----------------------------------------------------------
vertices <- c("OS", "RE", "CH")
geneSelect <- selectTopFeatures(rnaRaw, clusterVar = rnaCluster, 
                                vertices = vertices, nTop = 30)
head(geneSelect)

stats <- selectTopFeatures(rnaRaw, clusterVar = rnaCluster, 
                           vertices = vertices, nTop = 30, returnStats = TRUE)
head(stats)

## ----plotTernary, fig.height = 4, fig.width = 5-------------------------------
vt.tern <- c("OS", "RE", "CH")
gene.tern <- selectTopFeatures(rnaRaw, clusterVar = rnaCluster, vertices = vt.tern)
plotTernary(rnaRaw, clusterVar = rnaCluster, vertices = vt.tern, features = gene.tern)

## ----plotTernaryInteractive, fig.height = 4, fig.width = 5, eval = FALSE------
# plotTernary(rnaRaw, clusterVar = rnaCluster, vertices = vt.tern, features = gene.tern,
#             interactive = TRUE)

## ----loadVelo-----------------------------------------------------------------
data("rnaVelo")
print(paste0("Class of `rnaVelo`: ", class(rnaVelo), 
             ", dimension of `rnaVelo`: ", nrow(rnaVelo), " x ", ncol(rnaVelo)))

## ----ternVelo, fig.width = 5, fig.height = 4----------------------------------
plotTernary(rnaRaw, clusterVar = rnaCluster, vertices = vt.tern, 
            features = gene.tern, veloGraph = rnaVelo)

## ----splitCluster, fig.width = 7, fig.height = 7------------------------------
library(patchwork)
ternList <- plotTernary(rnaRaw, clusterVar = rnaCluster, vertices = vt.tern, 
                        features = gene.tern, 
                        byCluster = c("Stem", "RE", "ORT", "OS"))
print(names(ternList))
(ternList$Stem + ternList$RE) / (ternList$ORT + ternList$OS)

## ----ternVeloSplit, fig.width = 7, fig.height = 7-----------------------------
veloSplit <- plotTernary(rnaRaw, clusterVar = rnaCluster, vertices = vt.tern, 
                         features = gene.tern, veloGraph = rnaVelo, 
                         byCluster = c("Stem", "RE", "ORT", "OS"))
(veloSplit$Stem + veloSplit$RE) / (veloSplit$ORT + veloSplit$OS)

## ----plotQuaternary, fig.height = 4, fig.width = 6----------------------------
vt.quat <- c("OS", "RE", "CH", "ORT")
gene.quat <- selectTopFeatures(rnaRaw, clusterVar = rnaCluster, vertices = vt.quat)
plotQuaternary(rnaRaw, clusterVar = rnaCluster, vertices = vt.quat, 
               features = gene.quat, veloGraph = rnaVelo, interactive = FALSE)

## ----writeGIF, eval=FALSE, results="hide"-------------------------------------
# writeQuaternaryGIF(rnaRaw, clusterVar = rnaCluster, vertices = vt.quat,
#                    features = gene.quat, veloGraph = rnaVelo,
#                    width = 8, height = 5, res = 200)

