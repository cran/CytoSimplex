## ----setup, include=FALSE-----------------------------------------------------
library(rgl)
knitr::knit_hooks$set(webgl = hook_webgl)
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

## ----plotTernary, fig.height = 4, fig.width = 4-------------------------------
vt.tern <- c("OS", "RE", "CH")
gene.tern <- selectTopFeatures(rnaRaw, clusterVar = rnaCluster, vertices = vt.tern)
plotTernary(rnaRaw, clusterVar = rnaCluster, vertices = vt.tern, features = gene.tern)

## ----loadVelo-----------------------------------------------------------------
data("rnaVelo")
print(paste0("Class of `rnaVelo`: ", class(rnaVelo), 
             ", dimension of `rnaVelo`: ", nrow(rnaVelo), " x ", ncol(rnaVelo)))

## ----ternVelo, fig.width = 4, fig.height = 4----------------------------------
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

## ----plotQuaternary, fig.height = 4, fig.width = 4----------------------------
vt.quat <- c("OS", "RE", "CH", "ORT")
gene.quat <- selectTopFeatures(rnaRaw, clusterVar = rnaCluster, vertices = vt.quat)
plotQuaternary(rnaRaw, clusterVar = rnaCluster, vertices = vt.quat, 
               features = gene.quat, veloGraph = rnaVelo)

## ----plotQuaternaryRGL, webgl=TRUE--------------------------------------------
plotQuaternary(rnaRaw, clusterVar = rnaCluster, vertices = vt.quat, 
               features = gene.quat, veloGraph = rnaVelo, interactive = TRUE)

## ----writeGIF, results="hide"-------------------------------------------------
writeQuaternaryGIF(rnaRaw, clusterVar = rnaCluster, vertices = vt.quat, 
                   features = gene.quat, veloGraph = rnaVelo, 
                   gifPath = "rotating_tetra.gif")

## ----showGIF, echo = FALSE----------------------------------------------------
if (!knitr:::is_latex_output()) {
  knitr::include_graphics("rotating_tetra.gif")
}

## ----unlinkFile, include=FALSE------------------------------------------------
unlink("rotating_tetra.gif")

