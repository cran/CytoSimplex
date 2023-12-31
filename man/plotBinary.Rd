% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/binary.R
\name{plotBinary}
\alias{plotBinary}
\alias{plotBinary.default}
\alias{plotBinary.Seurat}
\alias{plotBinary.SingleCellExperiment}
\alias{plotBinary.simMat}
\title{Create binary plots}
\usage{
plotBinary(x, ...)

\method{plotBinary}{default}(
  x,
  clusterVar,
  vertices,
  features = NULL,
  byCluster = NULL,
  processed = FALSE,
  method = c("euclidean", "cosine", "pearson", "spearman"),
  force = FALSE,
  sigma = 0.08,
  scale = TRUE,
  dotColor = "grey60",
  returnData = FALSE,
  ...
)

\method{plotBinary}{Seurat}(
  x,
  layer = "counts",
  assay = NULL,
  clusterVar = NULL,
  processed = FALSE,
  ...
)

\method{plotBinary}{SingleCellExperiment}(x, assay.type = "counts", clusterVar = NULL, processed = FALSE, ...)

\method{plotBinary}{simMat}(
  x,
  dotSize = 0.6,
  dotColor = "grey60",
  densLinewidth = 0.8,
  labelColors = c("#3B4992FF", "#EE0000FF"),
  title = NULL,
  ...
)
}
\arguments{
\item{x}{Input data. Can be a \code{matrix} or \code{dgCMatrix} object with
cells as columns, a \code{Seurat} or \code{SingleCellExperiment} object.
"simMat" method takes intermediate values.}

\item{...}{Arguments passed to other methods.}

\item{clusterVar}{A vector/factor assigning the cluster variable to each
column of the matrix object. For "Seurat" method, \code{NULL} (default) for
\code{Idents(x)}, or a variable name in \code{meta.data} slot. For
"SingleCellExperiment" method, \code{NULL} (default) for \code{colLabels(x)},
or a variable name in \code{colData} slot.}

\item{vertices}{Vector of three unique cluster names that will be used for
plotting. Or a named list that groups clusters as three terminal vertices.
There must not be any overlap between groups.}

\item{features}{Valid matrix row subsetting index to select features for
similarity calculation. Default \code{NULL} uses all available features.}

\item{byCluster}{Default \code{NULL} to generate one plot with all cells.
Set \code{"all"} to split cells in plot by cluster and returns a list of
subplots for each cluster as well as the plot including all cells. Otherwise,
a vector of cluster names to generate a list of subplots for the specified
clusters.}

\item{processed}{Logical. Whether the input matrix is already processed.
\code{TRUE} will bypass internal preprocessing and input matrix will be
directly used for similarity calculation. Default \code{FALSE} and raw count
input is recommended. If missing in call, using \code{slot = "counts"} in
"Seurat" method or using \code{assay.type = "counts"} in
"SingleCellExperiment" method will force this argument to be \code{FALSE} and
others for \code{TRUE}.}

\item{method}{Similarity calculation method. Default \code{"euclidean"}.
Choose from \code{"euclidean"}, \code{"cosine"}, \code{"pearson"},
\code{"spearman"}.}

\item{force}{Whether to force calculate the similarity when more then 500
features are detected, which is generally not recommended. Default
\code{FALSE}.}

\item{sigma}{Gaussian kernel parameter that controls the effect of variance.
Only effective when using a distance metric (i.e. \code{method} is
\code{"euclidian"} or \code{"cosine"}). Larger value tighten the dot
spreading on figure. Default \code{0.08}.}

\item{scale}{Whether to min-max scale the distance matrix by clusters.
Default \code{TRUE}.}

\item{returnData}{Logical. Whether to return similarity data instead of
generating plot. Default \code{FALSE}.}

\item{layer}{For "Seurat" method, which layer of the assay to be used.
Default \code{"counts"}.}

\item{assay}{For "Seurat" method, the specific assay to get data from.
Default \code{NULL} to the default assay.}

\item{assay.type}{For "SingleCellExperiment" methods. Which assay to use for
calculating the similarity. Default \code{"counts"}.}

\item{dotSize, dotColor}{Dot aesthetics passed to
\code{\link[ggplot2]{geom_point}}. Default \code{0.6} and \code{"grey60"}.}

\item{densLinewidth}{Density plot line aesthetic. Default \code{0.8}.}

\item{labelColors}{Color of the axis lines and vertex labels. Default
\code{c("#3B4992FF", "#EE0000FF")} (blue and red).}

\item{title}{Title text of the plot. Default \code{NULL}.}
}
\value{
For 'simMat' method, a ggplot object. For other methods, a ggplot
object when \code{splitCluster = FALSE}, or a list of ggplot objects when
\code{splitCluster = TRUE}.
}
\description{
Create binary plots that show similarity between single cells and two
selected terminals in a barycentric coordinate. The two vertices are placed
at the left and right of a 2D plot where x-axis measures the similarity.
Y-axis is jittered for a clear view. A density (histogram) curve is added for
indicating the distribution.

See \code{\link{plotTernary}} manual for more details.
}
\examples{
gene <- selectTopFeatures(rnaRaw, rnaCluster, c("RE", "OS"))
plotBinary(rnaRaw, rnaCluster, c("RE", "OS"), gene)
\donttest{
# Seurat example
library(Seurat)
srt <- CreateSeuratObject(rnaRaw)
Idents(srt) <- rnaCluster
gene <- selectTopFeatures(srt, vertices = c("OS", "RE"))
plotBinary(srt, features = gene, vertices = c("OS", "RE"))
}
\donttest{
# SingleCellExperiment example
library(SingleCellExperiment)
sce <- SingleCellExperiment(assays = list(counts = rnaRaw))
colLabels(sce) <- rnaCluster
gene <- selectTopFeatures(sce, vertices = c("OS", "RE"))
plotBinary(sce, features = gene, vertices = c("OS", "RE"))
}
}
