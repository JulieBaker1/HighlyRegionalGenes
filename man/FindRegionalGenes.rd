\name{FindRegionalGenes}
\alias{FindRegionalFeatures}
\alias{FindRegionalGenes}
\title{Find Regional Genes}

\description{
  Identifies features that are regionally distributed.
}
\usage{
 obj <- FindRegionalGenes(obj,dims=1:10,nfeatures=2000,overlap_stop=0.75,snn=NULL,do_test,p_threshold)
}

\arguments{
\item{obj}{an object}
\item{dims}{Dimensions of reduction to use as input to build neighbor graph}
\item{nfeatures}{Number of features to select as top regional features}
\item{overlap_stop}{overlap that make the iteration stop}
\item{snn}{the neighborhood relationship that is calculated in advance}

}
\examples{
  pbmc <- ScaleData(pbmc, features = rownames(pbmc))
  pbmc <- RunPCA(pbmc, features = VariableFeatures(obj))
  ElbowPlot(obj)
  pbmc=FindRegionalGenes(pbmc,dims=1:10,nfeatures=2000)
}

