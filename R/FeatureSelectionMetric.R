###挑出的基因做完PCA前10个PCA的解释度
#'@export
#'@name PCA10_explanation
#'@title the top10 PCA explanation
#'@param obj an seurat object
#'@param gene_all_methods gene of all methods
#'
#'
#'
#'
PCA10_explanation <- function(obj,gene_all_methods){
    methods_name=colnames(gene_all_methods)[1:7]
    all.genes=rownames(obj)
    obj=ScaleData(obj,features = all.genes)
    explanation=c()
    res.pca=list()
    top10_explanation=c()
    for(i in 1:7){
      gene =rownames(gene_all_methods)[gene_all_methods[,i]<=2000]
      res.pca[[i]] <- (prcomp(t(GetAssayData(obj,slot = "scale.data")[gene,]), scale = FALSE)$sdev)^2
      top10_explanation[i]=sum(res.pca[[i]][1:10]/sum(res.pca[[i]]))
    }
    methods_rank=rank(-top10_explanation)
    names(methods_rank)=methods_name
    return(methods_rank)
}

###计算类与类之间的间距

correlation_cal <- function(method){
  cor_cos=matrix(NA,nrow = 2000,ncol = 2000)
  geneset=rownames(gene_all_methods)[gene_all_methods[,method]<=2000]
  for(i in 1:2000){
    for( j in i:2000){
      cor_cos[i,j]=sum(scale.data[geneset[i],]*scale.data[geneset[j],])/sqrt(sum(scale.data[geneset[i],]^2)*sum(scale.data[geneset[j],]^2))
    }
  }
  return(cor_cos)
}

#'@export
#'@name gene_correlation
#'@title gene_correlation
#'@param obj an seurat object
#'@param gene_all_methods gene of all methods
gene_correlation <- function(obj,gene_all_methods){
  cl <- makeCluster(getOption("cl.cores", 7))
  methods_name=colnames(gene_all_methods)[1:7]
  all.genes=rownames(obj)

  obj=ScaleData(obj,features = all.genes)
  scale.data=GetAssayData(obj,slot="scale.data")

  methods_names=colnames(gene_all_methods)[1:7]
  clusterExport(cl=cl, varlist=c("scale.data","gene_all_methods"), envir=environment())

  cor_cos[[N]] <- parLapply(cl,1:7,correlation_cal)

  names(cor_cos)=methods_name

  stopCluster(cl)

  number=c()
  for(i in 1:7){
    number[i]=sum(cor_cos[[i]][!is.na(cor_cos[[i]])]<0.1)
  }
  methods_rank=rank(number)
  names(methods_rank)=methods_names
  return(methods_rank)
}

##找出重要基因的fc
my_row_mean_aggregate <- function (mat, groups) {
  MAT <- as.matrix(mat)
  x <- split(seq(ncol(MAT)), groups)
  result <- sapply(x, function(a) rowMeans(MAT[, a]))
  return(result)
}
#'@export
#'@name clusterFC
#'@title clusterFC
#'@param obj an seurat object
#'@param gene_all_methods gene of all methods
clusterFC <- function(obj,gene_all_methods){
    methods_name=colnames(gene_all_methods)[1:7]
    gene_mean=my_row_mean_aggregate(exp(GetAssayData(obj,slot = "data")),factor(obj$type))

    mean_cluster_max=apply(gene_mean,1,max)
    mean_cluster_min=apply(gene_mean,1,min)
    mean_cluster_diff=log(mean_cluster_max)-log(mean_cluster_min)
    median_gene=c()
    for(i in 1:7){
      gene=rownames(gene_all_methods)[gene_all_methods[,i]<=2000]
      median_gene[i]=median(mean_cluster_diff[gene])
    }
    methods_rank=rank(-median_gene)
    names(methods_rank)=methods_name
    return(methods_rank)
}


