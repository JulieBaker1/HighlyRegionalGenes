PCA_SD <- function(method){
  gene =rownames(gene_all_methods)[gene_all_methods[,method]<=2000]
  res.pca<- (prcomp(t(scale.data[gene,]), scale = FALSE)$sdev)^2
  top10_explanation=sum(res.pca[1:10]/sum(res.pca))
}






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
  a=detectCores(logical = F)
  if(a<7){
    a=a-1
  }else{
    a=7
  }
  cl <- makeCluster(getOption("cl.cores", a))
  methods_name=colnames(gene_all_methods)
  all.genes=rownames(obj)
  obj=ScaleData(obj,features = all.genes)
  scale.data=GetAssayData(obj,slot="scale.data")

  clusterExport(cl=cl, varlist=c("scale.data","gene_all_methods"), envir=environment())

  top10_explanation <- parLapply(cl,1:length(methods_name),PCA_SD)
  stopCluster(cl)

  methods_rank=rank(-as.data.frame(top10_explanation))
  names(methods_rank)=methods_name
  return(methods_rank)
}

###计算类与类之间的间距

correlation_cal <- function(method){
  cor_cos=matrix(NA,nrow = 2000,ncol = 2000)
  geneset=rownames(gene_all_methods)[gene_all_methods[,method]<=2000]
  cor_cos_num=scale.data[geneset,] %*% t(scale.data[geneset,])
  cor_cos_deno=sqrt(diag(cor_cos_num)) %*% sqrt(t(diag(cor_cos_num)))
  cor_cos=cor_cos_num/cor_cos_deno
  return(cor_cos)
}

#'@export
#'@name gene_correlation
#'@title gene_correlation
#'@param obj an seurat object
#'@param gene_all_methods gene of all methods
gene_correlation <- function(obj,gene_all_methods){
  a=detectCores(logical = F)
  if(a<7){
    a=a-1
  }else{
    a=7
  }
  cl <- makeCluster(getOption("cl.cores", a))
  methods_name=colnames(gene_all_methods)
  all.genes=rownames(obj)

  obj=ScaleData(obj,features = all.genes)
  scale.data=GetAssayData(obj,slot="scale.data")

  clusterExport(cl=cl, varlist=c("scale.data","gene_all_methods"), envir=environment())

  cor_cos <- parLapply(cl,1:length(methods_name),correlation_cal)

  names(cor_cos)=methods_name

  stopCluster(cl)

  number=c()
  for(i in 1:length(methods_name)){
    number[i]=sum(cor_cos[[i]]<0.1)/2
  }
  methods_rank=rank(number)
  names(methods_rank)=methods_name
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
    methods_name=colnames(gene_all_methods)
    gene_mean=my_row_mean_aggregate((GetAssayData(obj,slot = "data")),factor(obj$type))

    mean_cluster_max=apply(gene_mean,1,max)
    mean_cluster_min=apply(gene_mean,1,min)
    mean_cluster_diff=log(mean_cluster_max+1)-log(mean_cluster_min+1)
    median_gene=c()
    for(i in 1:length(methods_name)){
      gene=rownames(gene_all_methods)[gene_all_methods[,i]<=2000]
      median_gene[i]=median(mean_cluster_diff[gene])
    }
    methods_rank=rank(-median_gene)
    names(methods_rank)=methods_name
    return(methods_rank)
}


