
#'@import Seurat
#'@name
#'FindRegionalGenes
#'@aliases
#'FindRegionalFeatures
#'FindRegionalGenes
#'@title
#'Find Regional Genes
#'@description
#'Identifies features that are regionally distributed.
#'@usage
#'obj <- FindRegionalGenes(obj,dims=1:10,nfeatures=2000,overlap_stop=0.75,max_iteration=10,snn=NULL,do_test,p_threshold,is.save=FALSE,dir=getwd())
#'@param obj an Seurat object
#'@param dims Dimensions of reduction to use as input to build neighbor graph
#'@param nfeatures Number of features to select as top regional features
#'@param overlap_stop overlap that make the iteration stop
#'@param snn the neighborhood relationship that is calculated in advance.
#'@param verbose to choose if show the iteration process plot.
#'@param max_iteration the max iteration times
#'@param is.save whether to save intermediate results
#'@param dir the save directory
#'@details the user should use all the features to run PCA and choose the dimensions used to for this function. Or at least run PCA so that the function can get the SNN.
#'@examples
#'pbmc <- ScaleData(pbmc, features = rownames(pbmc))
#'pbmc <- RunPCA(pbmc, features = VariableFeatures(obj))
#'ElbowPlot(obj)
#'pbmc=FindRegionalGenes(pbmc,dims=1:10,nfeatures=2000)

#'@export
FindRegionalGenes <- function(obj,dims=1:10,nfeatures=2000,overlap_stop=0.75,max_iteration=10,snn=NULL,do_test,p_threshold,verbose=TRUE,is.save=FALSE,dir = ""){
  if(is.save==TRUE){
    gene_all=list()
  }
  all.genes=rownames(obj)
  block.size=1000
  max.block <- ceiling(x = length(x = all.genes) / block.size)
  cell_num=dim(obj)[2]
  obj=ScaleData(obj,features = all.genes,verbose=FALSE)
  obj_data <- GetAssayData(object = obj,slot="scale.data")  #scale_data

  if(is.null(snn)){

    obj <- FindNeighbors(obj, dims =dims,verbose=FALSE)

    snn <- obj$RNA_snn

    diag(snn) <- 0


    if (verbose) {
      message("calculating gene score")
      pb <- txtProgressBar(min = 0, max = max.block, style = 3, file = stderr())
    }
    HRG_score=c()
    for(i in 1:max.block){
      my.inds <- ((block.size * (i - 1)):(block.size * i - 1)) + 1
      my.inds <- my.inds[my.inds <= length(x = all.genes)]
      for(index in my.inds){
        data_temp=as.matrix(obj_data[index,])
        HRG_score[index]=as.numeric(t(data_temp)%*%snn%*%(data_temp))
      }
      if(verbose){
        setTxtProgressBar(pb = pb, value = i)
      }
    }
    if (verbose) {
      close(con = pb)
    }

    names(HRG_score) <- all.genes
    feature_gene=names(sort(HRG_score,decreasing = TRUE))[1:nfeatures]
    gene_all[[1]]=feature_gene

    overlap=0
    count <- 0
    while((overlap<overlap_stop) & (count<max_iteration)){
      count=count+1
      obj=RunPCA(obj,features = feature_gene,verbose=FALSE)
      obj <- FindNeighbors(obj, dims =dims,verbose=FALSE)
      snn <- obj$RNA_snn
      diag(snn) <- 0

      if (verbose) {
        message("calculating gene score")
        pb <- txtProgressBar(min = 0, max = max.block, style = 3, file = stderr())
      }
      for(i in 1:max.block){
        my.inds <- ((block.size * (i - 1)):(block.size * i - 1)) + 1
        my.inds <- my.inds[my.inds <= length(x = all.genes)]
        for(index in my.inds){
          data_temp=as.matrix(obj_data[index,])
          HRG_score[index]=as.numeric(t(data_temp)%*%snn%*%(data_temp))
        }
        if(verbose){
          setTxtProgressBar(pb = pb, value = i)
        }
      }
      if (verbose) {
        close(con = pb)
      }
      feature_gene_new=names(sort(HRG_score,decreasing = TRUE))[1:nfeatures]
      overlap=length(intersect(feature_gene_new,feature_gene))/nfeatures

      if(verbose){
        message(paste0("overlap is ",overlap))
      }

      feature_gene=feature_gene_new
      gene_all[[count]]=feature_gene
    }
  } else{
    size=dim(snn)
    if(size[1]!=cell_num | size[2]!=cell_num){
      stop("the snn should be cellnumberxcellnumber matrix")
    }
    if (verbose) {
      message("calculating gene score")
      pb <- txtProgressBar(min = 0, max = max.block, style = 3, file = stderr())
    }
    HRG_score=c()

    for(i in 1:max.block){
      my.inds <- ((block.size * (i - 1)):(block.size * i - 1)) + 1
      my.inds <- my.inds[my.inds <= length(x = all.genes)]
      for(index in my.inds){
        data_temp=as.matrix(obj_data[index,])
        HRG_score[index]=as.numeric(t(data_temp)%*%snn%*%(data_temp))
      }
      if(verbose){
        setTxtProgressBar(pb = pb, value = i)
      }
    }
    if (verbose) {
      close(con = pb)
    }
    names(HRG_score) <- all.genes
  }

  HRG_rank=rank(-HRG_score)

  HRG_regional=(HRG_rank<=nfeatures)

  obj[["RNA"]]=AddMetaData(obj[["RNA"]],HRG_score,col.name="HRG.score")
  obj[["RNA"]]=AddMetaData(obj[["RNA"]],HRG_rank,col.name="HRG.rank")
  obj[["RNA"]]=AddMetaData(obj[["RNA"]],HRG_regional,col.name="HRG.regional")
  print(getwd())
  if(is.save==TRUE){
    names(gene_all)=c(1:length(gene_all))
    saveRDS(as.data.frame(gene_all),paste0(dir,"/gene_iteration.rds"))
  }
    return(obj)

}



#'@export
#'@name RegionalGenes
#'@title RegionalGenes
#'@aliases
#'RegionalGenes
#'RegionalFeatures
#'@usage RegionalGenes(obj)
#'@param obj an seurat object
#'@description Get and set regional feature information.
#'
#'

RegionalGenes <- function(obj){
  gene_metadata=obj@assays$RNA@meta.features
  features_num=sum(gene_metadata$HRG.regional)
  rank=gene_metadata$HRG.rank
  names(rank)=rownames(gene_metadata)
  return(names(sort(rank))[1:features_num])
}


