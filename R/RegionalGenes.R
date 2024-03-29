
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
#'@param obj An Seurat object
#'@param dims Dimensions of reduction to use as input to build neighbor graph
#'@param neigh_num Number of neighbors
#'@param nfeatures Number of features to select as top regional features
#'@param overlap_stop Overlap that make the iteration stop
#'@param snn The neighborhood relationship that is calculated in advance.
#'@param verbose To choose if show the iteration process plot.
#'@param max_iteration The max iteration times
#'@param is.save Whether to save intermediate results
#'@param dir The save directory
#'@details Users should use all the features to run PCA and choose the dimensions used to for this function. Or at least run PCA so that the function can get the SNN.
#'@examples
#'pbmc <- ScaleData(pbmc, features = rownames(pbmc))
#'pbmc <- RunPCA(pbmc, features = VariableFeatures(obj))
#'ElbowPlot(obj)
#'pbmc=FindRegionalGenes(pbmc,dims=1:10,nfeatures=2000)

#'@export
FindRegionalGenes <- function(obj,dims=1:10,nfeatures=2000,overlap_stop=0.75,max_iteration=10,snn=NULL,do_test,p_threshold,verbose=TRUE,neigh_num = 20,is.save=FALSE,dir = ""){
  if(is.save){
    gene_all=list()
  }
  all.genes=rownames(obj)
  block.size=1000
  max.block <- ceiling(x = length(x = all.genes) / block.size)
  cell_num=dim(obj)[2]
  obj=ScaleData(obj,features = all.genes,verbose=FALSE)
  obj_data <- GetAssayData(object = obj,slot="scale.data")  #scale_data

  if(is.null(snn)){

    obj <- FindNeighbors(obj, dims =dims,k.param = neigh_num,verbose=FALSE)

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
    if(is.save){
      gene_all[[1]]=feature_gene
    }


    overlap=0
    count <- 0
    while((overlap<overlap_stop) & (count<max_iteration)){
      count=count+1
      obj=RunPCA(obj,features = feature_gene,verbose=FALSE)
      obj <- FindNeighbors(obj, dims =dims,k.param = neigh_num,verbose=FALSE)
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
      if(is.save){
        gene_all[[count]]=feature_gene
      }
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

  obj[["RNA"]]=AddMetaData(obj[["RNA"]],HRG_score,col.name="HRG.score")
  obj[["RNA"]]=AddMetaData(obj[["RNA"]],HRG_rank,col.name="HRG.rank")
  if(is.save==TRUE){
    names(gene_all)=c(1:length(gene_all))
    saveRDS(as.data.frame(gene_all),paste0(dir,"/gene_iteration.rds"))
  }
    return(obj)

}

#'@import KneeArrower
#'@export
#'@name HRG_elbowplot
#'@title HRG_elbowplot
#'@aliases
#'HRG_elbowplot
#'@usage HRG_elbowplot(obj)
#'@param obj an Seurat object
#'@param method the method to define the knee point. Value can be "first" for first derivative cutoff or "curvature" for maximum curvature cutoff.
#'@description Find regional gene number
#'

HRG_elbowplot <- function(obj,method = "curvature"){
  score = obj[["RNA"]][["HRG.score"]]
  score = as.matrix(score)[,1]
  # names(score) = rownames(pbmc[["RNA"]][["HRG.score"]])
  score = sort(score,decreasing = TRUE)
  elbow_point  = KneeArrower::findCutoff(1:length(score),score,method)
  gene_num = floor(elbow_point$x)
  plot(1:length(score),score)
  points(gene_num,score[gene_num],col="red", pch = 19,lwd = 4)
  return(gene_num)
}





#'@export
#'@name RegionalGenes
#'@title RegionalGenes
#'@aliases
#'RegionalGenes
#'RegionalFeatures
#'@usage RegionalGenes(obj)
#'@param obj an seurat object
#'@param nfeatures Number of features to select as output regional features
#'@description Get and set regional feature information.
#'

RegionalGenes <- function(obj,nfeatures = 2000){
  gene_metadata=obj@assays$RNA@meta.features
  HRG_rank = gene_metadata["HRG.rank"]
  HRG_rank = as.matrix(HRG_rank)[,1]
  return(names(sort(HRG_rank))[1:nfeatures])
}


