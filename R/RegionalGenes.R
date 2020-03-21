
#'@param obj an Seurat object
#'@param dims the dimentions of PCA  used to construct SNN
#'@param nfeatures the number of features for downstream analysis
#'@param overlap_stop when the overlap come to this number, iterate will stop.
#'@param snn the neighborhood relationship that is calculated in advance.
#'@details the user should use all the features to run PCA and choose the dimensions used to for this function. Or at least run PCA so that the function can get the SNN.
#'
#'



FindRegionalGenes <- function(obj,dims=1:10,nfeatures=2000,overlap_stop=0.75,snn=NULL,do_test,p_threshold,verbose=TRUE){


  all.genes=rownames(obj)
  cell_num=dim(obj)[2]
  obj=ScaleData(obj,features = all.genes,verbose=FALSE)
  obj_data <- GetAssayData(object = obj,slot="scale.data")  #scale_data

  if(is.null(snn)){

    obj <- FindNeighbors(obj, dims =dims,verbose=FALSE)

    snn <- obj$RNA_snn

    diag(snn) <- 0

    HRG_score=c()
    for(index in 1:length(all.genes)){
      data_temp=as.matrix(obj_data[index,])
      HRG_score[index]=as.numeric(t(data_temp)%*%snn%*%(data_temp))
    }
    names(HRG_score) <- all.genes
    feature_gene=names(sort(HRG_score,decreasing = TRUE))[1:nfeatures]


    overlap=0
    gp <- ggplot()+geom_hline(yintercept=overlap_stop,  color = "red")+ylim(0,1)+scale_x_continuous(breaks=c(1:20))
    count <- 1
    while(overlap<overlap_stop){

      obj=RunPCA(obj,features = feature_gene,verbose=FALSE)
      obj <- FindNeighbors(obj, dims =dims,verbose=FALSE)
      snn <- obj$RNA_snn
      diag(snn) <- 0

      for(index in 1:length(all.genes)){
        data_temp=as.matrix(obj_data[index,])
        HRG_score[index]=as.numeric(t(data_temp)%*%snn%*%(data_temp))
      }
      feature_gene_new=names(sort(HRG_score,decreasing = TRUE))[1:nfeatures]
      overlap=length(intersect(feature_gene_new,feature_gene))/nfeatures

      if(verbose){
        print(overlap)
        gp <- gp + geom_point(data = data.frame(time=count,overlap=overlap) , aes(x=time , y=overlap))
        count=count+1
        print(gp)
      }

      feature_gene=feature_gene_new
    }
  } else{
    size=dim(snn)
    if(size[1]!=cell_num | size[2]!=cell_num){
      stop("the snn should be cellnumberxcellnumber matrix")
    }

    HRG_score=c()
    for(index in 1:length(all.genes)){
      data_temp=as.matrix(obj_data[index,])
      HRG_score[index]=as.numeric(t(data_temp)%*%snn%*%(data_temp))
    }
    names(HRG_score) <- all.genes
  }




  HRG_rank=rank(-HRG_score)

  HRG_regional=(HRG_rank<=nfeatures)

  obj[["RNA"]]=AddMetaData(obj[["RNA"]],HRG_score,col.name="HRG.score")
  obj[["RNA"]]=AddMetaData(obj[["RNA"]],HRG_rank,col.name="HRG.rank")
  obj[["RNA"]]=AddMetaData(obj[["RNA"]],HRG_regional,col.name="HRG.regional")
  return(obj)

}

RegionalGenes <- function(obj){
  gene_metadata=obj@assays$RNA@meta.features
  features_num=sum(gene_metadata$HRG.regional)
  rank=gene_metadata$HRG.rank
  names(rank)=rownames(gene_metadata)
  return(names(sort(rank))[1:features_num])
}


