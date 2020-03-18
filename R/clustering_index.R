


#' @inheritParams ARI
#' @param object an Seurat object that contains cluster result and the inherent category
#'
#' @param name the colname that contains the genes' inherent category of the cell metadata
#'
#' @return the value of an external index that evaluate the clustering result,range from -1 to 1, the larger, the better.
#'



ARI <- function(object,name){

  metadata=object@meta.data

  cell_type=metadata[[name]]
  cell_type=factor(cell_type,order=TRUE)
  class_number=length(levels(cell_type))

  seurat_clusters=metadata$seurat_clusters
  cluster_number=length(levels(seurat_clusters))


  final_result=data.frame(cell_type,seurat_clusters)

  cell_number=dim(metadata)[1]



  #####ARI######

  result_temp=final_result

  pair=array(0,dim=c(class_number,cluster_number))
  for(i in c(1:cell_number)){
    pair[result_temp[i,1],result_temp[i,2]]= pair[result_temp[i,1],result_temp[i,2]]+1
  }

  index=sum(choose(pair,2))
  a=rep(0,class_number)
  total_a=0
  b=rep(0,cluster_number)
  total_b=0
  for(i in c(1:class_number)){
    a[i]=sum(pair[i,])
    total_a=total_a+choose(a[i],2)}
  for(i in c(1:cluster_number)){
    b[i]=sum(pair[,i])
    total_b=total_b+choose(b[i],2)
  }
  expected_index=total_a*total_b/choose(sum(a),2)
  max_index=(total_a+total_b)/2

  ARI=(index-expected_index)/(max_index-expected_index)

  print(ARI)
}


#' @inheritParams NMI
#' @param object an Seurat object that contains cluster result and the inherent category
#'
#' @param name the colname that contains the genes' inherent category of the cell metadata
#'
#' @return the value of an external index that evaluate the clustering result,range from -1 to 1, the larger, the better.
#'


NMI <- function(object,name){

  metadata=object@meta.data

  cell_type=metadata[[name]]
  cell_type=factor(cell_type,order=TRUE) #?൱?ڸ?cell_typeһ????
  class_number=length(levels(cell_type))

  seurat_clusters=metadata$seurat_clusters
  cluster_number=length(levels(seurat_clusters))


  final_result=data.frame(cell_type,seurat_clusters)

  cell_number=dim(metadata)[1]





  result_temp=final_result

  pair=array(0,dim=c(class_number,cluster_number))
  for(i in c(1:cell_number)){
    pair[result_temp[i,1],result_temp[i,2]]= pair[result_temp[i,1],result_temp[i,2]]+1
  }

  pair=pair/sum(pair[,])
  a=rep(0,class_number)
  HX=0
  b=rep(0,cluster_number)
  HY=0
  for(i in c(1:class_number)){
    a[i]=sum(pair[i,])
    HX=HX+a[i]*log2(a[i])
  }
  HX=-HX
  for(i in c(1:cluster_number)){
    b[i]=sum(pair[,i])
    if(b[i]){
      HY=HY+b[i]*log2(b[i])
    }
  }
  HY=-HY
  MI=0
  for(i in c(1:class_number)){
    for(j in c(1:cluster_number)){
      if(pair[i,j]&a[i]&b[i]){
        MI=MI+pair[i,j]*log2(pair[i,j]/(a[i]*b[j]))
      }
    }
  }

  NMI=2*MI/(HX+HY)

  print(NMI)


}


#' @inheritParams DBI
#' @param object an Seurat object that contains cluster result
#'
#' @param PCnumber the number of PCs to use to calculate the distance of two cells.It is equal to the number of PCs to clustering.
#'
#' @return the value of internal index that evaluate the clustering result. The smaller, the better.


DBI <- function(object,PCnumber){


  PCA=object@reductions[["pca"]]@cell.embeddings[,1:PCnumber]

  seurat_clusters=object@meta.data$seurat_clusters

  eps=1e-8
  order_length=1:dim(PCA)[1]


  cluster_number=length(levels(seurat_clusters))#找到一共有多少个类
  distance=rep(0,cluster_number)#类内距离
  total_PCA=dim(PCA)[2]#一共有多少个PCA
  center=array(0,dim=c(cluster_number,total_PCA))
  for(cluster in c(0:(cluster_number-1))){
    order_index=order_length[seurat_clusters==cluster]
    total_number=length(order_index)
    PCA_cluster=PCA[order_index,]


    for(PCA_number in c(1:total_PCA)){
      center[cluster+1,PCA_number]=sum(PCA_cluster[,PCA_number])/(total_number+eps)
    }#计算中心点的坐标

    distance[cluster+1]=sum(dist(PCA_cluster,p=2))/choose(total_number,2)#计算类内距离

  }#找到属于第一类的索引
  distance_cluster_cluster=dist(center,p=2)
  for(cluster1 in c(1:cluster_number)){
    difference=0
    max_difference=rep(0,cluster_number)
    for(cluster2 in c(1:cluster_number)){
      if(cluster1!=cluster2){
        difference_temp=(distance[cluster1]+distance[cluster2])/(dist(center[c(cluster1,cluster2),],p=2))
        if(difference_temp>difference) difference=difference_temp
      }
    }
    max_difference[cluster1]=difference
  }
  DBI=sum(max_difference)/cluster_number

  print(DBI)
  #DBI越小越好
}


#' @inheritParams DVI
#' @param object an Seurat object that contains cluster result
#'
#' @param PCnumber the number of PCs to use to calculate the distance of two cells.It is equal to the number of PCs to clustering.
#'
#' @return the value of internal index that evaluate the clustering result. The bigger, the better.

DVI <- function(object,PCnumber){
  #此处存在一个小问题，到底该取多少个PCA计算距离

  PCA=object@reductions[["pca"]]@cell.embeddings[,1:PCnumber]

  seurat_clusters=object@meta.data$seurat_clusters

  eps=1e-8
  order_length=1:dim(PCA)[1]



  cluster_number=length(levels(seurat_clusters))#找到一共有多少个类
  PCA_cluster=list()
  distance=rep(0,cluster_number)#类内距离
  total_PCA=dim(PCA)[2]#一共有多少个PCA
  center=array(0,dim=c(cluster_number,total_PCA))#中心点的坐标
  total_number=rep(0,cluster_number)
  for(cluster in 0:(cluster_number-1)){
    order_index=order_length[seurat_clusters==cluster]
    total_number[cluster+1]=length(order_index)
    PCA_cluster[[cluster+1]]=PCA[order_index,]

    distance[cluster+1]=max(dist(PCA_cluster[[cluster+1]],p=2))#计算类内距离

  }#找到属于第一类的索引
  max_distance=max(distance)
  index=rep(0,choose(cluster_number,2))
  count=1
  min_distance=1000000
  for(cluster1 in 1:(cluster_number-1)){
    for(cluster2 in (cluster1+1):cluster_number){
      for(i in 1:dim(PCA_cluster[[cluster1]])[1]){
        for(j in 1:dim(PCA_cluster[[cluster2]])[1]){
          min_distance_temp=sqrt(sum((PCA_cluster[[cluster1]][i,]-PCA_cluster[[cluster2]][j,])^2))
          if(min_distance_temp<min_distance) min_distance=min_distance_temp
        }
      }
    }
  }
  DUNN=min_distance/max_distance

  print(DUNN)
  #DUNN越大越好
}




