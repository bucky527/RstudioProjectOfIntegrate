library(dplyr, quietly = T)
library(getopt, quietly = T)
library(Matrix, quietly = T)
library(future,quietly = T)
library(Seurat)
library(parallel)


ncores<-detectCores()-4

seuratdata<-CreateSeuratObject(counts = data,meta.data = metadata)
seuratdata<- FindVariableFeatures(seuratdata, selection.method = "vst",
                                  nfeatures = 2000, verbose = FALSE)
seuratdata.list<-SplitObject(seuratdata,split.by = "batch")

for (i in 1:length(seuratdata.list)) {
  seuratdata.list[[i]] <- NormalizeData(seuratdata.list[[i]], verbose = FALSE)
  #seuratdata.list[[i]] <- FindVariableFeatures(seuratdata.list[[i]], selection.method = "vst",
  #nfeatures = 2000, verbose = FALSE)
}
seurat.anchors<-FindIntegrationAnchors(seuratdata.list)
numb <- seurat.anchors@offsets
anchors <- seurat.anchors@anchors
anchors<-anchors[1:(nrow(anchors)/2),]
anchors<-anchors[order(anchors$score,decreasing = TRUE),]
#不加偏移量的anchors提取

#去除重复的anchors 对于重复的留取评分最高的anchors
#anchors<-anchors[!duplicated(anchors$cell1),]

dataset1<-seuratdata.list[[1]]@assays$RNA@counts[seuratdata.list[[1]]@assays$RNA@var.features,]
dataset2<-seuratdata.list[[2]]@assays$RNA@counts[seuratdata.list[[2]]@assays$RNA@var.features,]
Hdataset1<-seuratdata.list[[1]]@assays$RNA@counts
Hdataset2<-seuratdata.list[[2]]@assays$RNA@counts

dataset1<-as.matrix(dataset1)
dataset2<-as.matrix(dataset2)
Hdataset1<-as.matrix(Hdataset1)
Hdataset2<-as.matrix(Hdataset2)
anchorset1<-dataset1[,anchors$cell1]
anchorset2<-dataset2[,anchors$cell2]

sigma<-10
sigma<-2*(sigma^2)



# 初始化
cl <- makeCluster(ncores)
clusterExport(cl,varlist = list("anchorset1","dataset1","sigma"))

GuassianWeight<-parApply(cl=cl,X=dataset1, MARGIN=2,FUN =function(c,anch,edist=0){
  edist<-(anch-c)
  w<-matrix(1:ncol(edist),nrow=1)
  for(i in 1:ncol(edist)){w[i]<-(t(edist[,i])%*%edist[,i])}
  return(exp(-t(w)/sigma))
},anchorset1)


AnchorWeight<-parApply(cl,GuassianWeight, MARGIN = 2, FUN = function(c){
  c<-c/sum(c)
  return(c)
})




Aerror<-(anchorset1-anchorset2)

HAerror<-Hdataset1[,anchors$cell1]-Hdataset2[,anchors$cell2]
clusterExport(cl,varlist = list("AnchorWeight"))

correctError<-parApply(cl,HAerror,MARGIN = 1,FUN = function(a){
  return(a%*%AnchorWeight)
})


#归还线程
stopCluster(cl)

correctError<-t(correctError)
# correctError<-HAerror%*%AnchorWeight




integratedDataset1<-Hdataset1-correctError

integratedData<-cbind(integratedDataset1,Hdataset2)
integratedData[is.na(integratedData)]<-0

seuratInte<-IntegrateData(seurat.anchors,dims = 1:30)
seuratInte<-FindVariableFeatures(seuratInte,selection.method = "vst",nfeatures = 2000)
seuratInte<-ScaleData(seuratInte)%>%RunPCA()%>%RunUMAP(reduction="pca",dims=1:30)
seuratp<-DimPlot(seuratInte,reduction = "umap",group.by = "batch")


meta1<-data.frame(batch=seuratdata.list[[1]]@meta.data$batch)
meta2<-data.frame(batch=seuratdata.list[[2]]@meta.data$batch)
rownames(meta1)<-rownames(seuratdata.list[[1]]@meta.data)
rownames(meta2)<-rownames(seuratdata.list[[2]]@meta.data)

mergemeta<-rbind(meta1,meta2)

mydata<-CreateSeuratObject(integratedData,meta.data = mergemeta)
mydata<-FindVariableFeatures(mydata,selection.method = "vst",nfeatures = 2000)
mydata<-ScaleData(mydata)
mydata<-RunPCA(mydata)
mydata<-RunUMAP(mydata,reduction = "pca",dims =1:30 )

myp<-DimPlot(mydata,reduction = "umap",group.by = "batch")


#seuratdata<-ScaleData(seuratdata)%>%RunPCA()%>%RunUMAP(reduction="pca",dims=1:30)
#DimPlot(seuratdata,reduction = "umap",group.by = "batch")+DimPlot(seuratInte,reduction = "umap",group.by = "batch")
