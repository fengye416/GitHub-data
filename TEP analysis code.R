#加载R包
rm(list=ls())
options(stringsAsFactors = F) 
library(Seurat)
library(data.table)
library(dplyr)
#对原始文件进行处理，将原始文件分别整理为barcodes.tsv.gz，features.tsv.gz和matrix.mtx.gz并放到到各自的文件夹中。
#首先，将下载文件解压
setwd("C:\\Users\\Xu416\\Desktop\\data6")  
dir='GSE214607_RAW/' 
fs=list.files('GSE214607_RAW/','^GSM')
fs
library(tidyverse)
samples=str_split(fs,'_',simplify = T)[,1] #按“—”分开取第一列数据

##处理数据，将原始文件分别整理为barcodes.tsv.gz，features.tsv.gz和matrix.mtx.gz到各自的文件夹
#批量将文件名改为 Read10X()函数能够识别的名字
lapply(unique(samples),function(x){
  # x = unique(samples)[1]
  y=fs[grepl(x,fs)]
  folder=paste0("GSE214607_RAW/", paste(str_split(y[1],'_',simplify = T)[,1:2], collapse = "_")) #取出y[1]中的字符串，按照下划线进行分割，取出前两列，然后将这两列的值拼接成一个新字符串，用作创建子文件夹的名称
  dir.create(folder,recursive = T)
  #为每个样本创建子文件夹
  file.rename(paste0("GSE214607_RAW/",y[1]),file.path(folder,"barcodes.tsv.gz"))
  #重命名文件，并移动到相应的子文件夹里
  file.rename(paste0("GSE214607_RAW/",y[2]),file.path(folder,"features.tsv.gz"))
  file.rename(paste0("GSE214607_RAW/",y[3]),file.path(folder,"matrix.mtx.gz"))
})

#######3.数据读取与合并####
dir='GSE214607_RAW/'
samples=list.files( dir )
samples 
sceList = lapply(samples,function(pro){ 
  # pro=samples[1] 
  print(pro)  
  tmp = Read10X(file.path(dir,pro )) 
  if(length(tmp)==2){
    ct = tmp[[1]] 
  }else{ct = tmp}
  sce =CreateSeuratObject(counts =  ct ,
                          project =  pro  ,
                          min.cells = 5,
                          min.features = 300 )
  return(sce)#返回创建的Seurat对象，将其存储在sceList中。
}) 
#sceList实际上包含了10个样本的Seurat对象
View(sceList)
#查看其中一个
PT1 <- sceList[1]
#使用Seurat包的merge函数，将十个Seurat合并成一个Seurat对象
do.call(rbind,lapply(sceList, dim))
sce.all=merge(x=sceList[[1]],
              y=sceList[ -1 ],
              add.cell.ids = samples  ) 
names(sce.all@assays$RNA@layers)
#使用JoinLayers函数对layers进行合并
sce.all[["RNA"]]$counts 
# Alternate accessor function with the same result
LayerData(sce.all, assay = "RNA", layer = "counts")
#看看合并前后的sce变化
sce.all
sce.all <- JoinLayers(sce.all)
sce.all
#查看”sce.all“内部的一些信息，以此来检查数据是否完整
dim(sce.all[["RNA"]]$counts ) #dim() 函数用于获取矩阵的维度，即行数和列数
as.data.frame(sce.all@assays$RNA$counts[1:10, 1:2])
head(sce.all@meta.data, 10)
table(sce.all$orig.ident) #统计sce.all对象中每个细胞的原始标识符出现的次数，并返回一个表格。这有助于了解数据集中不同样本或生物学条件的分布情况
length(sce.all$orig.ident) #获取样本或个体的数量。

####4.添加meta.data分组信息###
#查看现有的meta.data信息
library(stringr)
phe = sce.all@meta.data
table(phe$orig.ident)
View(phe)
#函数 str_split 用于拆分字符串
phe$patient = str_split(phe$orig.ident,'[_]',simplify = T)[,2]
table(phe$patient)
#添加转移部位分组信息
phe$patient = phe$orig.ident
phe$patient = gsub("GSM\\d+_NC", "NC", phe$patient) 

#添加患者来源信息
phe$sample = phe$orig.ident
table(phe$sample)
phe$sample = gsub("GSM6613036_NC", "sample1", phe$sample) 
phe$sample = gsub("GSM6613037_NC", "sample2", phe$sample) 
phe$sample = gsub("GSM6613038_NC", "sample3", phe$sample) 
phe$sample = gsub("GSM6613039_NC", "sample4", phe$sample) 
phe$sample = gsub("GSM6613040_NC", "sample5", phe$sample) 
table(phe$sample)

sce.all@meta.data = phe
View(phe)
saveRDS(sce.all, file = "GSE189357.rds")

###(二)Seurat V5标准流程###
rm(list=ls())
options(stringsAsFactors = F) 
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(data.table)
library(dplyr)
#install.packages("scRNAstat")
library(scRNAstat)
setwd('C:\\Users\\Xu416\\Desktop\\data6')
sce.all <- readRDS("GSE189357.rds")
#sce.all <- readRDS("C:\\Users\\Xu416\\Desktop\\data6\\GSE189357.rds")
#创建了新的文件夹”1-QC“
dir.create("./1-QC")
setwd("./1-QC")
# 如果过滤的太狠，就需要去修改这个过滤代码
sce.all.filt = basic_qc(sce.all)
print(dim(sce.all))
print(dim(sce.all.filt))
##细胞减少了一点
#setwd('../')
#getwd()
fivenum(sce.all.filt$percent_ribo)
table(sce.all.filt$nFeature_RNA> 5)



#计算线粒体基因比例
mito_genes=rownames(sce.all)[grep("^MT-", rownames(sce.all),ignore.case = T)] 
print(mito_genes) #可能是13个线粒体基因，小鼠数据基因名为小写"^mt-"
#sce.all=PercentageFeatureSet(sce.all, "^MT-", col.name = "percent_mito")
#使用PercentageFeatureSet函数计算线粒体基因的百分比
sce.all=PercentageFeatureSet(sce.all, features = mito_genes, col.name = "percent_mito")
fivenum(sce.all@meta.data$percent_mito)

#计算核糖体基因比例
ribo_genes=rownames(sce.all)[grep("^Rp[sl]", rownames(sce.all),ignore.case = T)]
print(ribo_genes)
sce.all=PercentageFeatureSet(sce.all,  features = ribo_genes, col.name = "percent_ribo")
fivenum(sce.all@meta.data$percent_ribo)

#计算红血细胞基因比例
Hb_genes=rownames(sce.all)[grep("^Hb[^(p)]", rownames(sce.all),ignore.case = T)]
print(Hb_genes)
sce.all=PercentageFeatureSet(sce.all,  features = Hb_genes,col.name = "percent_hb")
fivenum(sce.all@meta.data$percent_hb)
head(sce.all@meta.data)

#可视化细胞的上述比例情况
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito",
           "percent_ribo", "percent_hb")
feats <- c("nFeature_RNA", "nCount_RNA")
p1=VlnPlot(sce.all, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 2) + 
  NoLegend()
p1 
w=length(unique(sce.all$orig.ident))/3+5;w
ggsave(filename="Vlnplot1.pdf",plot=p1,width = w,height = 5)
feats <- c("percent_mito", "percent_ribo", "percent_hb")
p2=VlnPlot(sce.all, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3, same.y.lims=T) + 
  scale_y_continuous(breaks=seq(0, 100, 5)) +
  NoLegend()
p2 
w=length(unique(sce.all$orig.ident))/2+5;w
ggsave(filename="Vlnplot2.pdf",plot=p2,width = w,height = 5)

p3=FeatureScatter(sce.all, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)
p3
ggsave(filename="Scatterplot.pdf",plot=p3)
#过滤指标1:最少表达基因数的细胞 and 最少表达细胞数的基因
#一般来说，在CreateSeuratObject的时候已经是进行了这个过滤操作，如果后期看到了自己的单细胞降维聚类分群结果很诡异，就可以回过头来看质量控制环节，先走默认流程即可。
#过滤指标1:最少表达基因数的细胞 and 最少表达细胞数的基因
#使用WhichCells函数从sce.all对象中选择表达超过500个基因的细胞
#选择在至少4个细胞中表达的基因
if(F){
  selected_c <- WhichCells(sce.all, expression = nFeature_RNA > 500)
  selected_f <- rownames(sce.all)[Matrix::rowSums(sce.all@assays$RNA$counts > 0 ) > 3]
  sce.all.filt <- subset(sce.all, features = selected_f, cells = selected_c) #创建一个新的SingleCellExperiment对象sce.all.filt，它只包含之前步骤中筛选出的基因和细胞
  dim(sce.all) #打印原始数据集
  dim(sce.all.filt) #筛选后数据集的维度
}

sce.all.filt =  sce.all #将原始数据集赋值给sce.all.filt
# par(mar = c(4, 8, 2, 1))
# 这里的C 这个矩阵，有一点大，可以考虑随抽样 
#从sce.all.filt中随机抽取100个细胞的RNA计数数据
C=subset(sce.all.filt,downsample=100)@assays$RNA$counts
dim(C)
C=Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100 #对矩阵C进行归一化处理，使得每个基因在所有细胞中的总表达量为100

most_expressed <- order(apply(C, 1, median), decreasing = T)[50:1] #找出表达量中位数最高的前50个基因

pdf("TOP50_most_expressed_gene.pdf",width=14)
#绘制一个水平的箱线图，展示前50个最表达基因在不同细胞中的表达水平
boxplot(as.matrix(Matrix::t(C[most_expressed, ])),
        cex = 0.1, las = 1, 
        xlab = "% total count per cell", 
        col = (scales::hue_pal())(50)[50:1], 
        horizontal = TRUE)
dev.off()
rm(C)

#过滤指标2:线粒体/核糖体基因比例(percent_mito < 25，percent_ribo > 3，percent_hb < 1)
selected_mito <- WhichCells(sce.all.filt, expression = percent_mito < 25)
selected_ribo <- WhichCells(sce.all.filt, expression = percent_ribo > 3)
selected_hb <- WhichCells(sce.all.filt, expression = percent_hb < 1 )
length(selected_hb)
length(selected_ribo)
length(selected_mito)
sce.all.filt <- subset(sce.all.filt, cells = selected_mito)
sce.all.filt <- subset(sce.all.filt, cells = selected_ribo)
sce.all.filt <- subset(sce.all.filt, cells = selected_hb)
dim(sce.all.filt) #获取矩阵的维度
table(sce.all.filt$orig.ident)
length(sce.all.filt$orig.ident)

#可视化过滤后的情况
feats <- c("nFeature_RNA", "nCount_RNA")
p1_filtered=VlnPlot(sce.all.filt, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 2) + 
  NoLegend()
w=length(unique(sce.all.filt$orig.ident))/3+5;w 
ggsave(filename="Vlnplot1_filtered.pdf",plot=p1_filtered,width = w,height = 5)

feats <- c("percent_mito", "percent_ribo", "percent_hb")
p2_filtered=VlnPlot(sce.all.filt, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3) + 
  NoLegend()
w=length(unique(sce.all.filt$orig.ident))/2+5;w 
ggsave(filename="Vlnplot2_filtered.pdf",plot=p2_filtered,width = w,height = 5) 


##2、Harmony整合多个单细胞样品
#单细胞测序常用的去批次算法有Harmony，CCA，RPCA,FastMNN,scVI等
#数据标准化--数据标准化的方法是通过对原始表达值进行对数转换"LogNormalize"，使其总体更加符合正态分布。
set.seed(10086)
table(sce.all.filt$orig.ident)
setwd('../')
getwd()
#创建了新的文件夹”2-harmony“
dir.create("./2-harmony")
setwd("./2-harmony")
sce.all.filt <- NormalizeData(sce.all.filt, 
                              normalization.method = "LogNormalize",
                              scale.factor = 1e4) 
pbmc.obj=sce.all.filt
# 默认情况下，使用每个数据集的 2，000 个基因
pbmc.obj <- FindVariableFeatures(pbmc.obj,
                                 selection.method = "vst",
                                 nfeatures = 2000)

# 获取前10个高变基因
top10.hgvs <- head(VariableFeatures(pbmc.obj), 10)
# [1] "PPBP"   "LYZ"    "S100A9" "IGLL5"  "GNLY"   "FTL"    "PF4"    "FTH1"   "GNG11" 
# [10] "S100A8"
top10.hgvs
p3 <- VariableFeaturePlot(pbmc.obj)
p3
p4 <- LabelPoints(plot = p3, points = top10.hgvs, repel = TRUE)
p4

#筛选高变基因-在某些细胞中高度表达，而在其他细胞中低度表达的基因。
sce.all.filt <- FindVariableFeatures(sce.all.filt)
p5 <- VariableFeaturePlot(sce.all.filt) 
p5
#数据归一化-将每个基因在所有细胞中的均值变为0，方差标为1，赋予每个基因在下游分析中同样的权重，不至于使高表达基因占据主导地位
sce.all.filt <- ScaleData(sce.all.filt)

#PCA线性降维
sce.all.filt <- RunPCA(sce.all.filt, features = VariableFeatures(object = sce.all.filt))

##可视化PCA结果
#图中通常会显示每个基因在不同主成分上的载荷值。载荷值可以是正的或负的，分别表示基因表达量与主成分的正相关或负相关。载荷值的大小表示了基因对主成分的贡献程度，载荷值越大，说明该基因在该主成分上的贡献越大
VizDimLoadings(sce.all.filt, dims = 1:2, reduction = "pca") #用于可视化每个主成分的基因载荷（即基因在主成分上的权重）
DimPlot(sce.all.filt, reduction = "pca") + NoLegend() #绘制降维后的二维散点图，展示不同细胞在主成分上的分布
DimHeatmap(sce.all.filt, dims = 1:12, cells = 500, balanced = TRUE) #用于绘制降维后的热图，展示前12个主成分的表达模式，并可以选择展示前500个细胞的数据
#查看每个纬度对数的贡献
ElbowPlot(sce.all.filt, ndims = 50, reduction = "pca") #默认为20
#得到每个PC的p值分布
sce.all.filt <- JackStraw(object = sce.all.filt, num.replicate = 100)
sce.all.filt <- ScoreJackStraw(object = sce.all.filt, dims = 1:20)
pdf(file="02.pcaJackStraw.pdf",width=8,height=6)
JackStrawPlot(object = sce.all.filt, dims = 1:20)
dev.off()

################### 5. 数据聚类 ###################
# 单细胞聚类
# 最近邻算法对细胞进行聚类
sce.all.filt <- FindNeighbors(sce.all.filt, reduction = "pca", dims = 1:20)
# 使用多分辨率模块优化算法，迭代将细胞进行分类
sce.all.filt <- FindClusters(sce.all.filt, resolution = 0.5)
head(sce.all.filt@meta.data)
table(sce.all.filt@meta.data$seurat_clusters)
######5.2 UMAP 可视化
# UMAP聚类可视化
sce.all.filt <- RunUMAP(sce.all.filt, reduction = "pca",
                    dims = 1:20, verbose = TRUE)
DimPlot(sce.all.filt, reduction = "umap", 
        label = TRUE, pt.size = 1.5)
head(sce.all.filt@reductions$umap@cell.embeddings)
#####5.3 tSNE 可视化
# tSNE聚类可视化
sce.all.filt <- RunTSNE(sce.all.filt, reduction = "pca", dims = 1:20)

DimPlot(sce.all.filt, reduction = "tsne", 
        label = TRUE, pt.size = 1.5)
# 自定义颜色， 聚类靠前颜色显著
# DimPlot(sce.all.filt, reduction = "tsne", cols = c("red", "blue", 3:10),
#         label = TRUE, pt.size = 1.5)
head(sce.all.filt@reductions$tsne@cell.embeddings)
saveRDS(sce.all.filt, file = "sce.all_int.rds")


#Harmony去批次
library(harmony)
seuratObj <- RunHarmony(sce.all.filt, "orig.ident")
names(seuratObj@reductions)
#使用UMAP/TSNE可视化Harmony去批次效果
seuratObj <- RunUMAP(seuratObj,  dims = 1:20, 
                     reduction = "harmony")
DimPlot(seuratObj,reduction = "umap",label=F ) 
seuratObj <- RunTSNE(seuratObj, dims = 1:20, 
                     reduction = "harmony")

DimPlot(seuratObj,reduction = "tsne",label=F ) 

#3 细胞聚类
#先运行FindNeighbors
sce.all.filt=seuratObj

sce.all.filt <- FindNeighbors(sce.all.filt, reduction = "harmony",
                              dims = 1:20) 
sce.all.filt.all=sce.all.filt
#再运行FindClusters
#设置不同的分辨率，观察分群效果(选择哪一个？)
for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5,0.8,1)) {
  sce.all.filt.all=FindClusters(sce.all.filt.all, #graph.name = "CCA_snn", 
                                resolution = res, algorithm = 1)
}
colnames(sce.all.filt.all@meta.data)
apply(sce.all.filt.all@meta.data[,grep("RNA_snn",colnames(sce.all.filt.all@meta.data))],2,table)

p1_dim=plot_grid(ncol = 3, DimPlot(sce.all.filt.all, reduction = "umap", group.by = "RNA_snn_res.0.01") + 
                   ggtitle("louvain_0.01"), DimPlot(sce.all.filt.all, reduction = "umap", group.by = "RNA_snn_res.0.1") + 
                   ggtitle("louvain_0.1"), DimPlot(sce.all.filt.all, reduction = "umap", group.by = "RNA_snn_res.0.2") + 
                   ggtitle("louvain_0.2"))
#ggsave(plot=p1_dim, filename="Dimplot_diff_resolution_low.pdf",width = 14)
ggsave(plot = p1_dim, filename = "Dimplot_diff_resolution_low.pdf", width = 14, height = 8)
p1_dim=plot_grid(ncol = 3, DimPlot(sce.all.filt.all, reduction = "umap", group.by = "RNA_snn_res.0.8") + 
                   ggtitle("louvain_0.8"), DimPlot(sce.all.filt.all, reduction = "umap", group.by = "RNA_snn_res.1") + 
                   ggtitle("louvain_1"), DimPlot(sce.all.filt.all, reduction = "umap", group.by = "RNA_snn_res.0.3") + 
                   ggtitle("louvain_0.3"))
#ggsave(plot=p1_dim, filename="Dimplot_diff_resolution_high.pdf",width = 18)
ggsave(plot = p1_dim, filename = "Dimplot_diff_resolution_high.pdf", width = 18, height = 10)
p2_tree=clustree(sce.all.filt.all@meta.data, prefix = "RNA_snn_res.")
#ggsave(plot=p2_tree, filename="Tree_diff_resolution.pdf")
ggsave(plot = p2_tree, filename = "Tree_diff_resolution.pdf", width = 18, height = 10)
table(sce.all.filt.all@active.ident) 
saveRDS(sce.all.filt.all, file = "sce.all_int.rds")
#setwd('../')


setwd("C:\\Users\\Xu416\\Desktop\\data6")
###(三)：细胞分群注释###
#选择0.5分辨率
rm(list=ls())
sce.all.int = readRDS('C:\\Users\\Xu416\\Desktop\\data6\\2-harmony\\sce.all_int.rds') 
sel.clust = "RNA_snn_res.0.5"
sce.all.int <- SetIdent(sce.all.int, value = sel.clust)
table(sce.all.int@active.ident) 
colnames(sce.all.int@meta.data) 

setwd('../')
getwd()
#创建了新的文件夹”3-Celltype“
dir.create("./3-Celltype")
setwd("./3-Celltype")
scRNA=sce.all.int
#绘制气泡图展示maker基因在各分群的表达
genes_to_check = c('KRT8','KRT18', 'CGA' , 'PKIB', 'PHLDA2', #Trophoblast
                   'MFAP4' , 'MDK', 'LUM','MEST' , 'TCF21', 'OLFML3',  'GPC3', #Mesenchymal
                   'EPCAM', 'FOXJ1', 'CAPS', #CEs
                   'PAX8', 'OVGP1',  #NCSE
                   'DCN' , 'COL1A1', 'PDGFRA','COL6A2' , 'COL1A2',  #Fibroblast
                   'POSTN', 'NR2F2','DES','ACTA2','MYH11', 'TAGLN',  #Myofibroblasts
                   'PDGFRB', 'CSPG4',  #SMCs
                   'KDR' ,'PECAM1', 'VWF', #Blood Endothelials
                   'PROX1', 'PDPN', #Lymphatic Endothelials
                   'SDC1', 'JCHAIN', 'MS4A1', #B cells
                   'KLRC1','RUNX3', 'CD3D','NKG7','GNLY', 'PTPRC',  #T/NK cells
                   'FOLR2', 'CD68', 'CD163','ITGAX', 'S100A12', 'S100A8', #Mo_Ma
                   'CST3','FCER1A')#DCs
p = DotPlot(scRNA, features = unique(genes_to_check),
            assay='RNA'  )  + coord_flip()
p

#再次查看分群TSNE图
mycolors <-c('#E64A35','#4DBBD4' ,'#01A187'  ,'#6BD66B','#3C5588'  ,'#F29F80'  ,
             '#8491B6','#91D0C1','#7F5F48','#AF9E85','#4F4FFF','#CE3D33',
             '#739B57','#EFE685','#446983','#BB6239','#5DB1DC','#7F2268','#800202','#D8D8CD'
)

library(RColorBrewer)
mycolors <- rep(mycolors, 3)

####方法一TSNE
tsne =DimPlot(scRNA, reduction = "tsne",cols = mycolors,pt.size = 0.8,
              group.by = "RNA_snn_res.0.5",label = T,label.box = T)
tsne

#结合DotPlot和TSNE图给各细胞群命名
####细胞生物学命名
celltype=data.frame(ClusterID=0:15,
                    celltype= 0:15) 
#这里强烈依赖于生物学背景，看dotplot的基因表达量情况来人工审查单细胞亚群名字。
celltype[celltype$ClusterID %in% c(0,1,3,4,5,6,7,9,10),2]='Trophoblast'
celltype[celltype$ClusterID %in% c(11),2]='T/NK cells'
celltype[celltype$ClusterID %in% c(2,8,14,15),2]='Mo_Ma'
celltype[celltype$ClusterID %in% c(12),2]='Mesenchymal'
celltype[celltype$ClusterID %in% c(13),2]='Fibroblast'
#将细胞类型信息添加到meta.data：
scRNA@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(scRNA@meta.data$RNA_snn_res.0.5 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(scRNA@meta.data$celltype)
#可视化细胞类型
th=theme(axis.text.x = element_text(angle = 45, 
                                    vjust = 0.5, hjust=0.5)) 
library(patchwork)


celltype_tsne =DimPlot(scRNA, reduction = "tsne",cols = mycolors,pt.size = 1,
                       group.by = "celltype",label = T)
celltype_tsne
celltype_umap =DimPlot(scRNA, reduction = "umap",cols = mycolors,pt.size = 1,
                       group.by = "celltype",label = T)
celltype_umap
#综合看样本类型、患者、细胞类型TSNE图
sample_tsne =DimPlot(scRNA, reduction = "tsne",cols = mycolors,pt.size = 0.2,
                     group.by = "sample") 
patient_tsne =DimPlot(scRNA, reduction = "tsne",cols = mycolors,pt.size = 0.2,
                      group.by = "patient") 
sample_tsne + patient_tsne+celltype_tsne 


####方法二UMAP
umap =DimPlot(scRNA, reduction = "umap",cols = mycolors,pt.size = 0.8,
              group.by = "RNA_snn_res.0.5",label = T,label.box = T) 
umap

sample_umap =DimPlot(scRNA, reduction = "umap",cols = mycolors,pt.size = 0.8,
                     group.by = "sample") 
sample_umap
umap+sample_umap

saveRDS(scRNA, "sce_celltype.rds")



####(四)：细胞分群可视化###
setwd("C:\\Users\\Xu416\\Desktop\\data6")
dir.create("4-plot")
setwd('4-plot/')
sce.all=readRDS( "../3-Celltype/sce_celltype.rds")
sce.all
###2.1 Featureplot可视化基因
#EPCAM,NKG7,LYZ,CD79A,CLDN5,DCN
#marker <- c('EPCAM','NKG7','LYZ','CD79A','CLDN5','DCN')
marker <- c('TP63','ITGA6','PEG10','PHLDA2','XAGE3','CDKN1C','KRT7','VGLL1','KRT18','XAGE2')
gene = marker
FeaturePlot(sce.all,features = marker,cols = c("lightgrey" ,"#4C4CB7"),ncol=5,raster=FALSE)
ggsave('FeaturePlot_marker.pdf',width = 12,height = 8)
####上图各个基因的颜色阈值范围不一致，我们可以调参，使其保持一致
p1 <- FeaturePlot(sce.all, features = marker, combine = FALSE,ncol=3,raster=FALSE )
#colours = c('lightgrey', "#DE1F1F")
fix.sc <- scale_color_gradientn( colours = c('lightgrey', "#4C4CB7"),  limits = c(0, 6))
#+NoLegend()+NoAxes()
p2 <- lapply(p1, function (x) x + fix.sc)
CombinePlots(p2)
#批量画基因,注意图例的范围不同：
FeaturePlot(sce.all, features =marker, 
            cols = c("lightgrey", '#4C4CB7'),
            ncol = 3 ) & NoLegend() & NoAxes() & theme(
              panel.border = element_rect(color = "black", size = 1)
            )
ggsave('FeaturePlot_marker1.pdf',width = 12,height = 8)
# 绘制小提琴图，使用调整后的celltype因子水平

VlnPlot(sce.all, features = marker, group.by = "celltype", ncol = 5)
ggsave('markers_top10_VlnPlot.pdf', width = 10, height = 7)
# 绘制山峦图
RidgePlot(sce.all, features = marker, group.by = "celltype", ncol = 5)
ggsave('markers_top10_RidgePlot.pdf', width = 10, height = 7)





#Featureplot还可以把两个基因画在同一个图中,看右上角可以发现黄色越深的地方两个基因叠加越多
FeaturePlot(sce.all, features = c('S100A9','S100A8'),
            cols = c("lightgrey", "green", "orange"),
            blend=T,blend.threshold=0)
#结合ggplot函数，我们还可以把三个基因画在同一个图中
#提取tsne坐标，并提取基因表达数据并与tsne坐标合并：
tsne_df <- as.data.frame(sce.all@reductions$umap@cell.embeddings)
tsne_df$cluster <- as.factor(sce.all$celltype)
head(tsne_df)
gene_df <- as.data.frame(GetAssayData(object = sce.all, slot = "data")[c('S100A9','S100A8','CXCL8'), ])
#ggplot绘制图形
library(ggnewscale)
merged_df <- merge(t(gene_df), tsne_df, by = 0, all = TRUE)
head(merged_df)
colnames(merged_df)
ggplot(merged_df, vars = c("umap_1", "umap_2", 'S100A9','S100A8','CXCL8'), aes(x = umap_1, y = umap_2, colour = S100A9)) +
  geom_point(size=0.3, alpha=1) +
  scale_colour_gradientn(colours = c("lightgrey", "green"), limits = c(0, 0.3), oob = scales::squish) +
  new_scale_color() +
  geom_point(aes(colour = S100A8), size=0.3, alpha=0.7) +
  scale_colour_gradientn(colours = c("lightgrey", "blue"), limits = c(0.1, 0.2), oob = scales::squish) +
  new_scale_color() +
  geom_point(aes(colour = CXCL8), size=0.3, alpha=0.1) +
  scale_colour_gradientn(colours = c("lightgrey", "red"), limits = c(0, 0.3), oob = scales::squish)+
  theme_classic()

#2.2 DoHeatmap绘制热图
#利用DoHeatmap函数绘制热图，可以展示不同细胞类型的top5 maker：
Idents(sce.all)
table(sce.all$celltype)
Idents(sce.all) = sce.all$celltype
sce1 = sce.all[, sce.all$celltype %in% c( 'Trophoblast', 'T/NK cells','Mo_Ma', 'Mesenchymal','Fibroblast' )]

if (!file.exists('sce.markers.csv')) {
  sce.markers <- FindAllMarkers(object = sce1, only.pos = TRUE, 
                                min.pct = 0.25, 
                                thresh.use = 0.25)
  write.csv(sce.markers,file='sce.markers.csv')
} else {
  
  sce.markers = read.csv('sce.markers.csv',row.names = 1)
}
library(dplyr) 
top5 <- sce.markers%>% group_by(cluster) %>% top_n(5, avg_log2FC)
#为了防止数据量太大不好出图，这里在每个亚群提取出来100个：
sce.Scale <- ScaleData(subset(sce1,downsample=100),
                       features = top5$gene )
DoHeatmap(sce.Scale,
          features = top5$gene ,
          # group.by = "celltype",
          assay = 'RNA', label = T)+
  scale_fill_gradientn(colors = c("white","grey","firebrick3"))

ggsave('markers_heatmap.pdf',width = 10,height = 7)

#2.3 DotPlot绘制气泡图
top5_dotplot <- DotPlot(sce.all, features = top5$gene)+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))
top5_dotplot
ggsave('markers_top5_dotplot.pdf',width = 10,height = 7)



# 提取每个亚群的TOP 10标记基因
top10 <- sce.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)

# 绘制热图
sce.Scale <- ScaleData(sce1, features = top10$gene)
DoHeatmap(sce.Scale,
          features = top10$gene,
          assay = 'RNA', label = T) +
  scale_fill_gradientn(colors = c("white", "grey", "firebrick3"))

ggsave('markers_top10_heatmap.pdf', width = 10, height = 7)
# 绘制气泡图
top10_dotplot <- DotPlot(sce.all, features = top10$gene) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))
top10_dotplot
ggsave('markers_top10_dotplot.pdf', width = 10, height = 7)
# top10_genes是包含TOP 10标记基因名称的向量
top10_genes <- c("PAGE4", "CGA", "PEG10", "PHLDA2", "XAGE3", "CDKN1C", "KRT7", "VGLL1", "KRT18", "XAGE2")
#PAGE4,CGA,PEG10,PHLDA2,XAGE3, CDKN1C, KRT7, VGLL1, KRT18, XAGE2
# 绘制小提琴图，使用调整后的celltype因子水平
#VlnPlot(sce1, features = top10_genes, group.by = "celltype", ncol = 5)
# 绘制小提琴图（去除黑点）
VlnPlot(sce.all, features = top10_genes, group.by = "celltype", ncol = 5, pt.size = 0)
ggsave('markers_top10_VlnPlot.pdf', width = 10, height = 7)
# 绘制山峦图
RidgePlot(sce1, features = top10_genes, group.by = "celltype", ncol = 5)
ggsave('markers_top10_RidgePlot.pdf', width = 10, height = 7)

setwd('../')
getwd()


###(五)：细胞比例###
#2.1 饼图-饼图可以直观展示组内各细胞比例的差异
rm(list=ls())
options(stringsAsFactors = F)
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(plotrix)
library(ggsci)
library(celldex)
library(singleseqgset)
library(devtools)
getwd()
dir.create("5-prop")
setwd('5-prop/')
sce.all=readRDS( "../3-Celltype/sce_celltype.rds")
sce.all
#绘制饼图
head(sce.all@meta.data)
table(sce.all$celltype) #创建了一个频率表，显示sce.all对象中celltype列的不同值及其出现次数
mynames <-   table(sce.all$celltype) %>% names() #将celltype列的不同值存储在变量mynames中。使用管道操作符%>%，它允许将前一个函数的输出作为下一个函数的输入
myratio <-  table(sce.all$celltype) %>% as.numeric() #将celltype列的频率表转换为数值向量，并存储在变量myratio中
pielabel <- paste0(mynames," (", round(myratio/sum(myratio)*100,2), "%)")

cols <-c('#E64A35','#4DBBD4' ,'#01A187','#6BD66B','#3C5588'  ,'#F29F80'  ,
         '#8491B6','#91D0C1','#7F5F48','#AF9E85','#4F4FFF','#CE3D33',
         '#739B57','#EFE685','#446983','#BB6239','#5DB1DC','#7F2268','#800202','#D8D8CD'
)
pie(myratio, labels=pielabel,
    radius = 1.0,clockwise=T,
    main = "celltype",col = cols)
# 首先，指定PDF文件的保存路径和文件名
pdf("Cell_Proportion_3D_Pie_Chart.pdf", width = 8, height = 6)
#绘制3D饼图
pie3D(myratio,labels = pielabel,explode = 0.1, 
      main = "Cell Proption",
      height = 0.3,
      labelcex  = 1)
# 最后，关闭PDF设备，这将保存并关闭文件
dev.off()
setwd('../')

###(六)：滋养层细胞亚群细分###
rm(list=ls())
options(stringsAsFactors = F)
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
dir.create("6-T")
setwd("6-T/") 
getwd()
set.seed(12345)
sce.all=readRDS( "../3-Celltype/sce_celltype.rds")
table(sce.all$celltype)
sce1 = sce.all[, sce.all$celltype %in% c( 'Trophoblast' )]


##2.2 降维分群聚类
LayerData(sce1, assay = "RNA", layer = "counts")
sce1 <- JoinLayers(sce1)
as.data.frame(sce1@assays$RNA$counts[1:10, 1:2])
head(sce1@meta.data, 10)
table(sce1$orig.ident) 
sce = sce1
sce <- NormalizeData(sce, normalization.method =  "LogNormalize",  
                     scale.factor = 1e4)
GetAssay(sce,assay = "RNA")
sce <- FindVariableFeatures(sce, 
                            selection.method = "vst", nfeatures = 2000)  
sce <- ScaleData(sce) 
sce <- RunPCA(object = sce, pc.genes = VariableFeatures(sce)) 

DimHeatmap(sce, dims = 1:12, cells = 100, balanced = TRUE)
ElbowPlot(sce) 
ggsave('ElbowPlot.pdf', width = 10, height = 7)
#选取0.1分辨率，对T细胞进行聚类分群
sce <- FindNeighbors(sce, dims = 1:15)
sce <- FindClusters(sce, resolution = 0.1)
table(sce@meta.data$RNA_snn_res.0.1)  
set.seed(321)
sce <- RunTSNE(object = sce, dims = 1:15, do.fast = TRUE)

sce <- RunUMAP(object = sce, dims = 1:5, do.fast = TRUE)

mycolors <-c('#E64A35','#4DBBD4' ,'#01A187'  ,'#3C5588'  ,'#F29F80'  ,
             '#8491B6','#91D0C1','#7F5F48','#AF9E85','#4F4FFF','#CE3D33',
             '#739B57','#EFE685','#446983','#BB6239','#5DB1DC','#7F2268','#6BD66B','#800202','#D8D8CD','pink'
)

p = DimPlot(sce,reduction = "tsne",label=T,cols = mycolors) 
p
ggsave('tsne.pdf', width = 10, height = 7)
p1 = DimPlot(sce,reduction = "umap",label=T,cols = mycolors) 
p1
ggsave('UMAP.pdf', width = 10, height = 7)

scRNA=sce
#绘制气泡图展示maker基因在各分群的表达
genes_to_check = c('TP63', 'ITGA6','KRT7', 'PARP1',   # VCTs
                   'ERVFRD-1', 'LGALS16', 'GDF15','TBX3',   # SCTs
                   'HLA-G', 'ASCL2', 'PLAC8' # EVTs
                   
)
p = DotPlot(scRNA, features = unique(genes_to_check),
            assay='RNA'  )  + coord_flip()
p
ggsave('marker_DotPlot.pdf', width = 10, height = 7)


#再次查看分群TSNE图
mycolors <-c('#E64A35','#4DBBD4' ,'#01A187'  ,'#6BD66B','#3C5588'  ,'#F29F80'  ,
             '#8491B6','#91D0C1','#7F5F48','#AF9E85','#4F4FFF','#CE3D33',
             '#739B57','#EFE685','#446983','#BB6239','#5DB1DC','#7F2268','#800202','#D8D8CD'
)

library(RColorBrewer)
mycolors <- rep(mycolors, 3)

####方法一TSNE
tsne =DimPlot(scRNA, reduction = "tsne",cols = mycolors,pt.size = 0.8,
              group.by = "RNA_snn_res.0.1",label = T,label.box = T)
tsne
ggsave('RNA_snn_res.0.1_tsne.pdf', width = 10, height = 7)
#结合DotPlot和TSNE图给各细胞群命名
####细胞生物学命名
celltype=data.frame(ClusterID=0:15,
                    celltype= 0:15) 
#这里强烈依赖于生物学背景，看dotplot的基因表达量情况来人工审查单细胞亚群名字。
celltype[celltype$ClusterID %in% c(1,5),2]='VCTs'
celltype[celltype$ClusterID %in% c(3),2]='SCTs'
celltype[celltype$ClusterID %in% c(0,2,4),2]='EVTs'
#将细胞类型信息添加到meta.data：
scRNA@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(scRNA@meta.data$RNA_snn_res.0.1 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(scRNA@meta.data$celltype)
#可视化细胞类型
th=theme(axis.text.x = element_text(angle = 45, 
                                    vjust = 0.5, hjust=0.5)) 
library(patchwork)


celltype_tsne =DimPlot(scRNA, reduction = "tsne",cols = mycolors,pt.size = 1,
                       group.by = "celltype",label = T)
celltype_tsne
ggsave('celltype_tsne.pdf', width = 10, height = 7)

celltype_umap =DimPlot(scRNA, reduction = "umap",cols = mycolors,pt.size = 1,
                       group.by = "celltype",label = T)
celltype_umap
ggsave('celltype_umap.pdf', width = 10, height = 7)
saveRDS(scRNA, "sce_celltype.rds")

####细胞分群可视化###
sce.all=scRNA
###2.1 Featureplot可视化基因
#EPCAM,NKG7,LYZ,CD79A,CLDN5,DCN
marker <- c('EPCAM','NKG7','LYZ','CD79A','CLDN5','DCN')
gene = marker
FeaturePlot(sce.all,features = marker,cols = c("lightgrey" ,"#DE1F1F"),ncol=3,raster=FALSE)
ggsave('FeaturePlot_marker.pdf',width = 12,height = 8)
####上图各个基因的颜色阈值范围不一致，我们可以调参，使其保持一致
p1 <- FeaturePlot(sce.all, features = marker, combine = FALSE,ncol=3,raster=FALSE )
#colours = c('lightgrey', "#DE1F1F")
fix.sc <- scale_color_gradientn( colours = c('lightgrey', "#DE1F1F"),  limits = c(0, 6))
#+NoLegend()+NoAxes()
p2 <- lapply(p1, function (x) x + fix.sc)
CombinePlots(p2)
#批量画基因,注意图例的范围不同：
FeaturePlot(sce.all, features =marker, 
            cols = c("lightgrey", 'red'),
            ncol = 3 ) & NoLegend() & NoAxes() & theme(
              panel.border = element_rect(color = "black", size = 1)
            )
#Featureplot还可以把两个基因画在同一个图中,看右上角可以发现黄色越深的地方两个基因叠加越多
FeaturePlot(sce.all, features = c('S100A9','S100A8'),
            cols = c("lightgrey", "green", "orange"),
            blend=T,blend.threshold=0)
#结合ggplot函数，我们还可以把三个基因画在同一个图中
#提取tsne坐标，并提取基因表达数据并与tsne坐标合并：
tsne_df <- as.data.frame(sce.all@reductions$umap@cell.embeddings)
tsne_df$cluster <- as.factor(sce.all$celltype)
head(tsne_df)
gene_df <- as.data.frame(GetAssayData(object = sce.all, slot = "data")[c('S100A9','S100A8','CXCL8'), ])
#ggplot绘制图形
library(ggnewscale)
merged_df <- merge(t(gene_df), tsne_df, by = 0, all = TRUE)
head(merged_df)
colnames(merged_df)
ggplot(merged_df, vars = c("umap_1", "umap_2", 'S100A9','S100A8','CXCL8'), aes(x = umap_1, y = umap_2, colour = S100A9)) +
  geom_point(size=0.3, alpha=1) +
  scale_colour_gradientn(colours = c("lightgrey", "green"), limits = c(0, 0.3), oob = scales::squish) +
  new_scale_color() +
  geom_point(aes(colour = S100A8), size=0.3, alpha=0.7) +
  scale_colour_gradientn(colours = c("lightgrey", "blue"), limits = c(0.1, 0.2), oob = scales::squish) +
  new_scale_color() +
  geom_point(aes(colour = CXCL8), size=0.3, alpha=0.1) +
  scale_colour_gradientn(colours = c("lightgrey", "red"), limits = c(0, 0.3), oob = scales::squish)+
  theme_classic()

#2.2 DoHeatmap绘制热图
#利用DoHeatmap函数绘制热图，可以展示不同细胞类型的top5 maker：
Idents(sce.all)
table(sce.all$celltype)
Idents(sce.all) = sce.all$celltype
sce1 = sce.all[, sce.all$celltype %in% c( 'VCTs', 'SCTs','EVTs' )]
# 预处理步骤（确保细胞类型顺序）
sce1$celltype <- factor(sce1$celltype, levels = c("VCTs", "SCTs", "EVTs"))
Idents(sce1) <- "celltype"

# 绘制dotplot
DotPlot(sce1,
        features = c('TP63', 'ITGA6','KRT7', 'PARP1', 'ERVFRD-1', 'LGALS16', 'GDF15','TBX3', 'HLA-G', 'ASCL2', 'PLAC8'),
        cols = c("lightgrey", "blue"),   # 颜色渐变设置
        dot.scale = 6,                   # 点的大小缩放
        cluster.idents = FALSE) +        # 关闭聚类保持原始顺序
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # x轴标签倾斜45度
  labs(x = "Genes", y = "Cell Types")    # 坐标轴标签
ggsave('markers_dotplot.pdf', width = 10, height = 7)




if (!file.exists('sce.markers.csv')) {
  sce.markers <- FindAllMarkers(object = sce1, only.pos = TRUE, 
                                min.pct = 0.25, 
                                thresh.use = 0.25)
  write.csv(sce.markers,file='sce.markers.csv')
} else {
  
  sce.markers = read.csv('sce.markers.csv',row.names = 1)
}
library(dplyr) 
# 提取每个亚群的TOP 10标记基因
top10 <- sce.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)

# 绘制热图
sce.Scale <- ScaleData(sce1, features = top10$gene)
DoHeatmap(sce.Scale,
          features = top10$gene,
          assay = 'RNA', label = T) +
  scale_fill_gradientn(colors = c("white", "grey", "firebrick3"))

ggsave('markers_top10_heatmap.pdf', width = 10, height = 7)

# 绘制气泡图
top10_dotplot <- DotPlot(sce.all, features = top10$gene) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))
top10_dotplot
ggsave('markers_top10_dotplot.pdf', width = 10, height = 7)

#VCTs
top10_genes <- c("PAGE4", "PEG10", "BEX3", "SMAGP", "XAGE3", "DUSP9", "ISYNA1", "VGLL1", "EFEMP1", "PEG3")
# 确保celltype列是因子类型，并设置因子水平以指定顺序
sce1$celltype <- factor(sce1$celltype, levels = c("VCTs", "SCTs", "EVTs"))
# 绘制小提琴图（去除黑点）
VlnPlot(sce1, features = top10_genes, group.by = "celltype", ncol = 5, pt.size = 0)
# 绘制小提琴图，使用调整后的celltype因子水平
#VlnPlot(sce1, features = top10_genes, group.by = "celltype", ncol = 5)
ggsave('markers_top10_VlnPlot_VCTs.pdf', width = 10, height = 7)
# 绘制山峦图
RidgePlot(sce1, features = top10_genes, group.by = "celltype", ncol = 5)
ggsave('markers_top10_RidgePlot_VCTs.pdf', width = 10, height = 7)

#SCTs
top10_genes <- c("CGA", "INSL4", "GADD45G", "GDF15", "CYP19A1", "HSD17B1", "TFPI2", "KISS1", "CSH2", "CSH1")
# 确保celltype列是因子类型，并设置因子水平以指定顺序
sce1$celltype <- factor(sce1$celltype, levels = c("VCTs", "SCTs", "EVTs"))
# 绘制小提琴图（去除黑点）
VlnPlot(sce1, features = top10_genes, group.by = "celltype", ncol = 5, pt.size = 0)
# 绘制小提琴图，使用调整后的celltype因子水平
#VlnPlot(sce1, features = top10_genes, group.by = "celltype", ncol = 5)
ggsave('markers_top10_VlnPlot_SCTs.pdf', width = 10, height = 7)
# 绘制山峦图
RidgePlot(sce1, features = top10_genes, group.by = "celltype", ncol = 5)
ggsave('markers_top10_RidgePlot_SCTs.pdf', width = 10, height = 7)

#EVTs
top10_genes <- c("FSTL3", "TPM1", "HPGD", "LAIR2", "IGFBP3", "FN1", "HLA-G", "NOTUM", "HTRA4", "SERPINE2")
# 确保celltype列是因子类型，并设置因子水平以指定顺序
sce1$celltype <- factor(sce1$celltype, levels = c("VCTs", "SCTs", "EVTs"))
# 绘制小提琴图（去除黑点）
VlnPlot(sce1, features = top10_genes, group.by = "celltype", ncol = 5, pt.size = 0)
# 绘制小提琴图，使用调整后的celltype因子水平
#VlnPlot(sce1, features = top10_genes, group.by = "celltype", ncol = 5)
ggsave('markers_top10_VlnPlot_EVTs.pdf', width = 10, height = 7)
# 绘制山峦图
RidgePlot(sce1, features = top10_genes, group.by = "celltype", ncol = 5)
ggsave('markers_top10_RidgePlot_EVTs.pdf', width = 10, height = 7)
saveRDS(sce, file = "T_sce_celltype.rds")
setwd("../")
