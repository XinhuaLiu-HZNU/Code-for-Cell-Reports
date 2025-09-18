### Single cell RNA-seq is mainly used for RNA quantification at the single cell level followed by cell clustering and annotation to compare the composition difference of cells among different groups.
### This pipeline contains the all analysis steps except the mapping and quantifiction by cellranger or other preprocessing softwares. Seurat and SingleR is the two main packages used in this pipeline.
### The analysis in this pipeline mainly including:
### 1. Filtering of low-quality cells and genes, including cells contain abnormally high MT gene, low gene number, and genes that expressed at low number cells.
### 2. Cell clustering. Determination of optimal resolution for cell clustering is critical (clustree).
### 3. Cell annotation. Combing SingleR and manual annotation result should be appreciated. The most heavy work is marker collection.
### 4. Comparison, GSEA.

########################## The following is the main body of this pipeline ####################################

### Load the required packages.
lapply(list("Seurat","SeuratObject","dplyr","reshape2","ggplot2","patchwork","SingleR"),require,character.only=TRUE)

### Quanlity control and clustering.
### (1) Rename the samples
scdir <- 'your_single_cell_data_path'
scfile <- list(scdir,pattern='*.tarx',recursive=TRUE,full.names=TRUE,include.dirs=TRUE)
pop_ <- function(x,sep=',',pos=0){
    l=unlist(strsplit(x,split=sep))
    if(pos==0){
        pos=length(l)
    }
    return(l[pos])
}
scsample <- unlist(lapply(scfile,pop_,sep='/',pos=0))
scsample <- unlist(lapply(scsample,pop_,sep='_',pos=1))

### (2) Read, integrate and normalize the data.
scdata.list <- list()
for(i in 1:length(scfile)){
    scobj = Read10X(scfile[i])  ####要关注features文件有几列，只有一列时要设置gene.column=1
    scobj = CreateSeuratObject(counts=scobj,min.cells=10,min.features=200,project=scsample[i])
    scdata.list[[scsample[i]]] = scobj
}
scdata.list

###features <- SelectIntegrationFeatures(object.list=scdata.list.int)
for (i in 1:6) {
    scdata.list[[i]] <- SCTransform(scdata.list[[i]]) %>%
    RunPCA()
}

anchors <- FindIntegrationAnchors(object.list=scdata.list)
scdata.list.int <- IntegrateData(anchorset=anchors)
scdata.list.int

### (3) Cell clustering
DefaultAssay(scdata.list.int) <- "integrated"
scdata.list.int <- ScaleData(scdata.list.int,verbose=FALSE)
scdata.list.int <- RunPCA(scdata.list.int,npcs=50,verbose=FALSE)
options(repr.plot.width=10,repr.plot.height=6)
ElbowPlot(scdata.list.int,ndims=50)

scdata.list.int <- RunUMAP(scdata.list.int,reduction='pca',dims=1:30)
scdata.list.int <- FindNeighbors(scdata.list.int,reduction='pca',dims=1:30)
scdata.list.int <- FindClusters(scdata.list.int,resolution=c(seq(0.1,2,0.1),seq(0.01,0.1,0.01)))

options(repr.plot.width=14,repr.plot.height=22)
library(clustree)
clustree(scdata.list.int@meta.data,prefix='integrated_snn_res.')

### (4) Annotation based on SingleR
library(celldex)
hpca <- celldex::HumanPrimaryCellAtlasData()
cell.ct <- GetAssayData(object=scdata.list.int[['RNA']],slot='counts')
dim(cell.ct)
pred_ident <- SingleR(test=cell.ct,ref=hpca,assay.type.test=1,labels=hpca$label.main)
scdata.list.int[['SingleR.labels']] <- pred_ident$labels

### (5) Cell composition plot based on the selected varaible
library(scales)
library(ggsci)
options(repr.plot.height=8,repr.plot.width=12)
md <- scdata.list.int@meta.data
md %>% group_by(orig.ident,SingleR.labels) %>% 
    summarise(num=n(),.groups='drop') %>%
    arrange(orig.ident,num) %>% 
    ggplot(aes(x=orig.ident,y=num,fill=SingleR.labels))+
    geom_bar(position='fill',stat='identity')+
    scale_fill_ucscgb()

### (6) Find markers for every cell cluster. Note: the optimal resolution determines the cluster numbers
options(repr.plot.height=10,repr.plot.width=15)
AllMarkers <- FindAllMarkers(scdata.list.int,group.by='integrated_snn_res.1.4',min.pct=0.1,logfc.threshold=0.25)
AllMarkers %>% group_by(cluster) %>% slice_max(n=10,avg_log2FC) -> Markers.top10
options(repr.plot.width=10,repr.plot.height=15)
DotPlot(scdata.list.int,features=unique(Markers.top10$gene))+theme(axis.text.x=element_text(vjust=0.5,hjust=0.5,angle=90))

### (7) Find differential genes for the interested variables in each cell type and conducted functional enrichment analysis.
library(pheatmap)
library(clusterProfiler)
library(EnsDb.Hsapiens.v86)
mysamples <- unique(md$orig.ident)
scdata.list <- list()
for(i in 1:length(mysamples)){
    scdata.list[[mysamples[i]]] = subset(scdata.list.int,orig.ident==mysamples[i])
}
myDEGs <- lapply(scdata.list,FindMarkers,min.pct=0.1,logfc.threshold=0.25,ident.1='acute',ident.2='recover')
myDEGs <- lapply(myDEGs,FUN=function(x){mutate(x,gene=rownames(x))})
do.call(rbind,myDEGs) %>% group_by(gene) %>%
    summarise(num=n(),.groups='drop') -> DiffNum
DEG3 <- DiffNum$gene[DiffNum$num>=3] %>% as.character()
DEG3LFC <- data.frame(row.names=DEG3,FNQ=myDEGs[['Fu.naiqiang']][DEG3,'avg_log2FC'],
    WJF=myDEGs[['Wang.jifang']][DEG3,'avg_log2FC'],
    YHJ=myDEGs[['Yang.huanjun']][DEG3,'avg_log2FC'],
    XXH=myDEGs[['Xu.xinhai']][DEG3,'avg_log2FC'],
    WGR=myDEGs[['Wang.guorong']][DEG3,'avg_log2FC'],
    WY=myDEGs[['Wang.yan']][DEG3,'avg_log2FC'])

pheatmap(DEG3LFC,color=colorRampPalette(c("blue","white","red"))(50))
DEG3 <- select(EnsDb.Hsapiens.v86,keys=DEG3,keytype='SYMBOL',columns='ENTREZID')
ek.DEG3 <- enrichKEGG(as.character(DEG3$ENTREZID))
dotplot(ek.DEG3,showCategory=20)

### (8) Gene Set Enrichment Analysis
library(fgsea)
library(msigdbr)
msigdbr(species="Homo sapiens",category="C2",subcategory="CP:KEGG") %>% 
    split(x=.$gene_symbol,f=.$gs_name) -> fgst.KEGG

DEG.agg <- FindMarkers(scdata.list.int,ident.1='acute',ident.2='recover',min.pct=0.1,logfc.threshold=0.25)
DEG.agg <- DEG.agg[order(DEG.agg$avg_log2FC,decreasing=T),]
FCgenelist <- DEG.agg$avg_log2FC
names(FCgenelist) <- rownames(DEG.agg)
DEG.agg.GSEA <- fgsea(pathways=fgst.KEGG,
    stats=FCgenelist,
    minSize=15,
    maxSize=500,
    nperm=10000)
topPathwayUp <- DEG.agg.GSEA[ES>0][head(order(pval),n=10),pathway]
topPathwayDown <- DEG.agg.GSEA[ES<0][head(order(pval),n=10),pathway]
topPathways <- c(topPathwayUp,topPathwayDown)
plotGseaTable(fgst.KEGG[topPathways],FCgenelist,DEG.agg.GSEA,gseaParam=0.5)
plotEnrichment(fgst.KEGG[['YourPathway']],FCgenelist)