# Single cell analysis pipeline
library(Seurat)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(ggsci)
library(gprofiler2)
library(DOSE)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(dplyr)


#______________________________________________________________________________________________________________________________#
# 1. Data import ####

old6.data <- Read10X(data.dir ="./raw_matrix/batch_6_redo/old6_redo/outs/filtered_feature_bc_matrix/")
oldex6.data <- Read10X(data.dir ="./raw_matrix/batch_6_redo/oldEx6_redo/outs/filtered_feature_bc_matrix/")

old6 <- CreateSeuratObject(counts  = old6.data, project = "old6",min.cells = 3)
oldex6 <- CreateSeuratObject(counts  = oldex6.data, project = "oldex6",min.cells = 3)

old6 <- RenameCells(old6, add.cell.id = "B6") # Adding the batch label to cell names incase duplicated barcode
oldex6 <- RenameCells(oldex6, add.cell.id = "B6") # Adding the batch label to cell names incase duplicated barcode

#______________________________________________________________________________________________________________________________#
# 2. Basic process
ovy6 <- merge(x=old6,y=oldex6)
ovy6 <- NormalizeData(ovy6)
ovy6 <- FindVariableFeatures(ovy6,nfeatures = 2000)

#______________________________________________________________________________________________________________________________#
# 3. QC and finish dimension reduction
mito.features <- grep(pattern = "^mt-", x = rownames(x = ovy6), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = ovy6, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = ovy6, slot = 'counts'))

ovy6[['percent.mito']] <- percent.mito # Adding to meta data

ovy6 <- subset(x = ovy6, subset =quantile(ovy6$nFeature_RNA,0.95,na.rm = T) > nFeature_RNA & nFeature_RNA > quantile(ovy6$nFeature_RNA,0.05,na.rm = T) & quantile(ovy6$nCount_RNA,0.95,na.rm = T)  > nCount_RNA & nCount_RNA > quantile(ovy6$nCount_RNA,0.05,na.rm = T) & percent.mito < 0.20 )

ovy6 <- ScaleData(ovy6,vars.to.regress = c("nCount_RNA"))
ovy6 <- RunPCA(object = ovy6, verbose = FALSE)
ovy6 <- RunTSNE(ovy6,dims = 1:30)
ovy6 <- FindNeighbors(object = ovy6, dims = 1:30)
ovy6 <- FindClusters(ovy6,resolution = 1.6)

#______________________________________________________________________________________________________________________________#
# 4. Cluster identification
Idents(ovy6) <- "RNA_snn_res.1.6"
pdf("./Batch6_Cluster_markers.pdf")
VlnPlot(object = ovy6, features = c("Kcnj8"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ovy6, features = c("Acta2"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ovy6, features = c("Cldn5"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ovy6, features = c("Ctss"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ovy6, features = c("Ntsr2"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ovy6, features = c("Syt1"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ovy6, features = c("Cdk1"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ovy6, features = c("Pdgfra"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ovy6, features = c("Cldn11"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ovy6, features = c("Sox11"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ovy6, features = c("Ccdc153"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ovy6, features = c("Ttr"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ovy6, features = c("Alas2"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ovy6, features = c("Pf4"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ovy6, features = c("Rax"),pt.size = 0,ncol = 1,y.max = 8)
VlnPlot(object = ovy6, features = c("Plac8"),pt.size = 0,ncol = 1,y.max = 8)
dev.off()

Idents(ovy6, cells = WhichCells(object = ovy6, idents = c("10")))<-"PC"
Idents(ovy6, cells = WhichCells(object = ovy6, idents = c("7","13")))<-"SMC"
Idents(ovy6, cells = WhichCells(object = ovy6, idents = c("0","1","2","4","5","11","17")))<-"EC"
Idents(ovy6, cells = WhichCells(object = ovy6, idents = c("3","14")))<-"MG"
# Idents(ovy6, cells = WhichCells(object = ovy6, idents = c("4","5","6","25")))<-"AC"
Idents(ovy6, cells = WhichCells(object = ovy6, idents = c("18")))<-"OPC"
# Idents(ovy6, cells = WhichCells(object = ovy6, idents = c("32")))<-"NRP"
Idents(ovy6, cells = WhichCells(object = ovy6, idents = c("16")))<-"OLG"
# Idents(ovy6, cells = WhichCells(object = ovy6, idents = c("24")))<-"imNeur"
# Idents(ovy6, cells = WhichCells(object = ovy6, idents = c("34")))<-"mNeur"
# Idents(ovy6, cells = WhichCells(object = ovy6, idents = c("23")))<-"EPC"
Idents(ovy6, cells = WhichCells(object = ovy6, idents = c("9")))<-"CPC"
# Idents(ovy6, cells = WhichCells(object = ovy6, idents = c("30","43")))<-"Hb_EC"
# Idents(ovy6, cells = WhichCells(object = ovy6, idents = c("17")))<-"MAC"
# Idents(ovy6, cells = WhichCells(object = ovy6, idents = c("38")))<-"TNC"
# Idents(ovy6, cells = WhichCells(object = ovy6, idents = c("33")))<-"MNC"
ovy6$Celltype <- Idents(ovy6)

#______________________________________________________________________________________________________________________________#
# 4. Subtype identification

Idents(ovy6) <- "Celltype"
EC <- subset(ovy6,idents = "EC")
EC <- as.SingleCellExperiment(EC)
saveRDS(EC,"./inter_data/EC.data")

EC.type <- readRDS("./cellassign_results/EC_results.rds")
EC.tmp<- as.character(EC.type)
names(EC.tmp)<- names(EC.type)

ovy6$Subtype <- EC.tmp # Send the results back to Seurat object

#______________________________________________________________________________________________________________________________#
# 5. Find the DEGs

Idents(ovy6) <- "Subtype"

A1.data <- subset(ovy6, idents = "A1") # Extract data into individual obj ovy6ts
A2.data <- subset(ovy6, idents = "A2")
V.data <- subset(ovy6, idents = "V")
AV.data <- subset(ovy6, idents = "AV")
Cap.data <- subset(ovy6, idents = "Cap")
VCap.data <- subset(ovy6, idents = "VCap")
Idents(ovy6) <- "Celltype"
EC.data <- subset(ovy6,idents = "EC")

Idents(A1.data) <- "orig.ident" # Set the idents back to orig.ident
Idents(A2.data) <- "orig.ident"
Idents(V.data) <- "orig.ident"
Idents(AV.data) <- "orig.ident"
Idents(Cap.data) <- "orig.ident"
Idents(VCap.data) <- "orig.ident"
Idents(EC.data) <- "orig.ident"

EC.ove.diff <- FindMarkers(object = EC.data, ident.1 = "old6", ident.2 = "oldex6",  assay ="RNA",verbose = T,test.use = "MAST",logfc.threshold =0)
A1.ove.diff <- FindMarkers(object = A1.data, ident.1 = "old6", ident.2 = "oldex6", logfc.threshold =0, assay = "RNA",test.use = "MAST")
A2.ove.diff <- FindMarkers(object = A2.data, ident.1 = "old6", ident.2 = "oldex6", logfc.threshold =0, assay = "RNA",test.use = "MAST")
V.ove.diff <- FindMarkers(object = V.data, ident.1 = "old6", ident.2 = "oldex6", logfc.threshold =0, assay = "RNA",test.use = "MAST")
AV.ove.diff <- FindMarkers(object = AV.data, ident.1 = "old6", ident.2 = "oldex6", logfc.threshold =0, assay = "RNA",test.use = "MAST")
Cap.ove.diff <- FindMarkers(object = Cap.data, ident.1 = "old6", ident.2 = "oldex6", logfc.threshold =0, assay = "RNA",test.use = "MAST")
VCap.ove.diff <- FindMarkers(object = VCap.data, ident.1 = "old6", ident.2 = "oldex6", logfc.threshold =0, assay = "RNA",test.use = "MAST")
EC.ove6.dif <- list(EC=EC.ove.diff, A1=A1.ove.diff,A2=A2.ove.diff,Cap=Cap.ove.diff,VCap=VCap.ove.diff,V=V.ove.diff, AV=AV.ove.diff)


# Cell propotion

Idents(ovy6) <- "orig.ident"
ovy.old <- subset(ovy6,cells = WhichCells(object = ovy6,idents = c('old6')))
ovy.young <- subset(ovy6,cells = WhichCells(object = ovy6,idents = c('oldex6')))
subtype.old <- as.data.frame(table(ovy.old$Subtype)) %>% as_tibble() %>% filter(Var1 %in% c("A1","A2","Cap","VCap","V","AV"))
subtype.old <- subtype.old[c(1,2,4,6,5,3),] %>% mutate(prop= round(Freq/sum(.$Freq)*100,3)) %>% mutate(lab.ypos = round(cumsum(prop) - 0.5*prop,2))
subtype.young <- as.data.frame(table(ovy.young$Subtype)) %>% as_tibble() %>% filter(Var1 %in% c("A1","A2","Cap","VCap","V","AV"))
subtype.young <- subtype.young[c(1,2,4,6,5,3),] %>% mutate(prop= round(Freq/sum(.$Freq)*100,3)) %>% mutate(lab.ypos = round(cumsum(prop) - 0.5*prop,2))
subtype.old$Var1 <- c("aEC1","aEC2","CapEC","vCapEC","vEC","avEC")
subtype.young$Var1 <- c("aEC1","aEC2","CapEC","vCapEC","vEC","avEC")
subtype.old$Var1 <- factor(subtype.old$Var1, levels = c("CapEC","vCapEC","vEC","avEC","aEC1","aEC2"))
subtype.young$Var1 <- factor(subtype.young$Var1, levels = c("CapEC","vCapEC","vEC","avEC","aEC1","aEC2"))


# Export the DEG for comparision

Idents(ovy6) <- "Celltype"
EC.data <- subset(ovy6, idents = "EC")
MG.data <- subset(ovy6, idents = "MG")
SMC.data <- subset(ovy6, idents = "SMC")
PC.data <- subset(ovy6, idents = "PC")
OPC.data <- subset(ovy6, idents = "18")
CPC.data <- subset(ovy6, idents = "9")
OLG.data <- subset(ovy6, idents = "16")

Idents(EC.data) <- "orig.ident" # Set the idents back to orig.ident
Idents(MG.data) <- "orig.ident"
Idents(SMC.data) <- "orig.ident"
Idents(PC.data) <- "orig.ident"
Idents(OPC.data) <- "orig.ident"
Idents(CPC.data) <- "orig.ident"
Idents(OLG.data) <- "orig.ident"

EC.ove.diff <- FindMarkers(object = EC.data, ident.1 = "old6", ident.2 = "oldex6",  assay ="RNA",verbose = T,test.use = "MAST",logfc.threshold =0)
MG.ove.diff <- FindMarkers(object = MG.data, ident.1 = "old6", ident.2 = "oldex6",  assay ="RNA",verbose = T,test.use = "MAST",logfc.threshold =0)
SMC.ove.diff <- FindMarkers(object = SMC.data, ident.1 = "old6", ident.2 = "oldex6",  assay ="RNA",verbose = T,test.use = "MAST",logfc.threshold =0)
PC.ove.diff <- FindMarkers(object = PC.data, ident.1 = "old6", ident.2 = "oldex6",  assay ="RNA",verbose = T,test.use = "MAST",logfc.threshold =0)
OPC.ove.diff <- FindMarkers(object = OPC.data, ident.1 = "old6", ident.2 = "oldex6",  assay ="RNA",verbose = T,test.use = "MAST",logfc.threshold =0)
CPC.ove.diff <- FindMarkers(object = CPC.data, ident.1 = "old6", ident.2 = "oldex6",  assay ="RNA",verbose = T,test.use = "MAST",logfc.threshold =0)
OLG.ove.diff <- FindMarkers(object = OLG.data, ident.1 = "old6", ident.2 = "oldex6",  assay ="RNA",verbose = T,test.use = "MAST",logfc.threshold =0)

ovy6.deg <- list(EC.ove=EC.ove.diff, MG.ove=MG.ove.diff, SMC.ove=SMC.ove.diff, PC.ove=PC.ove.diff, OPC.ove=OPC.ove.diff,
                 CPC.ove=CPC.ove.diff, OLG.ove=OLG.ove.diff)

saveRDS(ovy6.deg, file = "./ovy6Degs_updated.rds")
