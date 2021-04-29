# Reanalysis the EC batch5 data for Glp1r project
library(Seurat)
library(tidyverse)
library(devtools)
library(ggpmisc)
library(ggpubr)
library(ggrepel)
source_gist("524eade46135f6348140") # Load the add formular to plot function
#______________________________________________________________________________________________________________________________#
# 1. Data import####

old5.data <- Read10X(data.dir ="~/bioinfo/Richard/NVUB5/Raw/old5/")
young5.data <- Read10X(data.dir ="~/bioinfo/Richard/NVUB5/Raw/young5/")
oldEx5.data <- Read10X(data.dir ="~/bioinfo/Richard/NVUB5/Raw/oldEx5/")

old5 <- CreateSeuratObject(counts  = old5.data, project = "old5",min.cells = 5)
young5 <- CreateSeuratObject(counts  = young5.data, project = "young5",min.cells = 5)
oldEx5 <- CreateSeuratObject(counts = oldEx5.data, project = 'oldEx5', min.cells = 5)

old5 <- RenameCells(old5, add.cell.id = "B5O") 
young5 <- RenameCells(young5, add.cell.id = "B5Y")
oldEx5 <- RenameCells(oldEx5, add.cell.id = "B5E") 


ovy5 <- merge(x = old5,y = list(oldEx5,young5) ,project = "glp1r")

ovy5 <- PercentageFeatureSet(ovy5, pattern = "^mt-", col.name = "percent.mt")
ovy5 <- subset(x = ovy5, subset = quantile(ovy5$nCount_RNA,0.95,na.rm = T)  > nCount_RNA & nCount_RNA > quantile(ovy5$nCount_RNA,0.05,na.rm = T) & percent.mt < 20 )

ovy5 <- SCTransform(ovy5,verbose = T)

ovy5 <- RunPCA(ovy5)
ovy5 <- RunUMAP(ovy5,dims = 1:30,umap.method = "umap-learn")
ovy5 <- FindNeighbors(ovy5,dims = 1:30)
#ovy5 <- FindClusters(ovy5,resolution = 1.5)
ovy5 <- FindClusters(ovy5,resolution = 1.4)
#ovy5 <- FindClusters(ovy5,resolution = 1.2)
#______________________________________________________________________________________________________________________________#
# 3. Celltype identify ####
Idents(ovy5) <- "SCT_snn_res.1.4"
Idents(ovy5, cells = WhichCells(object = ovy5, idents = c("21")))<-"PC"
Idents(ovy5, cells = WhichCells(object = ovy5, idents = c("15","27")))<-"SMC"
Idents(ovy5, cells = WhichCells(object = ovy5, idents = c("1","2","4","14","23")))<-"EC"
Idents(ovy5, cells = WhichCells(object = ovy5, idents = c("6","11","7","8","12","19","29")))<-"MG"
Idents(ovy5, cells = WhichCells(object = ovy5, idents = c("3","0","5","38")))<-"AC"
Idents(ovy5, cells = WhichCells(object = ovy5, idents = c("9","31")))<-"OPC"
Idents(ovy5, cells = WhichCells(object = ovy5, idents = c("40")))<-"OLG"
Idents(ovy5, cells = WhichCells(object = ovy5, idents = c("24")))<-"imNeur"
Idents(ovy5, cells = WhichCells(object = ovy5, idents = c("46","48")))<-"mNeur"
Idents(ovy5, cells = WhichCells(object = ovy5, idents = c("33")))<-"NRP"
Idents(ovy5, cells = WhichCells(object = ovy5, idents = c("36")))<-"EPC"
Idents(ovy5, cells = WhichCells(object = ovy5, idents = c("16","20","30")))<-"CPC"
Idents(ovy5, cells = WhichCells(object = ovy5, idents = c("32")))<-"Hb_EC"
Idents(ovy5, cells = WhichCells(object = ovy5, idents = c("22")))<-"MAC"
Idents(ovy5, cells = WhichCells(object = ovy5, idents = c("52")))<-"TNC"
Idents(ovy5, cells = WhichCells(object = ovy5, idents = c("45")))<-"MNC"
ovy5$Celltype <- Idents(ovy5)
ovy5 <- subset(ovy5,idents = c("PC","SMC","EC","MG","AC","OPC","NRP","OLG","EPC","TNC","imNeur","mNeur","CPC","Hb_EC","MAC","MNC"))

ovy5.less <- subset(ovy5, cells= sample(colnames(ovy5),4000)) #Downsample to 6000 cells
#______________________________________________________________________________________________________________________________#
# 4. DEG calculation ####
Idents(ovy5) <- "Celltype"
PC.data <- subset(ovy5,idents = c("PC"))
SMC.data <- subset(ovy5,idents = c("SMC"))
EC.data <- subset(ovy5,idents = c("EC"))
MG.data <- subset(ovy5,idents = c("MG"))
AC.data <- subset(ovy5,idents = c("AC"))
OPC.data <- subset(ovy5,idents = c("OPC"))
OLG.data <- subset(ovy5,idents = c("OLG"))
imNeur.data <- subset(ovy5,idents = c("imNeur"))
mNeur.data <- subset(ovy5,idents = c("mNeur"))
CPC.data <- subset(ovy5,idents = c("CPC"))
Hb.data <- subset(ovy5,idents = c("Hb_EC"))
MAC.data <- subset(ovy5,idents = c("MAC"))
MNC.data <- subset(ovy5,idents = c("MNC"))

celltype.list <- list(PC.data,SMC.data,EC.data,MG.data,AC.data,OPC.data,OLG.data,
                      imNeur.data,mNeur.data,CPC.data,Hb.data,MAC.data,MNC.data)

ovy.deg <- list() # old vs young
for (i in 1:length(celltype.list)) {
  Idents(celltype.list[[i]]) <- "orig.ident"
  ovy.deg[[i]] <- FindMarkers(celltype.list[[i]],ident.1 = "old5",ident.2 = "young5",verbose = T,test.use = "MAST",logfc.threshold =0)
}
names(ovy.deg) <- c("PC","SMC","EC","MG","AC","OPC","OLG","imNeur","mNeur","CPC","Hb_EC","MAC","MNC")


evo.deg <- list() # old with exenatide vs old
for (i in 1:length(celltype.list)) {
  Idents(celltype.list[[i]]) <- "orig.ident"
  evo.deg[[i]] <- FindMarkers(celltype.list[[i]],ident.1 = "oldEx5",ident.2 = "old5",verbose = T,test.use = "MAST",logfc.threshold =0)
}

names(evo.deg) <- c("PC","SMC","EC","MG","AC","OPC","OLG","imNeur","mNeur","CPC","Hb_EC","MAC","MNC")


# Export the deg list

deg.expt <- function(obj,path){
  for (i in 1:length(obj)) {
    write_csv(as_tibble(obj[[i]],rownames="Gene"), paste0(path,names(obj)[i],"_deg_all.csv"))
  }
  
}


deg.path <- "~/bioinfo/Richard/NVUB5/Text/DEG_list/OVY/"
deg.expt(ovy.deg,deg.path)

deg.path <- "~/bioinfo/Richard/NVUB5/Text/DEG_list/EVO/"
deg.expt( .deg,deg.path)

# 4.1 save some DEGs as Gene Analytic input

siggene.path <- "~/bioinfo/Richard/NVUB5/Text/DEG_sig(GeneAnalytics)/Top500_new/"
x <- evo.deg
y <- ovy.deg
for (i in 1:length(x)) {
  df1 <- x[[i]]
  df2 <- y[[i]]
  its.genes <- intersect(rownames(df1),rownames(df2))
  sig1=(df1[its.genes, "p_val_adj"]<0.05)*1
  sig2=(df2[its.genes, "p_val_adj"]<0.05)*2
  sig=sig1+sig2
  df <- tibble(evo=df1[its.genes,"avg_logFC"],ovy=df2[its.genes,"avg_logFC"],
               sig=sig ,row.names = its.genes)
  df <- df %>% filter(sig>0) # Filter the unsignificant genes
  new.levels <- c("evo_sig","ovy_sig","both_sig")
  df$sig <- factor(new.levels[df$sig],levels = c("evo_sig","ovy_sig","both_sig"))
  #df <- df %>% filter(sig!="evo_sig")
  df$dir <- ifelse(sign(df$evo*df$ovy)==1,"pos","neg")
  df <- df %>% filter(dir=="neg") #%>% filter((abs(ovy)>quantile(abs(df$ovy),0.6))&(abs(evo)>quantile(abs(df$evo),0.6)))
  ovy300 <- top_n(df,500,abs(ovy))
  evo300 <- top_n(df,500,abs(evo))
  df <- intersect(ovy300$row.names,evo300$row.names)
  write_csv(enframe(df),paste0(siggene.path,names(x)[i],"_fctop500.csv"))
} # New filter method by top rank of FC


#______________________________________________________________________________________________________________________________#
# 5. AC subpopulation analysis ####
DefaultAssay(AC.data) <- "RNA"
AC.data[["SCT"]] <- NULL # Redo the normalization and scale

AC.data <- SCTransform(AC.data)

AC.data <- RunPCA(AC.data)
AC.data <- RunUMAP(AC.data,dims = 1:50,umap.method = "umap-learn")
AC.data <- FindNeighbors(AC.data,dims = 1:50)
AC.data <- FindClusters(AC.data,resolution = 1)

Idents(AC.data) <- "SCT_snn_res.1"
Idents(AC.data, cells = WhichCells(object = AC.data, idents = c("5","15","17","7","11","12","13"),invert = T))<-"Gfap_low"
Idents(AC.data, cells = WhichCells(object = AC.data, idents = c("5","15","17")))<-"DAA"
Idents(AC.data, cells = WhichCells(object = AC.data, idents = c("7","11","12","13")))<-"Gfap_high"

AC.data$DAA <- Idents(AC.data)

DAA.num <- as_tibble(table(AC.data[[c("DAA","orig.ident")]]))
DAA.num %>% group_by(orig.ident) %>% mutate(sum=sum(n)) %>% mutate(prop=n/sum) # Calculate the proportion of DAA in individual groups

# 4.1 Diffusion map analysis ####

# 4.2 Seperate the AC into regional subpopulation and find DAA
Idents(AC.data) <- "SCT_snn_res.1"

Idents(AC.data, cells = WhichCells(object = AC.data, idents = c("1")))<-"ACNT1"
Idents(AC.data, cells = WhichCells(object = AC.data, idents = c("3")))<-"ACNT2"

Idents(AC.data, cells = WhichCells(object = AC.data, idents = c("11","12","15")))<-"Radial glia"

Idents(AC.data, cells = WhichCells(object = AC.data, idents = c("17"))) <- "ACMB"

Idents(AC.data, cells = WhichCells(object = AC.data, idents = c("2","5")))<-"ACTE1"
Idents(AC.data, cells = WhichCells(object = AC.data, idents = c("0","6","4","13","7")))<-"ACTE2"


AC.data$Lin <- Idents(AC.data)

AC.data <- subset(AC.data,idents = c("ACNT1","ACNT2","ACTE1","ACTE2"))

Rad.bc <- enframe(AC.data$Lin,name = "bc", value = "Lin") %>% filter( Lin == "Radial glia") %>% pull(bc)

Clean.AC <- enframe(AC.data$Lin,name = "bc", value = "Lin") %>% filter( Lin != "Radial glia") %>% pull(bc)

AC.data <- subset(ovy5,idents = c("AC"))

AC.bc <- colnames(AC.data)[!(colnames(AC.data) %in% Clean.AC)]

AC.data.clean <- subset(AC.data, cells = AC.bc,invert=T) # Remove the radial cells and contaminations
# 4.3 Remove the radial glia AC and redo some analysis

# 4.3.1 The data for F1.

# 4.3.2 Redo the DEG analysis

AC.data <- subset(AC.data, cells = AC.bc,invert=T) # Remove the radial cells and contaminations

celltype.list[[5]] <- AC.data
Idents(celltype.list[[5]]) <- "orig.ident"
ovy.deg$AC <- FindMarkers(celltype.list[[5]],ident.1 = "old5",ident.2 = "young5",verbose = T,test.use = "MAST",logfc.threshold =0)
evo.deg$AC <- FindMarkers(celltype.list[[5]],ident.1 = "oldEx5",ident.2 = "old5",verbose = T,test.use = "MAST",logfc.threshold =0)

Idents(AC.data.clean) <- "orig.ident"
AC.ovy.clean <- FindMarkers(AC.data.clean,ident.1 = "old5",ident.2 = "young5",verbose = T,test.use = "MAST",logfc.threshold =0)
AC.evo.clean <- FindMarkers(AC.data.clean,ident.1 = "oldEx5",ident.2 = "old5",verbose = T,test.use = "MAST",logfc.threshold =0)

# 4.3.3 AC subpopulation DEGs

AC.sub.data <- list()
AC.sub.ident <- c("ACNT1","ACNT2","ACTE1","ACTE2")

for (i in 1:length(AC.sub.ident)) {
  Idents(AC.data) <- "Lin"
  AC.sub.data[[i]] <- subset(AC.data,idents = AC.sub.ident[i])
}
names(AC.sub.data) <- c("ACNT1","ACNT2","ACTE1","ACTE2")


ovy.AC.sub.deg <- list()
for (i in 1:length(AC.sub.data)) {
  Idents(AC.sub.data[[i]]) <- "orig.ident"
  ovy.AC.sub.deg[[i]] <- FindMarkers(AC.sub.data[[i]],ident.1 = "old5",ident.2 = "young5",verbose = T,test.use = "MAST",logfc.threshold =0)
}
names(ovy.AC.sub.deg) <- c("ACNT1","ACNT2","ACTE1","ACTE2")

evo.AC.sub.deg <- list()
for (i in 1:length(AC.sub.data)) {
  Idents(AC.sub.data[[i]]) <- "orig.ident"
  evo.AC.sub.deg[[i]] <- FindMarkers(AC.sub.data[[i]],ident.1 = "oldEx5",ident.2 = "old5",verbose = T,test.use = "MAST",logfc.threshold =0)
}
names(evo.AC.sub.deg) <- c("ACNT1","ACNT2","ACTE1","ACTE2")



# 5. Microglia analysis (from https://www.nature.com/articles/s41586-020-2496-1)

admicro.gene <- read_csv("./Outsource/AD_signatures/Amit_expression_mic3tomic1.csv")

adgene3 <- admicro.gene %>% filter(admicro.gene$`up/down`>0) %>% top_n(100,`-log10(p-value) (Microglia1 vs. Microglia3)`) %>%  pull(X1)
adgene3 <- c("Cst7","Lpl",adgene3)
adgene3 <- adgene3[adgene3 %in% rownames(MG.data)]
#adgene1 <- admicro.gene %>% filter(admicro.gene$`up/down`<0) %>% pull(X1)

MG.data <- subset(ovy5,idents = c("MG"))
DefaultAssay(MG.data) <- "RNA"
MG.data[["SCT"]] <- NULL # Redo the normalization and scale

MG.data.SCT <- SCTransform(MG.data)

MG.data.SCT <- RunPCA(MG.data.SCT)
MG.data.SCT <- RunUMAP(MG.data.SCT,dims = 1:50,umap.method = "umap-learn")
# MG.data.SCT <- FindNeighbors(MG.data.SCT,dims = 1:50)
# MG.data.SCT <- FindClusters(MG.data.SCT,resolution = 1)

#MG.data <- CellCycleScoring(MG.data, g2m.features = adgene3, s.features = adgene1)

MG.data.SCT <- AddModuleScore(
  object = MG.data.SCT,
  features = list(c(adgene3)),
  ctrl = 100,
  name = 'AD_top50'
)

MG.data$AD_top50 <- MG.data.SCT$AD_top501
# 6. SMC subcluster analysis (by Cellranger)

SMC.data <- subset(ovy5,idents = c("SMC"))

DefaultAssay(SMC.data) <- "RNA"
SMC.data[["SCT"]] <- NULL # Redo the normalization and scale

SMC.data.SCT <- SCTransform(SMC.data)

SMC.data <- as.SingleCellExperiment(SMC.data)
saveRDS(SMC.data,file = "./Inter_data/Cellranger/SMC_SCT.rds")


SMC.celltype <- readRDS("./Inter_data/Cellranger/SMCfit_res.rds")

SMC.data$subtype <- SMC.celltype
PC.data <- subset(ovy5,idents = c("PC"))
PC.bc <- rep("PC", dim(PC.data)[1])
names(PC.bc) <- rownames(PC.data[["Celltype"]])

SMCPC.data <- subset(ovy5,idents = c("PC","SMC"))

SMCPC.data$Subtype <- c(SMC.celltype,PC.bc )

# 6.1 SMC subtypes markers

DefaultAssay(SMCPC.data) <- "RNA"
SMCPC.data[["SCT"]] <- NULL
SMCPC.data <- SCTransform(SMCPC.data)
Idents(SMCPC.data) <- "Subtype"
SMCPC.mrks.raw <- FindAllMarkers(SMCPC.data, test.use = "MAST",logfc.threshold = 0.2)

# 6.2 SMC subtypes DEGs

Idents(SMCPC.data) <- "Subtype"

SMCPC.idents <- c("aSMC","aaSMC","vSMC","PC")
SMCPC.list <- list()
for (i in 1:4) {
  idt <- SMCPC.idents[i]
  SMCPC.list[[i]] <- subset(SMCPC.data, idents = idt)
}
names(SMCPC.list) <- c("aSMC","aaSMC","vSMC","PC")


SMCPC.ovy.deg <- list()
for(i in 1:4){
  Idents(SMCPC.list[[i]]) <- "orig.ident"
  SMCPC.ovy.deg[[i]] <- FindMarkers(SMCPC.list[[i]],ident.1 = "old5",ident.2 = "young5",verbose = T,test.use = "MAST",logfc.threshold =0)
}
names(SMCPC.ovy.deg) <-  c("aSMC","aaSMC","vSMC","PC")

SMCPC.evo.deg <- list()
for(i in 1:4){
  Idents(SMCPC.list[[i]]) <- "orig.ident"
  SMCPC.evo.deg[[i]] <- FindMarkers(SMCPC.list[[i]],ident.1 = "oldEx5",ident.2 = "old5",verbose = T,test.use = "MAST",logfc.threshold =0)
}
names(SMCPC.ovy.deg) <-  c("aSMC","aaSMC","vSMC","PC")

#Extra analysis

# The aged vs young data

celltype.subset <- list(AC.data,OPC.data,MG.data,MAC.data,OLG.data,SMC.data,PC.data,EC.data)

evy.deg <- list() # old with exenatide vs young
for (i in 1:length(celltype.subset)) {
  Idents(celltype.subset[[i]]) <- "orig.ident"
  evy.deg[[i]] <- FindMarkers(celltype.subset[[i]],ident.1 = "oldEx5",ident.2 = "young5",verbose = T,test.use = "MAST",logfc.threshold =0)
}

names(evy.deg) <- c("AC","OPC","MG","MAC","OLG","SMC","PC","EC")

b5.bar.ct <- names(evy.deg)
for(i in 1:length(b5.bar.ct)){
  if(i == 1){
    ct <- b5.bar.ct[i]
    ovy.deg.nm <- dim(ovy.deg[[ct]][(ovy.deg[[ct]]$p_val_adj<0.05),])[1]
    evy.deg.nm <- dim(evy.deg[[ct]][(evy.deg[[ct]]$p_val_adj<0.05),])[1]
    deg.num <- tibble(ovy.deg=ovy.deg.nm, evy.deg=evy.deg.nm)  
    next()
  }
  ct <- b5.bar.ct[i]
  ovy.deg.nm <- dim(ovy.deg[[ct]][(ovy.deg[[ct]]$p_val_adj<0.05),])[1]
  evy.deg.nm <- dim(evy.deg[[ct]][(evy.deg[[ct]]$p_val_adj<0.05),])[1]
  deg.num <- bind_rows(deg.num, tibble(ovy.deg=ovy.deg.nm, evy.deg=evy.deg.nm))
}

deg.num <- mutate(deg.num, celltype=b5.bar.ct, .before = 1)

deg.num$evy.deg <- deg.num$evy.deg*-1

deg.num$celltype <- factor(deg.num$celltype, levels = rev(deg.num$celltype))
