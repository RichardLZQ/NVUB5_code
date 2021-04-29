# The code generating NVUB5 project figures

# Fig 1a. The cell type umap
Idents(ovy5) <- "Celltype"
ovy5.less <- subset(ovy5, idents = c("AC","OPC","MG","MAC","OLG","PC","SMC","EC")) 
ovy5.less <- subset(ovy5, cells= sample(colnames(ovy5.less),6000)) #Downsample to 6000 cells
ovy5.less$Celltype <- droplevels(ovy5.less$Celltype)
ovy5.less$Celltype <- factor(ovy5.less$Celltype, levels = c("EC","OLG","OPC","AC","MG","MAC","SMC","PC"))
DimPlot(ovy5.less,group.by = "Celltype",label=T)

# Fig 1b. UMAP seperated by sample ids
Idents(ovy5.less) <- "orig.ident"
DimPlot(ovy5.less, cells.highlight = WhichCells(ovy5.less, idents = "old5"),pt.size = 0.2,sizes.highlight = 0.3)
DimPlot(ovy5.less, cells.highlight = WhichCells(ovy5.less, idents = "oldEx5"),pt.size = 0.2,sizes.highlight = 0.3)
DimPlot(ovy5.less, cells.highlight = WhichCells(ovy5.less, idents = "young5"),pt.size = 0.2,sizes.highlight = 0.3)

# Fig 1b. The cell number
Orig.num <- enframe(ovy5$orig.ident,"cellbc","orig")
celltype.num <- enframe(ovy5$Celltype,"cellbc","celltype")

# Fig 1c. The reversal effect coordinate plot
degdot.final <- function(x,y,ct.list){ # Visualized the intersec items between x and y by LogFC in a coordinate colored by adjPvalue
  for (i in 1:length(ct.list)) {
    ct <- ct.list[i]
    df1 <- x[[ct]]
    df2 <- y[[ct]]
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
    formula <- y ~ x
    print(nrow(df))
    print(ggplot(df, aes(x=ovy,y=evo))+ geom_point(size=2,alpha = 0.7,colour="#325d81")+
            stat_smooth_func(geom="text",method="lm",hjust=0,vjust=1,parse=T,fullrange = T)+
            geom_smooth(method = "lm", se = FALSE,colour="#b3b3b3",alpha=0.8)+
            ggtitle(paste0(ct.list[i],"_Celltype"))+theme_bw()+
            stat_fit_glance(method = "lm",
                            method.args = list(formula = formula),
                            label.x = "right",
                            label.y = "bottom",
                            aes(label = paste("italic(P)*\"-value = \"*",
                                              signif(..p.value.., digits = 4), sep = "")),
                            parse = TRUE)+
            scale_x_continuous(breaks=seq(-1,1,0.25),limits = c(-1, 1))+ scale_y_continuous(breaks=seq(-0.75,0.75,0.25),limits = c(-0.8, 0.8))+
            geom_hline(yintercept = 0) +geom_vline(xintercept = 0))
  }
} 

ct.list <- c("AC","OPC","MG","MAC","OLG","SMC","PC","EC")
pdf("./Revision/AC_revised.pdf")
degdot.final(evo.deg,ovy.deg,ct.list = ct.list)
dev.off()


# Fig 1d. The slope and reversal propotion line plot
slope.list <- c()
rev.porp <- c()
cor.dataframe <- data.frame()

x <- evo.deg
y <- ovy.deg
for (i in 1:length(ct.list)) { # Calculate the linear regression and retrive slope 
  ct <- ct.list[i]
  df1 <- x[[ct]]
  df2 <- y[[ct]]
  its.genes <- intersect(rownames(df1),rownames(df2))
  sig1=(df1[its.genes, "p_val_adj"]<0.05)*1
  sig2=(df2[its.genes, "p_val_adj"]<0.05)*2
  sig=sig1+sig2
  df <- tibble(evo=df1[its.genes,"avg_logFC"],ovy=df2[its.genes,"avg_logFC"], 
               sig=sig ,row.names = its.genes)
  df <- df %>% filter(sig>0) # Filter the unsignificant genes
  new.levels <- c("evo_sig","ovy_sig","both_sig")
  df$sig <- factor(new.levels[df$sig],levels = c("evo_sig","ovy_sig","both_sig"))
  df <- df %>% filter(sig!="evo_sig")
  df$dir <- ifelse(sign(df$evo*df$ovy)==1,"pos","neg")
  rev.porp <- c(rev.porp,(sum(df$dir=="neg")/nrow(df)))
  #df <- df %>% filter(ovy>-1 &ovy<1) %>% filter(evo> -0.75 &evo<0.75) #To make sure manually calculated number is same with automatical results
  linearMod <- lm(evo ~ ovy, data=df)
  cor.dataframe <- rbind(cor.dataframe,t(as.data.frame(linearMod$coefficients)))
  slope.list <- c(slope.list,linearMod$coefficients[2])
  ifelse(i==1,
         slop.prop <- tibble(Celltype=ct.list[i],Slope=linearMod$coefficients[2],Propotion=(sum(df$dir=="neg")/nrow(df))),
         slop.prop <- bind_rows(slop.prop, tibble(Celltype=ct.list[i],Slope=linearMod$coefficients[2],Propotion=(sum(df$dir=="neg")/nrow(df)))))
}

rownames(cor.dataframe) <-  c("PC","SMC","EC","MG","AC","OPC","OLG","MAC")
names(slope.list) <-  c("PC","SMC","EC","MG","AC","OPC","OLG","MAC")
slope.list <- enframe(slope.list, name = "Celltype","Slope")
slope.list <- slope.list[c(5,6,4,8,7,2,1,3),]
line.data <- left_join(slope.list,rev.porp,by="Celltype")
line.data$Celltype <- factor(line.data$Celltype,levels = line.data$Celltype)
#line.data$Slope <- line.data$Slope*-1
#line.data$Slope[9] <- 0 # Force negative to 0
#line.data$bk_remove[9] <- 0
line.data$group[4] <- "Glia"
line.data <- line.data[,-5]

line.data2 <- gather(line.data, "Type","Value",-c("Celltype","group")) 
#line.data2[9:,4] <- line.data2[9:24,4]*0.5 # Reshapre the data to fit "overlap" plot
line.data2$group2 <- c(pull(line.data2,group)[1:8],paste0(pull(line.data2,group)[1:8],"2"))
#line.data2$Value[1:8] <- ((line.data2$Value[1:8]-0.1)/0.9)*0.3+0.7
line.data2$Value[1:8] <- line.data2$Value[1:8]/2
line.data2$Value[9:16] <-  line.data2$Value[9:16]-0.5
  
ggplot(data=line.data2, aes(x=Celltype, y=Value, group=group2))+
  geom_line(aes(colour=group,linetype=Type))+
geom_point(aes(colour=group,shape=Type))+scale_y_continuous(
  name = "Pct of DEGs reversed", limits = c(0,0.5),
  sec.axis = sec_axis( trans=~./2, name="Negative slope",breaks = seq(-1,-0.2,0.2))# Add a second axis and specify its features
) 

# Fig 2a. This figure was finished by other software

# Fig 2b. The bubble plot of selected genes from pathway results
# 1. Microglia
mg.marker <- read_csv("./Outsource/Selected_gene/MG_final.csv")
marker.type <- mg.marker$Type
mg.marker <- mg.marker[!duplicated(mg.marker$Gene),]

mgmak.ovy <- as_tibble(ovy.deg$MG,rownames = "Gene")
mgmak.evo <- as_tibble(evo.deg$MG,rownames = "Gene")

mgmak.ovy <- inner_join(mgmak.ovy,mg.marker,by = "Gene")
mgmak.evo <- inner_join(mgmak.evo,mg.marker,by = "Gene")

mgmak.ovy <- mgmak.ovy %>% select(c(1,3,6)) %>% mutate(Type="OVY") %>% filter(p_val_adj!="NA")
mgmak.evo <- mgmak.evo %>% select(c(1,3,6)) %>% mutate(Type="EVO") %>% filter(p_val_adj!="NA")

mgmak.dat <- bind_rows(mgmak.ovy,mgmak.evo)
mgmak.dat$p_val_adj <- ifelse(mgmak.dat$p_val_adj==0,1e-320,mgmak.dat$p_val_adj)
mgmak.dat <- mgmak.dat %>% mutate(`log2pval` = log2(p_val_adj)) %>% mutate(`-log2pval` = ifelse(log2pval<(-200),200,(-1*log2pval)))
mgmak.dat$avg_logFC <- ifelse(mgmak.dat$avg_logFC>0.3,0.3,mgmak.dat$avg_logFC)
mgmak.dat$avg_logFC <- ifelse(mgmak.dat$avg_logFC<(-0.3),(-0.3),mgmak.dat$avg_logFC)

mgmak.dat$Type <- factor(mgmak.dat$Type,levels = c("OVY","EVO"))
mgmak.dat$Gene <- factor(mgmak.dat$Gene,levels = rev(mg.marker$Gene))

ggplot(mgmak.dat, aes(x=Type, y=Gene,size=`-log2pval`)) + geom_point(aes(fill=avg_logFC),colour="black",shape=21,stroke=1)+ # Dot plot
  geom_point(aes(colour=avg_logFC))+scale_size_continuous(range = c(3,10),limits = c(0,200))+ 
  scale_color_gradient2(midpoint = 0, low = "steelblue", mid = "white",
                        high = "#E64B35", space = "Lab",limits=c(-0.31,0.31))


# 2. Astrocyte
ac.marker <- read_csv("./Outsource/Selected_gene/AC_final.csv")
ac.marker <- ac.marker[!duplicated(ac.marker$Gene),]

acmak.ovy <- as_tibble(ovy.deg$AC,rownames = "Gene")
acmak.evo <- as_tibble(evo.deg$AC,rownames = "Gene")

acmak.ovy <- inner_join(acmak.ovy,ac.marker,by = "Gene")
acmak.evo <- inner_join(acmak.evo,ac.marker,by = "Gene")

acmak.ovy <- acmak.ovy %>% select(c(1,3,6)) %>% mutate(Type="OVY") %>% filter(p_val_adj!="NA")
acmak.evo <- acmak.evo %>% select(c(1,3,6)) %>% mutate(Type="EVO") %>% filter(p_val_adj!="NA")

acmak.dat <- bind_rows(acmak.ovy,acmak.evo)
acmak.dat$p_val_adj <- ifelse(acmak.dat$p_val_adj==0,1e-320,acmak.dat$p_val_adj)
acmak.dat <- acmak.dat %>% mutate(`log2pval` = log2(p_val_adj)) %>% mutate(`-log2pval` = ifelse(log2pval<(-200),200,(-1*log2pval)))
acmak.dat$avg_logFC <- ifelse(acmak.dat$avg_logFC<(-0.3),(-0.3),acmak.dat$avg_logFC)
acmak.dat$avg_logFC <- ifelse(acmak.dat$avg_logFC>(0.3),(0.3),acmak.dat$avg_logFC)

acmak.dat$Type <- factor(acmak.dat$Type,levels = c("OVY","EVO"))
acmak.dat$Gene <- factor(acmak.dat$Gene, levels = rev(ac.marker$Gene))

ggplot(acmak.dat, aes(x=Type, y=Gene,size=`-log2pval`)) + geom_point(aes(fill=avg_logFC),colour="black",shape=21,stroke=1)+ # Dot plot
geom_point(aes(colour=avg_logFC))+scale_size_continuous(range = c(3,10),limits = c(0,200))+
scale_color_gradient2(midpoint = 0, low = "steelblue", mid = "white",
                    high = "#E64B35", space = "Lab",limits=c(-0.31,0.31))


# 3. SMC
smc.marker <- read_csv("./Outsource/Selected_gene/SMC_final.csv")
smc.marker <- smc.marker[!duplicated(smc.marker$Gene),]

smcmak.ovy <- as_tibble(ovy.deg$SMC,rownames = "Gene")
smcmak.evo <- as_tibble(evo.deg$SMC,rownames = "Gene")

smcmak.ovy <- inner_join(smcmak.ovy,smc.marker,by = "Gene")
smcmak.evo <- inner_join(smcmak.evo,smc.marker,by = "Gene")

smcmak.ovy <- smcmak.ovy %>% select(c(1,3,6)) %>% mutate(Type="OVY") %>% filter(p_val_adj!="NA")
smcmak.evo <- smcmak.evo %>% select(c(1,3,6)) %>% mutate(Type="EVO") %>% filter(p_val_adj!="NA")

smcmak.dat <- bind_rows(smcmak.ovy,smcmak.evo)
smcmak.dat <- smcmak.dat %>% mutate(`log2pval` = log2(p_val_adj)) %>% mutate(`-log2pval` = ifelse(log2pval<(-200),200,(-1*log2pval)))
smcmak.dat$avg_logFC <- ifelse(smcmak.dat$avg_logFC>0.3,0.3,smcmak.dat$avg_logFC)
smcmak.dat$avg_logFC <- ifelse(smcmak.dat$avg_logFC<(-0.3),(-0.3),smcmak.dat$avg_logFC)

smcmak.dat$Type <- factor(smcmak.dat$Type,levels = c("OVY","EVO"))
smcmak.dat$Gene <- factor(smcmak.dat$Gene, levels = rev(smc.marker$Gene))

ggplot(smcmak.dat, aes(x=Type, y=Gene,size=`-log2pval`)) + geom_point(aes(fill=avg_logFC),colour="black",shape=21,stroke=1)+ # Dot plot
  geom_point(aes(colour=avg_logFC))+scale_size_continuous(range = c(3,10),limits = c(0,200))+
  scale_color_gradient2(midpoint = 0, low = "steelblue", mid = "white",
                        high = "#E64B35", space = "Lab",limits=c(-0.31,0.31))

#Fig 2c. The feature plot of AD signature in MG and separated by sample groups.
MG.data.less <- subset(MG.data, cells = sample(colnames(MG.data),6000))
pdf("./Figures/Figure_materials/MG_ADscore_UMAP.pdf",width = 5.5, height = 5.5)
Idents(MG.data.less) <- "orig.ident"
FeaturePlot(MG.data.less,features = "AD_top50",cells = WhichCells(MG.data.less,idents = "young5")) +expand_limits(colour = seq(0, 1.5, by = 0.1))+
  scale_x_continuous(limits = c(6,24),breaks=seq(6,24,2))+
  scale_y_continuous(limits = c(-2,10),breaks=seq(-2,10,2))+theme(legend.position = "none")+ggtitle("")
FeaturePlot(MG.data.less,features = "AD_top50",cells = WhichCells(MG.data.less,idents = "old5")) +expand_limits(colour = seq(0, 1.5, by = 0.1))+
  scale_x_continuous(limits = c(6,24),breaks=seq(6,24,2))+
  scale_y_continuous(limits = c(-2,10),breaks=seq(-2,10,2))+theme(legend.position = "none")+ggtitle("")
FeaturePlot(MG.data.less,features = "AD_top50",cells = WhichCells(MG.data.less,idents = "oldEx5")) +expand_limits(colour = seq(0, 1.5, by = 0.1))+
  scale_x_continuous(limits = c(6,24),breaks=seq(6,24,2))+
  scale_y_continuous(limits = c(-2,10),breaks=seq(-2,10,2))+theme(legend.position = "none")+ggtitle("")
dev.off()

#Fig 2d. The feature plot of AD signature in MG and separated by sample groups.
pdf("./Figures/Figure_materials/MD_ADScore_with_median_05quatile",width = 5.5, height = 5.5)
VlnPlot(MG.data,features = "AD_top50",group.by = "orig.ident",pt.size = 0,y.max=1)+
  geom_boxplot(width=0.2,outlier.size=0)
dev.off()

#Fig S1. The the expression of marker genes in UMAP and violin plot
pdf("./Figures/Figure_materials/Cluster_markers_new.pdf") # Post-processed by Adobe
VlnPlot(object = ovy5.clean, features = c("Kcnj8"),pt.size = 0,ncol = 1,y.max = 8,idents = c("AC","OPC","MG","MAC","OLG","SMC","PC","EC"))
VlnPlot(object = ovy5.clean, features = c("Acta2"),pt.size = 0,ncol = 1,y.max = 8,idents = c("AC","OPC","MG","MAC","OLG","SMC","PC","EC"))
VlnPlot(object = ovy5.clean, features = c("Cldn5"),pt.size = 0,ncol = 1,y.max = 8,idents = c("AC","OPC","MG","MAC","OLG","SMC","PC","EC"))
VlnPlot(object = ovy5.clean, features = c("Ctss"),pt.size = 0,ncol = 1,y.max = 8,idents = c("AC","OPC","MG","MAC","OLG","SMC","PC","EC"))
VlnPlot(object = ovy5.clean, features = c("Ntsr2"),pt.size = 0,ncol = 1,y.max = 8,idents = c("AC","OPC","MG","MAC","OLG","SMC","PC","EC"))
VlnPlot(object = ovy5.clean, features = c("Pdgfra"),pt.size = 0,ncol = 1,y.max = 8,idents = c("AC","OPC","MG","MAC","OLG","SMC","PC","EC"))
VlnPlot(object = ovy5.clean, features = c("Cldn11"),pt.size = 0,ncol = 1,y.max = 8,idents = c("AC","OPC","MG","MAC","OLG","SMC","PC","EC"))
VlnPlot(object = ovy5.clean, features = c("Pf4"),pt.size = 0,ncol = 1,y.max = 8,idents = c("AC","OPC","MG","MAC","OLG","SMC","PC","EC"))
dev.off()

ovy5.less <- subset(ovy5.less, cells = AC.bc,invert=T) # Remove the radial cells and contaminations
pdf("./Figures/Figure_materials/Cluster_markers_UMAP.pdf")
FeaturePlot(ovy5.less,features = "Kcnj8",label = T,order=T,cells = WhichCells(ovy5.less, idents = c("AC","OPC","MG","MAC","OLG","SMC","PC","EC")))+expand_limits(colour = seq(0, 5, by = 1))
FeaturePlot(ovy5.less,features = "Acta2",label = T,order=T,cells = WhichCells(ovy5.less, idents = c("AC","OPC","MG","MAC","OLG","SMC","PC","EC")),max.cutoff = 5)+expand_limits(colour = seq(0, 5, by = 1))
FeaturePlot(ovy5.less,features = "Cldn5",label = T,order=T,cells = WhichCells(ovy5.less, idents = c("AC","OPC","MG","MAC","OLG","SMC","PC","EC")))+expand_limits(colour = seq(0, 5, by = 1))
FeaturePlot(ovy5.less,features = "Ctss",label = T,order=T,cells = WhichCells(ovy5.less, idents = c("AC","OPC","MG","MAC","OLG","SMC","PC","EC")))+expand_limits(colour = seq(0, 5, by = 1))
FeaturePlot(ovy5.less,features = "Ntsr2",label = T,order=T,cells = WhichCells(ovy5.less, idents = c("AC","OPC","MG","MAC","OLG","SMC","PC","EC")))+expand_limits(colour = seq(0, 5, by = 1))
FeaturePlot(ovy5.less,features = "Pdgfra",label = T,order=T,cells = WhichCells(ovy5.less, idents = c("AC","OPC","MG","MAC","OLG","SMC","PC","EC")))+expand_limits(colour = seq(0, 5, by = 1))
FeaturePlot(ovy5.less,features = "Cldn11",label = T,order=T,cells = WhichCells(ovy5.less, idents = c("AC","OPC","MG","MAC","OLG","SMC","PC","EC")))+expand_limits(colour = seq(0, 5, by = 1))
FeaturePlot(ovy5.less,features = "Pf4",label = T,order=T,cells = WhichCells(ovy5.less, idents = c("AC","OPC","MG","MAC","OLG","SMC","PC","EC")))+expand_limits(colour = seq(0, 5, by = 1))
dev.off()

#Fig S2. The DEG number of selected cell types in two groups 

pdf("./Revision/Rev1/DEG_barplot_withoutFC.pdf")
ggplot(data = deg.num)+geom_bar(aes(x=celltype,y=ovy.deg),position = "stack", stat="identity",fill="#E64B35",width =  0.75)+
  geom_bar(aes(x=celltype,y=evy.deg),position = "stack", stat="identity",fill="steelblue",width = 0.75)+
  coord_flip()
dev.off()

#Fig S3. Regional AC subtypes
pdf("./Figures/Figure_materials/Subpop_markers_UMAP_extra_orig.pdf")
DimPlot(ovy5.AC,label = F,group.by = "Lin",)+theme(legend.position = "none")+ ggtitle("")+
  scale_x_continuous(limits = c(-2,8),breaks=seq(-2,8,2))+
  scale_y_continuous(limits = c(0,12),breaks=seq(3,9,3))
FeaturePlot(ovy5.AC,features = "Gfap",coord.fixed = F)+expand_limits(colour = seq(0, 5, by = 1))+theme(legend.position = "none")+ggtitle("") +
  scale_x_continuous(limits = c(-2,8),breaks=seq(-2,8,2))+
  scale_y_continuous(limits = c(0,12),breaks=seq(3,9,3))
FeaturePlot(ovy5.AC,features = "Hopx")+expand_limits(colour = seq(0, 5, by = 1))+theme(legend.position = "none")+ggtitle("") +
  scale_x_continuous(limits = c(-2,8),breaks=seq(-2,8,2))+
  scale_y_continuous(limits = c(0,12),breaks=seq(3,9,3))
FeaturePlot(ovy5.AC,features = "Atp1b1")+expand_limits(colour = seq(0, 5, by = 1))+theme(legend.position = "none")+ggtitle("") +
  scale_x_continuous(limits = c(-2,8),breaks=seq(-2,8,2))+
  scale_y_continuous(limits = c(0,12),breaks=seq(3,9,3))
FeaturePlot(ovy5.AC,features = "Ppp1r3g")+expand_limits(colour = seq(0, 5, by = 1))+theme(legend.position = "none")+ggtitle("") +
  scale_x_continuous(limits = c(-2,8),breaks=seq(-2,8,2))+
  scale_y_continuous(limits = c(0,12),breaks=seq(3,9,3))
FeaturePlot(ovy5.AC,features = "Agt")+expand_limits(colour = seq(0, 5, by = 1))+theme(legend.position = "none")+ggtitle("") +
  scale_x_continuous(limits = c(-2,8),breaks=seq(-2,8,2))+
  scale_y_continuous(limits = c(0,12),breaks=seq(3,9,3))
FeaturePlot(ovy5.AC,features = "Slc6a9")+expand_limits(colour = seq(0, 5, by = 1))+theme(legend.position = "none")+ggtitle("") +
  scale_x_continuous(limits = c(-2,8),breaks=seq(-2,8,2))+
  scale_y_continuous(limits = c(0,12),breaks=seq(3,9,3))
FeaturePlot(ovy5.AC,features = "Aldh1a1")+expand_limits(colour = seq(0, 5, by = 1))+theme(legend.position = "none")+ggtitle("") +
  scale_x_continuous(limits = c(-2,8),breaks=seq(-2,8,2))+
  scale_y_continuous(limits = c(0,12),breaks=seq(3,9,3))
FeaturePlot(ovy5.AC,features = "Itih3")+expand_limits(colour = seq(0, 5, by = 1))+theme(legend.position = "none")+ggtitle("") +
  scale_x_continuous(limits = c(-2,8),breaks=seq(-2,8,2))+
  scale_y_continuous(limits = c(0,12),breaks=seq(3,9,3))
dev.off()

Idents(AC.data) <- "Lin"
ht.genes <- FindAllMarkers(AC.data, test.use = "MAST",logfc.threshold = 0.2)
dim(ht.genes)

ht.genes.raw <- FindAllMarkers(AC.data, test.use = "MAST",logfc.threshold = 0.2)
ht.genes <- ht.genes.raw %>% as_tibble(rownames = "Gene") %>% group_by(cluster) %>% top_n(15,avg_logFC)
ht.genes$cluster <- factor(ht.genes$cluster, levels = c("ACTE1","ACTE2","ACNT1","ACNT2"))
ht.genes <- ht.genes %>% dplyr::arrange(cluster,.by_group= T)
AC.data.less <- subset(AC.data, cells= sample(colnames(AC.data),4000))
DoHeatmap(AC.data.less,features = ht.genes$gene,disp.min = -3, disp.max = 3)

#Fig S4. Consistency of differential expressions in AC subtype

cor.ident <- t(combn(1:4,2)) # Generate the permutation combine numbers

x <- ovy.AC.sub.deg
pdf("./Figures/Figure_materials/subpop_coordi_ovy_updated.pdf.pdf")
for (i in 1:6) {
  ident.par <- cor.ident[i,]
  df1 <- x[[ident.par[1]]]
  df2 <- x[[ident.par[2]]]
  its.genes <- intersect(rownames(df1),rownames(df2))
  sig1=(df1[its.genes, "p_val_adj"]<0.05)*1
  sig2=(df2[its.genes, "p_val_adj"]<0.05)*2
  sig=sig1+sig2
  df <- tibble(evo=df1[its.genes,"avg_logFC"],ovy=df2[its.genes,"avg_logFC"],
               sig=sig ,row.names = its.genes)
  df <- df %>% filter(sig>0) # Filter the unsignificant genes
  formula <- y ~ x
  
  print(ggplot(df, aes(x=evo,y=ovy))+ geom_point(size=2,alpha = 0.7,colour="#325d81")+
          xlab(AC.sub.ident[ident.par[1]]) + ylab(AC.sub.ident[ident.par[2]]) +
          stat_smooth_func(geom="text",method="lm",hjust=0,vjust=1,parse=T,fullrange = T)+
          geom_smooth(method = "lm", se = FALSE,colour="#b3b3b3",alpha=0.8)+
          ggtitle(paste0(AC.sub.ident[ident.par[1]],"_vs_",AC.sub.ident[ident.par[2]]))+theme_bw()+
          stat_fit_glance(method = "lm",
                          method.args = list(formula = formula),
                          label.x = "right",
                          label.y = "bottom",
                          aes(label = paste("italic(P)*\"-value = \"*",
                                            signif(..p.value.., digits = 4), sep = "")),
                          parse = TRUE)+
          scale_x_continuous(breaks=seq(-1,1,0.25),limits = c(-1, 1))+ scale_y_continuous(breaks=seq(-1,1,0.25),limits = c(-1, 1))+
          geom_hline(yintercept = 0) +geom_vline(xintercept = 0))
}
dev.off()

x <- evo.AC.sub.deg
pdf("./Figures/Figure_materials/subpop_coordi_evo_updated.pdf.pdf")
for (i in 1:6) {
  ident.par <- cor.ident[i,]
  df1 <- x[[ident.par[1]]]
  df2 <- x[[ident.par[2]]]
  its.genes <- intersect(rownames(df1),rownames(df2))
  sig1=(df1[its.genes, "p_val_adj"]<0.05)*1
  sig2=(df2[its.genes, "p_val_adj"]<0.05)*2
  sig=sig1+sig2
  df <- tibble(evo=df1[its.genes,"avg_logFC"],ovy=df2[its.genes,"avg_logFC"],
               sig=sig ,row.names = its.genes)
  df <- df %>% filter(sig>0) # Filter the unsignificant genes
  formula <- y ~ x
  
  print(ggplot(df, aes(x=evo,y=ovy))+ geom_point(size=2,alpha = 0.7,colour="#325d81")+
          xlab(AC.sub.ident[ident.par[1]]) + ylab(AC.sub.ident[ident.par[2]]) +
          stat_smooth_func(geom="text",method="lm",hjust=0,vjust=1,parse=T,fullrange = T)+
          geom_smooth(method = "lm", se = FALSE,colour="#b3b3b3",alpha=0.8)+
          ggtitle(paste0(AC.sub.ident[ident.par[1]],"_vs_",AC.sub.ident[ident.par[2]]))+theme_bw()+
          stat_fit_glance(method = "lm",
                          method.args = list(formula = formula),
                          label.x = "right",
                          label.y = "bottom",
                          aes(label = paste("italic(P)*\"-value = \"*",
                                            signif(..p.value.., digits = 4), sep = "")),
                          parse = TRUE)+
          scale_x_continuous(breaks=seq(-1,1,0.25),limits = c(-1, 1))+ scale_y_continuous(breaks=seq(-1,1,0.25),limits = c(-1, 1))+
          geom_hline(yintercept = 0) +geom_vline(xintercept = 0))
}
dev.off()

x <- evo.AC.sub.deg
y <- ovy.AC.sub.deg
pdf("./Figures/Figure_materials/subpop_reversal.pdf")
for (i in 1:length(AC.sub.ident)) {
  ct <- AC.sub.ident[i]
  df1 <- x[[ct]]
  df2 <- y[[ct]]
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
  formula <- y ~ x
  # print(nrow(df))
  # return(nrow(df))
  print(ggplot(df, aes(x=ovy,y=evo))+ geom_point(size=2,alpha = 0.7,colour="#325d81")+
          stat_smooth_func(geom="text",method="lm",hjust=0,vjust=1,parse=T,fullrange = T)+
          geom_smooth(method = "lm", se = FALSE,colour="#b3b3b3",alpha=0.8)+
          ggtitle(paste0(ct,"_Subtype"))+theme_bw()+
          stat_fit_glance(method = "lm",
                          method.args = list(formula = formula),
                          label.x = "right",
                          label.y = "bottom",
                          aes(label = paste("italic(P)*\"-value = \"*",
                                            signif(..p.value.., digits = 4), sep = "")),
                          parse = TRUE)+
          scale_x_continuous(breaks=seq(-1,1,0.25),limits = c(-1, 1))+ scale_y_continuous(breaks=seq(-1,1,0.25),limits = c(-1, 1))+
          geom_hline(yintercept = 0) +geom_vline(xintercept = 0))
}
dev.off()

# Fig S5. Mural cell subtypes

pdf("./Figures/Figure_materials/SMCPC_marker_UMAP_updated1.pdf",width = 5.5, height = 5.5)
FeaturePlot(SMCPC.data, features = "Cnn1")+expand_limits(colour = seq(0, 4, by = 0.5))+theme(legend.position = "none")+ggtitle("") +
  scale_x_continuous(limits = c(-9,9),breaks=seq(-9,9,3))+
  scale_y_continuous(limits = c(-2,3),breaks=seq(-2,3,1))
FeaturePlot(SMCPC.data, features = "Acta2")+expand_limits(colour = seq(0, 4, by = 0.5))+theme(legend.position = "none")+ggtitle("") +
  scale_x_continuous(limits = c(-9,9),breaks=seq(-9,9,3))+
  scale_y_continuous(limits = c(-2,3),breaks=seq(-2,3,1))
FeaturePlot(SMCPC.data, features = "Ccnd1")+expand_limits(colour = seq(0, 4, by = 0.5))+theme(legend.position = "none")+ggtitle("") +
  scale_x_continuous(limits = c(-9,9),breaks=seq(-9,9,3))+
  scale_y_continuous(limits = c(-2,3),breaks=seq(-2,3,1))
FeaturePlot(SMCPC.data, features = "Kcnj8")+expand_limits(colour = seq(0, 4, by = 0.5))+theme(legend.position = "none")+ggtitle("") +
  scale_x_continuous(limits = c(-9,9),breaks=seq(-9,9,3))+
  scale_y_continuous(limits = c(-2,3),breaks=seq(-2,3,1))
dev.off()


SMCPC.mrks <- SMCPC.mrks.raw %>% as_tibble(rownames = "Gene") %>% group_by(cluster) %>% top_n(15,avg_logFC)
SMCPC.mrks <- SMCPC.mrks %>% dplyr::arrange(cluster,.by_group= T)
DoHeatmap(SMCPC.data,features = SMCPC.mrks$gene,disp.min = -3, disp.max = 3)

# Fig S6. Mural cell reversal after treatment

x <- SMCPC.evo.deg
y <- SMCPC.ovy.deg
pdf("./Figures/Figure_materials/SMCPC_subtype_reversal.pdf")
for (i in 1:length(SMCPC.idents)) {
  ct <- SMCPC.idents[i]
  df1 <- x[[ct]]
  df2 <- y[[ct]]
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
  formula <- y ~ x
  # print(nrow(df))
  # return(nrow(df))
  print(ggplot(df, aes(x=ovy,y=evo))+ geom_point(size=2,alpha = 0.7,colour="#325d81")+
          stat_smooth_func(geom="text",method="lm",hjust=0,vjust=1,parse=T,fullrange = T)+
          geom_smooth(method = "lm", se = FALSE,colour="#b3b3b3",alpha=0.8)+
          ggtitle(paste0(ct,"_Subtype"))+theme_bw()+
          stat_fit_glance(method = "lm",
                          method.args = list(formula = formula),
                          label.x = "right",
                          label.y = "bottom",
                          aes(label = paste("italic(P)*\"-value = \"*",
                                            signif(..p.value.., digits = 4), sep = "")),
                          parse = TRUE)+
          scale_x_continuous(breaks=seq(-1,1,0.25),limits = c(-1, 1))+ scale_y_continuous(breaks=seq(-1,1,0.25),limits = c(-1, 1))+
          geom_hline(yintercept = 0) +geom_vline(xintercept = 0))
}
dev.off()

cor.ident <- t(combn(1:4,2)) # Generate the permutation combine numbers
x <- SMCPC.ovy.deg
pdf("./Figures/Figure_materials/subpop_coordi_ovy.pdf")
for (i in 1:6) {
  ident.par <- cor.ident[i,]
  df1 <- x[[ident.par[1]]]
  df2 <- x[[ident.par[2]]]
  its.genes <- intersect(rownames(df1),rownames(df2))
  sig1=(df1[its.genes, "p_val_adj"]<0.05)*1
  sig2=(df2[its.genes, "p_val_adj"]<0.05)*2
  sig=sig1+sig2
  df <- tibble(evo=df1[its.genes,"avg_logFC"],ovy=df2[its.genes,"avg_logFC"],
               sig=sig ,row.names = its.genes)
  df <- df %>% filter(sig>0) # Filter the unsignificant genes
  formula <- y ~ x
  print(ggplot(df, aes(x=evo,y=ovy))+ geom_point(size=2,alpha = 0.7,colour="#325d81")+
          xlab(SMCPC.idents[ident.par[1]]) + ylab(SMCPC.idents[ident.par[2]]) +
          stat_smooth_func(geom="text",method="lm",hjust=0,vjust=1,parse=T,fullrange = T)+
          geom_smooth(method = "lm", se = FALSE,colour="#b3b3b3",alpha=0.8)+
          ggtitle(paste0(SMCPC.idents[ident.par[1]],"_vs_",SMCPC.idents[ident.par[2]]))+theme_bw()+
          stat_fit_glance(method = "lm",
                          method.args = list(formula = formula),
                          label.x = "right",
                          label.y = "bottom",
                          aes(label = paste("italic(P)*\"-value = \"*",
                                            signif(..p.value.., digits = 4), sep = "")),
                          parse = TRUE)+
          scale_x_continuous(breaks=seq(-1,1,0.25),limits = c(-1, 1))+ scale_y_continuous(breaks=seq(-1,1,0.25),limits = c(-1, 1))+
          geom_hline(yintercept = 0) +geom_vline(xintercept = 0))
}
dev.off()

pdf("./Figures/Figure_materials/subpop_coordi_evo.pdf")
for (i in 1:6) {
  ident.par <- cor.ident[i,]
  df1 <- x[[ident.par[1]]]
  df2 <- x[[ident.par[2]]]
  its.genes <- intersect(rownames(df1),rownames(df2))
  sig1=(df1[its.genes, "p_val_adj"]<0.05)*1
  sig2=(df2[its.genes, "p_val_adj"]<0.05)*2
  sig=sig1+sig2
  df <- tibble(evo=df1[its.genes,"avg_logFC"],ovy=df2[its.genes,"avg_logFC"],
               sig=sig ,row.names = its.genes)
  df <- df %>% filter(sig>0) # Filter the unsignificant genes
  formula <- y ~ x
  print(ggplot(df, aes(x=evo,y=ovy))+ geom_point(size=2,alpha = 0.7,colour="#325d81")+
          xlab(AC.sub.ident[ident.par[1]]) + ylab(AC.sub.ident[ident.par[2]]) +
          stat_smooth_func(geom="text",method="lm",hjust=0,vjust=1,parse=T,fullrange = T)+
          geom_smooth(method = "lm", se = FALSE,colour="#b3b3b3",alpha=0.8)+
          ggtitle(paste0(AC.sub.ident[ident.par[1]],"_vs_",AC.sub.ident[ident.par[2]]))+theme_bw()+
          stat_fit_glance(method = "lm",
                          method.args = list(formula = formula),
                          label.x = "right",
                          label.y = "bottom",
                          aes(label = paste("italic(P)*\"-value = \"*",
                                            signif(..p.value.., digits = 4), sep = "")),
                          parse = TRUE)+
          scale_x_continuous(breaks=seq(-1,1,0.25),limits = c(-1, 1))+ scale_y_continuous(breaks=seq(-1,1,0.25),limits = c(-1, 1))+
          geom_hline(yintercept = 0) +geom_vline(xintercept = 0))
}
dev.off()

# Fig S7. Further validation 

ct.num <- tibble(celltype=c("EC","MG","SMC","CPC","PC","OLG","OPC"),
                 b6=c(2901,636,381,186,149,69,40))

ct.num <- mutate(ct.num, prop= b6/sum(ct.num$b6)*100)

pdf(file = "./Revision/Rev1/B6_Celltype_bar_chart.pdf")
ggbarplot(data = ct.num, x = "celltype", y = "prop", fill= "celltype",palette = c( "#4DBBD5FF", "#E64B35FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF"))
dev.off()

library(ggrepel)
dot.labels <- c(mg.marker$Gene, smc.marker$Gene)
degdot.t500pvallabel <- function(x,y){ # Visualized the intersec items between x and y by LogFC in a coordinate colored by adjPvalue
  for (i in 1:length(x)) {
    df1 <- x[[i]]
    df2 <- y[[i]]
    its.genes <- intersect(rownames(df1),rownames(df2))
    sig1=(df1[its.genes, "p_val_adj"]<0.05)*1
    sig2=(df2[its.genes, "p_val_adj"]<0.05)*2
    sig=sig1+sig2
    df <- tibble(batch5_ovy=df1[its.genes,"avg_logFC"],batch6_evo=df2[its.genes,"avg_logFC"]*-1,
                 sig=sig ,row.names = its.genes)
    df <- df %>% filter(sig>0) # Filter the unsignificant genes
    df$pval <- ifelse(df$sig==1,1,-1*log10(df2[its.genes,"p_val_adj"]))
    new.levels <- c("batch5_ovy","batch6_evo","both_sig")
    df$sig <- factor(new.levels[df$sig],levels = c("batch5_ovy","batch6_evo","both_sig"))
    path.gene <- read_csv(paste0(siggene.path,ct.list[i],"_fctop500.csv")) 
    path.gene <- path.gene$value
    #df <- df %>% filter(sig!="evo_sig")
    formula <- y ~ x
    # print(nrow(df))
    # return(nrow(df))
    df$overlap <- ifelse(df$row.names %in% path.gene, "Overlap","None")
    df <- df %>% filter(overlap=="Overlap")
    df$lab <- ifelse(df$row.names %in% dot.labels, df$row.names, '')
    df$sig <- ifelse(df$row.names %in% dot.labels, "sig_sig",new.levels[df$sig])
    print(nrow(df))
    print(ggplot(df, aes(x=batch5_ovy,y=batch6_evo))+ geom_point(alpha = 0.7,aes(size=pval,colour=sig))+
            scale_size_continuous(range = c(1,10),limits = c(1,100))+ scale_color_manual(values = c("batch5_ovy" = "#4885b8", "both_sig" = "#325d81"))+
            geom_text_repel(aes(label = lab),size = 3.5)+
            stat_smooth_func(geom="text",method="lm",hjust=0,vjust=1,parse=T,fullrange = T)+
            geom_smooth(method = "lm", se = FALSE,colour="#b3b3b3",alpha=0.8)+
            ggtitle(paste0(ct.list[i],"_Celltype"))+theme_bw()+
            stat_fit_glance(method = "lm",
                            method.args = list(formula = formula),
                            label.x = "right",
                            label.y = "bottom",
                            aes(label = paste("italic(P)*\"-value = \"*",
                                              signif(..p.value.., digits = 4), sep = "")),
                            parse = TRUE)+
            scale_x_continuous(breaks=seq(-1,1,0.25),limits = c(-1, 1))+ scale_y_continuous(breaks=seq(-1,1,0.25),limits = c(-1, 1))+
            geom_hline(yintercept = 0) +geom_vline(xintercept = 0))
  }
} 

pdf(file = "./Revision/Rev1/B5OVY_B6EVO_T500_p05_psize_label.pdf",width = 8, height = 7)
degdot.t500pvallabel(ovy5.deg,ovy6.deg)
dev.off()

ovy124.deg <- readRDS("./ovy124_deg.rds")
ove124.deg <- readRDS("./ove124_deg.rds")

ovy124.celltype.list <- c("AC","OPC","MG","MAC","OLG","SMC","PC","EC")

degdot.final <- function(x,y,ct.list){ # Visualized the intersec items between x and y by LogFC in a coordinate colored by adjPvalue
  for (i in 1:length(ct.list)) {
    ct <- ct.list[i]
    df1 <- x[[ct]]
    df2 <- y[[ct]]
    its.genes <- intersect(rownames(df1),rownames(df2))
    sig1=(df1[its.genes, "p_val_adj"]<0.05)*1
    sig2=(df2[its.genes, "p_val_adj"]<0.05)*2
    sig=sig1+sig2
    df <- tibble(ovy=df1[its.genes,"avg_logFC"],evo=df2[its.genes,"avg_logFC"]*-1, 
                 sig=sig ,row.names = its.genes)
    df <- df %>% filter(sig>0) # Filter the unsignificant genes
    new.levels <- c("ovy_sig","evo_sig","both_sig")
    df$sig <- factor(new.levels[df$sig],levels = c("ovy_sig","evo_sig","both_sig"))
    #df <- df %>% filter(sig!="evo_sig")
    formula <- y ~ x
    #print(nrow(df))
    # return(nrow(df))
    print(ggplot(df, aes(x=ovy,y=evo))+ geom_point(size=2,alpha = 0.7,colour="#325d81")+
            stat_smooth_func(geom="text",method="lm",hjust=0,vjust=1,parse=T,fullrange = T)+
            geom_smooth(method = "lm", se = FALSE,colour="#b3b3b3",alpha=0.8)+
            ggtitle(paste0(ct.list[i],"_Celltype"))+theme_bw()+
            stat_fit_glance(method = "lm",
                            method.args = list(formula = formula),
                            label.x = "right",
                            label.y = "bottom",
                            aes(label = paste("italic(P)*\"-value = \"*",
                                              signif(..p.value.., digits = 4), sep = "")),
                            parse = TRUE)+
            scale_x_continuous(breaks=seq(-1,1,0.25),limits = c(-1, 1))+ scale_y_continuous(breaks=seq(-1,1,0.25),limits = c(-1, 1))+
            geom_hline(yintercept = 0) +geom_vline(xintercept = 0))
  }
} 

pdf(file = "./Revision/Rev1/B124OVY_EVO_p005.pdf")
degdot.final(ovy124.deg,ove124.deg,ovy124.celltype.list)
dev.off()
