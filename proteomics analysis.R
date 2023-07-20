# 
install.packages("EnrichmentMap")

if (!require("BiocManager", quietly = TRUE))
  install.packages("EnrichmentMap")

BiocManager::install("clusterProfiler")

library(RColorBrewer)
library(UpSetR)
library(VennDiagram)
library(dplyr)
library(eulerr)
library(factoextra)
library(forcats)
library(gg.gap)
library(ggVolcano)
library(ggalt)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggsignif)
library(hrbrthemes)
library(limma)
library(openxlsx)
library(pheatmap)
library(reshape2) 
library(scales)
library(tidyr)
library(tidyverse)          
library(venn)
library(venneuler)
library(viridis)
library(xlsx)
library(DESeq2)
# library(clusterProfiler)
library(org.Hs.eg.db)
# library(EnrichmentMap)
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(readr)
library("genefilter")
library(pheatmap)
# library(DO.db)
# require(DOSE)
# library(AnnotationHub)
library(DOSE)
library(GO.db)
library(topGO)
library(GSEABase)
library(clusterProfiler)
############## Data load ###########################################################################################

setwd('D:/yunxia/proteomics data analysis/')
gd1_long <- read.xlsx('test_data.xlsx', 2, detectDates = TRUE,fillMergedCells = TRUE)
gd1_long[,1:5]

############## Start: Data distribution ###########################################################################################
plot_data <- gd1_long[,-3]

plot_data[1:5,1:5]

plot_data01 <- melt(plot_data,id = c("uniprot","batch"))

head(plot_data01)


p1 <- ggplot(plot_data01,aes(x=uniprot,y=value,fill=batch))+
  geom_boxplot(width=0.6,alpha=0.8)

p1

t<-unclass(Sys.time())
filename = paste("Data distribution.png",sep = "")
setwd('D:/yunxia/proteomics data analysis/Result_png/')
png(file=filename, width = 2000, height = 1400,res = 300)
print(p1)
dev.off()

############## End: Data distribution ###########################################################################################

############## Start: PCA ###########################################################################################

gd1_long_num <- gd1_long[,-c(1,2,3)]
gd1_long_num_t <- t(gd1_long_num)

dim(data_all)
gd1_long_num_t_no_na <- t(na.omit(gd1_long_num_t))
dim(gd1_long_num_t_no_na)

data_all <- cbind(gd1_long[,c(1,2,3)],gd1_long_num_t_no_na)


datanames <- unique(data_all[,'label'])

data_combine <- NULL

num_nm = 2
for (num_nm in c(1:length(datanames))){
  
  datanm <- datanames[num_nm]
  dataraw <- data_all[data_all[,'label'] == datanm,]

  data_re <- transform(dataraw, uniprot = sub("[12]$", "", uniprot))[,-c(2,3)]
  
  
  
  data_arrag <- aggregate(. ~ uniprot, data = data_re, mean)
  data_arrag[,c(1:10)]
  data_combine <- rbind(data_combine,data_arrag)
  data_combine[,c(1:10)]
  
}


gd1_long.pca <- prcomp(data_combine[,-1])          



P <- fviz_pca_ind(gd1_long.pca,          
             label = "PCA",         
             habillage = data_combine$uniprot,          
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),          
             addEllipses = TRUE        
)


filename = paste("PCA for samples.png",sep = "")
setwd('D:/yunxia/proteomics data analysis/Result_png/')
png(file=filename, width = 2000, height = 2100,res = 300)
print(P)
dev.off()

############## End: PCA ###########################################################################################

############## Start: calculate ratios ###########################################################################################

gd1_long_num <- gd1_long[,-c(1,2,3)]
gd1_long_num_t <- t(gd1_long_num)

dim(data_all)
gd1_long_num_t_no_na <- t(na.omit(gd1_long_num_t))
dim(gd1_long_num_t_no_na)

data_all <- cbind(gd1_long[,c(1,2,3)],gd1_long_num_t_no_na)


datanames <- unique(data_all[,'label'])

data_combine <- NULL

num_nm = 1
for (num_nm in c(1:length(datanames))){
  
  datanm <- datanames[num_nm]
  dataraw <- data_all[data_all[,'label'] == datanm,]
  
  data_re <- transform(dataraw, uniprot = sub("[12]$", "", uniprot))[,-c(2,3)]

  # get the average of expanded T cells
  data_arrag <- aggregate(. ~ uniprot, data = data_re, mean)
  data_arrag[,c(1:10)]
  
  data_ratio <- NULL
  for (nrownum in c(1:5)){
    data_new <- dataraw[nrownum,-c(1,2,3)] / data_arrag[2,-1]
    data_ratio <- rbind(data_ratio,data_new)
    
  }
  
  data_ratio <- cbind(dataraw[,c(1,2,3)],data_ratio)
  
  
  data_combine <- rbind(data_combine,data_ratio)
  data_combine[,c(1:10)]
  
  
}

# get the ratios of each experimental channel to the average of expanded channels within each donor
write.csv(data_combine,file = "D:/yunxia/proteomics data analysis/Result_data/the ratios of each experimental channel.csv")

# get the mean ratio between the donors
data_ratios <- data_combine[,-c(1,2)]
data_mean_ration[,c(1:10)]
data_mean_ration <- aggregate(. ~ label, data = data_ratios, mean)

write.csv(data_mean_ration,file = "D:/yunxia/proteomics data analysis/Result_data/the mean ratio between the donors.csv")

############## End: calculate ratios ###########################################################################################

############## Start: differential expression analysis ###########################################################################################

# get protein expression only expanded T cells and activated T cells
gd1_long_num <- gd1_long[-c(1,6,11,16),-c(1,2,3)] 
gd1_long_num_t <- t(gd1_long_num)
dim(gd1_long_num_t)

# remove the protein that has one miss value
gd1_long_num_t_no_na <- na.omit(gd1_long_num_t)
head(gd1_long_num_t_no_na)
dim(gd1_long_num_t_no_na)

# get sample label
sample_data <-  gd1_long[-c(1,6,11,16),c(1,2,3)] 
sample_data_label <- transform(sample_data, uniprot = sub("[12]$", "", uniprot))


# log2 transformation
exp6 <- log2(gd1_long_num_t_no_na) #log
exp6[exp6 == -Inf] <- 0 #将log化后的负无穷值替换为0

# make sample group label
design=model.matrix(~0+factor(sample_data_label$uniprot))
colnames(design)=levels(factor(sample_data_label$uniprot))
row.names(design) <- colnames(exp6)
design


# Difference comparison matrix
contrast_matrix=makeContrasts(paste0(c('act_','exp_'),collapse = '-'),levels = design)
contrast_matrix  #-1 and 1 mean that expanded T cells is used to be compared, activated T cells is used to be compared

# ----------------------------------------------------------------------------------------
# step 1:lmFit
fit=lmFit(exp6,design)
fit2=contrasts.fit(fit,contrast_matrix)


# step 2:eBayes
fit3=eBayes(fit2)


# step 3:topTable
tempoutput=topTable(fit3,coef = 1,n=Inf)
DEG=na.omit(tempoutput)  #The difference analysis matrix is obtained, focusing on logFC and p-values
head(DEG)  
write.csv(DEG,'D:/yunxia/proteomics data analysis/Result_data/DEG_M.csv')

# ----------------------------------------------------------------------------------------
# Select a threshold for the difference multiple

logFC_cutoff=1
adj_P_Val <- 0.05

# Label the genes in the differential analysis matrix data
DEG$change=as.factor(ifelse(DEG$adj.P.Val < adj_P_Val & abs(DEG$logFC)>logFC_cutoff,
                            ifelse(DEG$logFC>logFC_cutoff,'UP','DOWN'),'NOT'))
head(DEG)

nrow(DEG[DEG$change=='UP',])  
nrow(DEG[DEG$change=='DOWN',])  

thttile=paste0('\nCutoff for log2FC is ',round(logFC_cutoff,1),
               '\nThe number of up gene is ',nrow(DEG[DEG$change=='UP',]),
               '\nThe number of down gene is ',nrow(DEG[DEG$change=='DOWN',]))
thttile  

g=ggplot(data = DEG,aes(x=logFC,y=-log10(adj.P.Val),color=change))+
  geom_point(alpha=0.4,size=1.75)+
  theme_set(theme_set(theme_bw(base_size = 20)))+
  xlab('log2 fold change')+ylab('-log10 adj.P.Val')+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(),
                   axis.line = element_line(colour = 'black'))+
  ggtitle(thttile)+theme(plot.title = element_text(size = 15,hjust = 0.5))+
  scale_color_manual(values = c('blue','black','red'))
print(g)  

setwd('D:/yunxia/proteomics data analysis/Result_png/')
filename = paste("Volcano plot for samples.png",sep = "")
png(file=filename, width = 4000, height = 3000,res = 300)
print(g)
dev.off()

# Differential expression gene expression data
head(DEG)
rownames(DEG[DEG$change=='UP' | DEG$change=='DOWN',])
DEG_MATRIX=exp6[row.names(exp6) %in% rownames(DEG[DEG$change=='UP' | DEG$change=='DOWN',]),]
dim(DEG_MATRIX)  
head(DEG_MATRIX)

# Up-regulated gen  e expression data
UP_DEG=exp6[row.names(exp6) %in% rownames(DEG[DEG$change=='UP',]),]
dim(UP_DEG)  

# Down-regulated gene expression data
DOWN_DEG=exp6[row.names(exp6) %in% rownames(DEG[DEG$change=='DOWN',]),]
dim(DOWN_DEG)  

write.csv(DEG_MATRIX,'D:/yunxia/proteomics data analysis/Result_data/DEG_MATRIX.csv')
write.csv(UP_DEG,'D:/yunxia/proteomics data analysis/Result_data/UP_DEG.csv')
write.csv(DOWN_DEG,'D:/yunxia/proteomics data analysis/Result_data/DOWN_DEG.csv')

############## End: differential expression analysis ###########################################################################################

############## Start: perform GO-term enrichment analysis on the changing proteins ###########################################################################################


gene_da<-rownames(DEG_MATRIX)
head(gene_da)
ego<-enrichGO(gene = gene_da,keyType = "UNIPROT",
              OrgDb = org.Hs.eg.db,
              ont="all",
              pAdjustMethod = "BH",
              universe = background_uniprot,
              pvalueCutoff = 0.05,
              qvalueCutoff = 0.1,
              readable =TRUE
)
head(ego)
dim(ego)

write.csv(ego,file="D:/yunxia/proteomics data analysis/Result_data/EnrichGO_result.csv")


# bar plot
p <- barplot(ego, showCategory=20,drop=T)
setwd('D:/yunxia/proteomics data analysis/Result_png/')
filename = paste("Bar plot for enrichGO.png",sep = "")
png(file=filename, width = 4000, height = 3000,res = 300)
print(p)
dev.off()

# dot plot
p <- dotplot(ego,showCategory=20) 
setwd('D:/yunxia/proteomics data analysis/Result_png/')
filename = paste("Dot plot for enrichGO.png",sep = "")
png(file=filename, width = 4000, height = 3000,res = 300)
print(p)
dev.off()

############## Start: perform GO-term enrichment analysis on the changing proteins ###########################################################################################




