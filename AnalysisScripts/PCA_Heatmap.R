##Log: mRNA PCs + miRNA BS/miRNA
##Implementation of Option 5

##First, derive miRNA signatures, based on both numBS and Expr
library(survival)
library(survminer)
library(stringr)
source("/ihome/hpark/yub20/SUR3UTR/GDC_TCGA_FPKM_analysis/independence/functions_independence.R")

##Read in data
mRNA_dir <- "/ihome/hpark/yub20/SUR3UTR/data/TCGA_harmonized/FPKM"
surv_dir <- "/ihome/hpark/yub20/SUR3UTR/new/survival"
APA4_dir <- "/ihome/hpark/yub20/SUR3UTR/pythons/TCGA_FPKM/TC3A/result_TC3A"
Ctrl4_dir <- "/ihome/hpark/yub20/SUR3UTR/pythons/TCGA_FPKM/TC3A/result_TC3A"
PDUI_dir <- "/ihome/hpark/yub20/SUR3UTR/new/data"


cancerType1 <- c("BRCA","SKCM","LUAD","HNSC","KIRC", "LGG","LUSC","OV","STAD","UCEC")

cancertype <- cancerType1[1]
print(paste0("Start Processing ", cancertype))
data_list <- readin_clinical(cancertype)
mRNA.tab <- data_list$mRNA
numBS.tab <- data_list$numBS
miRNA.tab <- data_list$miRNA
PDUI.tab <- data_list$PDUI
cancer.tab <- data.frame(patID=mRNA.tab$patID, cancer=rep(cancertype, nrow(mRNA.tab)))

rbind_table <- function(tab, tmp){
  tab <- tab[,colnames(tab) %in% colnames(tmp)]
  tmp <- tmp[,colnames(tmp) %in% colnames(tab)]
  tmp <- tmp[,match(colnames(tab), colnames(tmp))]
  return(rbind(tab, tmp))
}

for (i1 in 2:length(cancerType1)){
  cancertype <- cancerType1[i1]
  print(paste0("Start Processing ", cancertype))
  data_list <- readin_clinical(cancertype)
  
  mRNA.tmp <- data_list$mRNA
  numBS.tmp <- data_list$numBS
  miRNA.tmp <- data_list$miRNA
  PDUI.tmp <- data_list$PDUI
  cancer.tmp <- data.frame(patID=mRNA.tmp$patID, cancer=rep(cancertype, nrow(mRNA.tmp)))
  
  mRNA.tab <- rbind_table(mRNA.tab, mRNA.tmp)
  numBS.tab <- rbind_table(numBS.tab, numBS.tmp)
  miRNA.tab <- rbind_table(miRNA.tab, miRNA.tmp)
  PDUI.tab <- rbind_table(PDUI.tab, PDUI.tmp)
  cancer.tab <- rbind(cancer.tab, cancer.tmp)
}

groups <- as.factor(cancer.tab$cancer)

## calculate within variance/ total variance
mRNA.mat <- mRNA.tab[,-1]
numBS.mat <- numBS.tab[,-1]
miRNA.mat <- miRNA.tab[,-1]
PDUI.mat <- PDUI.tab[,-1]
names(groups) <- cancer.tab$patID
between_rate <- numeric(length = 4)
within_rate <- numeric(length = 4)
names(between_rate) <- c("mRNA","numTS","expr","PDUI")
names(within_rate) <- c("mRNA","numTS","expr","PDUI")
# mRNA
center <- sapply(base::split(mRNA.mat, groups), colMeans)
ss <- function(x) sum(scale(x, scale = FALSE)^2)
betweenss <- ss(t(center)[groups,])

withinss <- sapply(base::split(mRNA.mat, groups), ss)
tot.withinss <- sum(withinss) # or  resid <- x - fitted(cl); ss(resid)

totss <- ss(mRNA.mat)

between_rate[1] <- betweenss/totss
within_rate[1] <- tot.withinss/totss

# numBS
center <- sapply(base::split(numBS.mat, groups), colMeans)
betweenss <- ss(t(center)[groups,])

withinss <- sapply(base::split(numBS.mat, groups), ss)
tot.withinss <- sum(withinss) # or  resid <- x - fitted(cl); ss(resid)

totss <- ss(numBS.mat)

between_rate[2] <- betweenss/totss
within_rate[2] <- tot.withinss/totss

# expr
center <- sapply(base::split(miRNA.mat, groups), colMeans)
betweenss <- ss(t(center)[groups,])

withinss <- sapply(base::split(miRNA.mat, groups), ss)
tot.withinss <- sum(withinss) # or  resid <- x - fitted(cl); ss(resid)

totss <- ss(miRNA.mat)

between_rate[3] <- betweenss/totss
within_rate[3] <- tot.withinss/totss

# PDUI
center <- sapply(base::split(PDUI.mat, groups), colMeans)
betweenss <- ss(t(center)[groups,])

withinss <- sapply(base::split(PDUI.mat, groups), ss)
tot.withinss <- sum(withinss) # or  resid <- x - fitted(cl); ss(resid)

totss <- ss(PDUI.mat)

between_rate[4] <- betweenss/totss
within_rate[4] <- tot.withinss/totss

plot_tab <- data.frame(feature = names(between_rate), rate = between_rate)
plot_tab$feature <- factor(plot_tab$feature, levels = c("numTS", "expr", "PDUI", "mRNA"))
ggplot(plot_tab, aes(x = feature, y=rate)) + geom_bar(stat="identity", fill = "black") +theme_bw()

## drawing PCA
library(factoextra)
library(dplyr)
#library(ggbiplot)
data.frame(libSize = rowSums(mRNA.tab[,-1]),groups) %>% group_by(groups) %>% dplyr::summarise(mean(libSize))
data.frame(libSize = rowSums(numBS.tab[,-1]),groups) %>% group_by(groups) %>% dplyr::summarise(mean(libSize))
data.frame(libSize = rowSums(miRNA.tab[,-1]),groups) %>% group_by(groups) %>% dplyr::summarise(mean(libSize))
data.frame(libSize = rowSums(PDUI.tab[,-1]),groups) %>% group_by(groups) %>% dplyr::summarise(mean(libSize))


library(RColorBrewer)
display.brewer.pal(n = 10, name = 'Set3')
color_vec <- brewer.pal(n = 10, name = 'Set3')
names(color_vec) <- names(table(groups))


mRNA.pca <- prcomp(mRNA.tab[,-1],center = T, scale. = F)
fviz_eig(mRNA.pca)
mRNA_plot <- fviz_pca_ind(mRNA.pca, geom="point", pointsize = 3,
                          col.ind = groups,legend.title = "Cancer",
                          label = "none", title = "PCA based on mRNA expression")
mRNA_plot <- fviz_pca_ind(mRNA.pca, geom="point", pointsize = 1.5,alpha.ind = 0.8,
                          col.ind = groups,legend.title = "Cancer",
                          label = "none",axes.linetype=NA, invisible="quali")+labs(title =element_blank(), x = "PC1 (19.8%)", y = "PC2 (8.3%)")+ theme_minimal()+
  scale_shape_manual(values=rep(19,10)) + theme_classic()+ theme(panel.grid.minor = element_blank(),axis.text.x = element_blank(),
                                                                 axis.text.y = element_blank(),axis.ticks = element_blank())+
  scale_y_continuous(limits = c(-120, 120))+
  scale_x_continuous(limits = c(-60, 160))
mRNA_plot+ scale_color_manual(values=c("#faa0d4", "#da950a", "#e85bf5" ,"#a6a90a", "#9994ff", "#41b90a", "#0ab4f7", "#1b6df2", "#fa766d","orange"))
fviz_pca_var(mRNA.pca)
fviz_pca_biplot(mRNA.pca, repel = TRUE, col.ind = groups,legend.title = "Cancer",ellipse.type = "confidence",
                title = "PCA based on mRNA expression", select.var = list(contrib = 10))

numBS.pca <- prcomp(na.omit(numBS.tab)[,-1],center = T, scale. = F)
fviz_eig(numBS.pca)
numBS_plot <- fviz_pca_ind(numBS.pca, geom="point", pointsize = 1.5,alpha.ind = 0.8,
                           col.ind = groups[!is.na(rowSums(numBS.tab[,-1]))],legend.title = "Cancer",
                           label = "none",axes.linetype=NA, invisible="quali")+labs(title =element_blank(), x = "PC1 (58.8%)", y = "PC2 (9.2%)")+ theme_minimal()+
  scale_shape_manual(values=rep(19,10)) + theme_classic()+ theme(panel.grid.minor = element_blank(),axis.text.x = element_blank(),
                                                                 axis.text.y = element_blank(),axis.ticks = element_blank())+
  scale_y_continuous(limits = c(-12, 28), breaks = c(-10, 0, 10,20))#+scale_color_manual(values=c("#006400", "#00008b", "#b03060", "#ff4500", "#ffd700", "#7cfc00", "#00ffff", "#ff00ff", "#6495ed", "#ffdab9"))
numBS_plot + scale_color_manual(values=c("#faa0d4", "#da950a", "#e85bf5" ,"#a6a90a", "#9994ff", "#41b90a", "#0ab4f7", "#1b6df2", "#fa766d","orange"))
numBS_plot + scale_color_manual(values=color_vec)
numBS_plot
fviz_pca_biplot(numBS.pca, repel = TRUE, col.ind = groups,legend.title = "Cancer",ellipse.type = "confidence",
                title = "PCA based on number of target sites", select.var = list(contrib = 10))

coefVar <- apply(numBS.tab[,-1], 2, function(x){sd(x)/mean(x)})
sort(coefVar, decreasing = T)[1:20]
#data_list$mRNA <- data_list$mRNA[, c("patID", names(sort(coefVar, decreasing = T)[1:2000]))]


miRNA.pca <- prcomp(na.omit(miRNA.tab)[,-1],center = T, scale. = F)
fviz_eig(miRNA.pca)
miRNA_plot <- fviz_pca_ind(miRNA.pca,
                           col.ind = groups[!is.na(rowSums(numBS.tab[,-1]))],legend.title = "Cancer",ellipse.type = "confidence",
                           label = "none", title = "PCA based on miRNA expression")
miRNA_plot <- fviz_pca_ind(miRNA.pca, geom="point", pointsize = 1.5,alpha.ind = 0.8,
                           col.ind = groups[!is.na(rowSums(numBS.tab[,-1]))],legend.title = "Cancer",
                           label = "none",axes.linetype=NA, invisible="quali")+labs(title =element_blank(), x = "PC1 (35.1%)", y = "PC2 (11.7%)")+ theme_minimal()+
  scale_shape_manual(values=rep(19,10)) + theme_classic()+ theme(panel.grid.minor = element_blank(),axis.text.x = element_blank(),
                                                                 axis.text.y = element_blank(),axis.ticks = element_blank())+
  scale_y_continuous(limits = c(-27, 49))+
  scale_x_continuous(limits = c(-55, 38))
miRNA_plot+ scale_color_manual(values=c("#faa0d4", "#da950a", "#e85bf5" ,"#a6a90a", "#9994ff", "#41b90a", "#0ab4f7", "#1b6df2", "#fa766d","orange"))

fviz_pca_biplot(miRNA.pca, repel = TRUE, col.ind = groups,legend.title = "Cancer",ellipse.type = "confidence",
                title = "PCA based on miRNA expression", select.var = list(contrib = 10))

PDUI.pca <- prcomp(PDUI.tab[,-1],center = T, scale. = F)
fviz_eig(PDUI.pca)
PDUI_plot <- fviz_pca_ind(PDUI.pca,
                          col.ind = groups,legend.title = "Cancer",ellipse.type = "confidence",
                          label = "none", title = "PCA based on PDUI")
PDUI_plot
PDUI_plot <- fviz_pca_ind(PDUI.pca, geom="point", pointsize = 1.5,alpha.ind = 0.8,
                          col.ind = groups,legend.title = "Cancer",
                          label = "none",axes.linetype=NA, invisible="quali")+labs(title =element_blank(), x = "PC1 (39.6%)", y = "PC2 (14.3%)")+ theme_minimal()+
  scale_shape_manual(values=rep(19,10)) + theme_classic()+ theme(panel.grid.minor = element_blank(),axis.text.x = element_blank(),
                                                                 axis.text.y = element_blank(),axis.ticks = element_blank())+
  scale_y_continuous(limits = c(-2.1, 2.8))+
  scale_x_continuous(limits = c(-2.2, 3.9))
PDUI_plot + scale_color_manual(values=c("#faa0d4", "#da950a", "#e85bf5" ,"#a6a90a", "#9994ff", "#41b90a", "#0ab4f7", "#1b6df2", "#fa766d","orange"))
fviz_pca_biplot(PDUI.pca, repel = TRUE, col.ind = groups,legend.title = "Cancer",ellipse.type = "confidence",
                title = "PCA based on PDUI", select.var = list(contrib = 10))

## Distance, mRNA
## mahanobis distance between centroids
mRNA.rotate <- scale(mRNA.pca$x)
mRNA.rotate <- as.data.frame(mRNA.rotate)
identical(names(groups), rownames(mRNA.rotate))
mRNA.rotate$groups <- groups
mRNA.rotate.center <- mRNA.rotate %>% group_by(groups) %>% 
  summarise_at(vars(matches("PC")), list(name = mean)) 
mRNA.rotate.center <- as.data.frame(mRNA.rotate.center)
rownames(mRNA.rotate.center) <- mRNA.rotate.center[,1]
mRNA.rotate.center <- mRNA.rotate.center[,-1]
mRNA.maha <- dist(mRNA.rotate.center, method = "euclidean")
mean(mRNA.maha)
## mean pairwise distance maha: 8.253436 
mRNA.maha <-as.data.frame(as.matrix(mRNA.maha))
write.csv(mRNA.maha,file = "mRNA.maha.csv")
## euclidean distance of first two PCs
mRNA.2PC <- dist(mRNA.rotate.center[,1:2], method = "euclidean")
mean(mRNA.2PC)
## mean pairwise distance 2PC: 1.524494 
mRNA.2PC <-as.data.frame(as.matrix(mRNA.2PC))
write.csv(mRNA.2PC,file = "mRNA.2PC.csv")

## Distance, numBS
## mahanobis distance between centroids
numBS.rotate <- scale(numBS.pca$x)
numBS.rotate <- as.data.frame(numBS.rotate)
identical(names(groups[!is.na(rowSums(numBS.tab[,-1]))]), rownames(numBS.rotate))
numBS.rotate$groups <- groups[!is.na(rowSums(numBS.tab[,-1]))]
numBS.rotate.center <- numBS.rotate %>% group_by(groups) %>% 
  summarise_at(vars(matches("PC")), list(name = mean)) 
numBS.rotate.center <- as.data.frame(numBS.rotate.center)
rownames(numBS.rotate.center) <- numBS.rotate.center[,1]
numBS.rotate.center <- numBS.rotate.center[,-1]
numBS.maha <- dist(numBS.rotate.center, method = "euclidean")
mean(numBS.maha)
## mean pairwise distance maha: 7.500251 
numBS.maha <-as.data.frame(as.matrix(numBS.maha))
write.csv(numBS.maha,file = "numBS.maha.csv")
## euclidean distance of first two PCs
numBS.2PC <- dist(numBS.rotate.center[,1:2], method = "euclidean")
mean(numBS.2PC)
## mean pairwise distance 2PC: 1.66492 
numBS.2PC <-as.data.frame(as.matrix(numBS.2PC))
write.csv(numBS.2PC,file = "numBS.2PC.csv")

## Distance, miRNA
## mahanobis distance between centroids
miRNA.rotate <- scale(miRNA.pca$x)
miRNA.rotate <- as.data.frame(miRNA.rotate)
identical(names(groups[!is.na(rowSums(miRNA.tab[,-1]))]), rownames(miRNA.rotate))
miRNA.rotate$groups <- groups[!is.na(rowSums(miRNA.tab[,-1]))]
miRNA.rotate.center <- miRNA.rotate %>% group_by(groups) %>% 
  summarise_at(vars(matches("PC")), list(name = mean)) 
miRNA.rotate.center <- as.data.frame(miRNA.rotate.center)
rownames(miRNA.rotate.center) <- miRNA.rotate.center[,1]
miRNA.rotate.center <- miRNA.rotate.center[,-1]
miRNA.maha <- dist(miRNA.rotate.center, method = "euclidean")
mean(miRNA.maha)
## mean pairwise distance maha: 6.391537 
miRNA.maha <-as.data.frame(as.matrix(miRNA.maha))
write.csv(miRNA.maha,file = "miRNA.maha.csv")
## euclidean distance of first two PCs
miRNA.2PC <- dist(miRNA.rotate.center[,1:2], method = "euclidean")
mean(miRNA.2PC)
## mean pairwise distance 2PC: 1.544147
miRNA.2PC <-as.data.frame(as.matrix(miRNA.2PC))
write.csv(miRNA.2PC,file = "miRNA.2PC.csv")

## Distance, PDUI
## mahanobis distance between centroids
PDUI.rotate <- scale(PDUI.pca$x)
PDUI.rotate <- as.data.frame(PDUI.rotate)
identical(names(groups), rownames(PDUI.rotate))
PDUI.rotate$groups <- groups
PDUI.rotate.center <- PDUI.rotate %>% group_by(groups) %>% 
  summarise_at(vars(matches("PC")), list(name = mean)) 
PDUI.rotate.center <- as.data.frame(PDUI.rotate.center)
rownames(PDUI.rotate.center) <- PDUI.rotate.center[,1]
PDUI.rotate.center <- PDUI.rotate.center[,-1]
PDUI.maha <- dist(PDUI.rotate.center, method = "euclidean")
mean(PDUI.maha)
## mean pairwise distance maha: 6.196824 
PDUI.maha <-as.data.frame(as.matrix(PDUI.maha))
write.csv(PDUI.maha,file = "PDUI.maha.csv")
## euclidean distance of first two PCs
PDUI.2PC <- dist(PDUI.rotate.center[,1:2], method = "euclidean")
mean(PDUI.2PC)
## mean pairwise distance 2PC: 1.354765
PDUI.2PC <-as.data.frame(as.matrix(PDUI.2PC))
write.csv(PDUI.2PC,file = "PDUI.2PC.csv")


## heatmap
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
Heatmap(as.matrix(numBS.tab[,-1]), col = colorRamp2(c(-100, 0, 100), brewer.pal(n=3, name="RdBu")), 
        column_order = colnames(numBS.tab[,-1])[order(colSums(numBS.tab[,-1]>0), decreasing = T)],
        show_row_names = F, show_column_names = F)

t_numBS.tab <- t(numBS.tab[,-1])

library(RColorBrewer)
display.brewer.pal(n = 10, name = 'Set3')
color_vec <- brewer.pal(n = 10, name = 'Set3')
#color_vec <- c("#ff68be", "#da950a", "#e970f4" ,"#a6a90a", "#9994ff", "#41b90a", "#0ab4f7", "#0ac282", "#fa766d","#0ac1c7")
color_vec <- c("#faa0d4", "#da950a", "#e85bf5" ,"#a6a90a", "#9994ff", "#41b90a", "#0ab4f7", "#1b6df2", "#fa766d","orange")
names(color_vec) <- names(table(groups))

ha <- HeatmapAnnotation(
  cancer= groups,
  col = list(cancer = color_vec)
)


scale_rows <- function (x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m)/s)
}

scaled_t_numBS.tab <- scale_rows(t_numBS.tab)

Heatmap(scaled_t_numBS.tab, col = colorRamp2(c(5, 0, -5), c("red","white","blue")), top_annotation = ha,
        show_row_names = F, show_column_names = F, cluster_columns = T,clustering_method_rows = "average", clustering_method_columns="average")

scaled_t_miRNA.tab <- scale_rows(t(miRNA.tab[,-1]))
Heatmap(scaled_t_miRNA.tab, col = colorRamp2(c(5, 0, -5), c("red","white","blue")), top_annotation = ha,
        show_row_names = F, show_column_names = F, cluster_columns = T,clustering_method_rows = "average", clustering_method_columns="average")

scaled_t_mRNA.tab <- scale_rows(t(mRNA.tab[,-1]))
Heatmap(scaled_t_mRNA.tab, col = colorRamp2(c(5, 0, -5), c("red","white","blue")), top_annotation = ha,
        show_row_names = F, show_column_names = F, cluster_columns = T,clustering_method_rows = "average", clustering_method_columns="average")
