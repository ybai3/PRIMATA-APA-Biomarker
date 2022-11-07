setwd("~/SUR3UTR/GDC_TCGA_FPKM_analysis/externalValid")
load("immune_BT_LUAD_immune.Rdata")
library(dplyr)

## perform hypergeometric test
PS_quantile <- 0.5
CCLE_miRNA <- rownames(Data$result_miRNA_case)[Data$result_miRNA_case$MiRNA_case1_percentage >= quantile(Data$result_miRNA_case$MiRNA_case1_percentage, PS_quantile)]
TCGA_miRNA <- rownames(Data$result_miRNA_ctrl)[Data$result_miRNA_ctrl$MiRNA_ctrl1_percentage >= quantile(Data$result_miRNA_ctrl$MiRNA_ctrl1_percentage, PS_quantile)]
s_sample <- sum(CCLE_miRNA %in% TCGA_miRNA)
s_TCGA <- length(TCGA_miRNA)
f_TCGA <- nrow(Data$result_miRNA_ctrl) - length(TCGA_miRNA)
s_size <- length(CCLE_miRNA)
hyper_pval <- phyper(s_sample-1, s_TCGA, f_TCGA, s_size, lower.tail= FALSE)

library(VennDiagram)
temp <- venn.diagram(
  x = list(CCLE_miRNA, TCGA_miRNA),
  category.names = c("Seo.etal." , "TCGA LUAD"),
  filename = NULL,
  output=F
)
grid.draw(temp)
library(grDevices)
pdf(file="Immune.pdf")
grid.draw(temp)
dev.off()

load("proliferation_BT_LUAD_prolif.Rdata")
## perform hypergeometric test 
PS_quantile <- 0.5
CCLE_miRNA <- rownames(Data$result_miRNA_case)[Data$result_miRNA_case$MiRNA_case1_percentage >= quantile(Data$result_miRNA_case$MiRNA_case1_percentage, PS_quantile)]
TCGA_miRNA <- rownames(Data$result_miRNA_ctrl)[Data$result_miRNA_ctrl$MiRNA_ctrl1_percentage >= quantile(Data$result_miRNA_ctrl$MiRNA_ctrl1_percentage, PS_quantile)]
s_sample <- sum(CCLE_miRNA %in% TCGA_miRNA)
s_TCGA <- length(TCGA_miRNA)
f_TCGA <- nrow(Data$result_miRNA_ctrl) - length(TCGA_miRNA)
s_size <- length(CCLE_miRNA)
hyper_pval <- phyper(s_sample-1, s_TCGA, f_TCGA, s_size, lower.tail= FALSE)


temp <- venn.diagram(
  x = list(CCLE_miRNA, TCGA_miRNA),
  category.names = c("Seo.etal." , "TCGA LUAD"),
  filename = NULL,
  output=F
)
grid.draw(temp)
library(grDevices)
pdf(file="Prolif.pdf")
grid.draw(temp)
dev.off()
