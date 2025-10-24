##################################################
## Project:  DNA methylation
## Script purpose:
## Date: 2024/10/10
## Author: 
## Version: 2
##################################################

library(ggplot2)
library(tidyr)
library(dplyr)
library(ggfortify)  
library(rstatix)
library(ggprism)
library(ggpubr)
library(FactoMineR)
library(factoextra)  
library(RColorBrewer)
library(ggsci)


data <- read.table("input_file", header = TRUE, stringsAsFactors = FALSE)

data <- data %>%
  mutate(across(5:ncol(data), as.numeric))

colnames(data)[1:4] <- c("chromosome", "start", "end", "strand")
sample_cols <- colnames(data)[5:ncol(data)]  

sample_groups <- c(
  "2cell.1.ratio" = "KO 2cell",
  "2cell.2.ratio" = "KO 2cell",
  "WT.2cell.1.ratio" = "WT 2cell",
  "WT.2cell.2.ratio" = "WT 2cell"
)

data_long <- data1 %>%
  pivot_longer(cols = all_of(sample_cols), names_to = "sample", values_to = "ratio")

data_long$group <- sample_groups[data_long$sample]
data_long$ratio <- as.numeric(data_long$ratio)
typeof(data_long$ratio)

meth_stat <- data_long %>% 
  wilcox_test(ratio~group) %>% 
  adjust_pvalue(method = "bonferroni") %>% 
  add_significance() %>% 
  mutate(y.position=1.05)

new_order <- c("WT 2cell", "KO 2cell")

data_long$group <- factor(data_long$group, levels = new_order)

ggviolin_plot <- ggviolin(data_long, x="group", y="ratio", xlab="", ylab="CG methylation of 1Mb bins", fill = "group", palette = c("grey", "#FC4E07"))+ 
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) + ggtitle("Two cell methylation") +
  scale_y_continuous(labels = scales::percent) +
  stat_pvalue_manual(meth_stat %>% filter(group1!="IgG"), hide.ns = FALSE, label = "p={p.adj}" )+
  theme(legend.position = "none",axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(output_violin, plot = ggviolin_plot, width = 3, height = 5, units = "in")

################################# PCA #####################################
dat.pca <- PCA(t(data1[,c(5:8)]) , graph = F) 
percentVar <- get_eig(dat.pca); percentVar[1,'variance.percent']
fviz_screeplot(dat.pca, addlabels = TRUE)
var <- get_pca_var(dat.pca) 
head(sort(var$contrib[,2],decreasing = T)) 
fviz_contrib(dat.pca, choice ="var", axes = 1:2, top = 10)
var_contrib <- as.data.frame(subset(var$contrib,select=c(1:2)))
var_contrib$sum <- var_contrib$Dim.1*percentVar[1,'variance.percent']*0.01+var_contrib$Dim.2*percentVar[2,'variance.percent']*0.01
head(rownames(var_contrib)[order(var_contrib$sum,decreasing = T)],10)

pca <- fviz_pca_ind(dat.pca,
                    title = "Principal Component Analysis",
                    legend.title = "Groups",
                    geom.ind = c("point", "text"), #"point", "text"
                    pointsize = 1.5,
                    label="none",
                    # labelsize = 4,
                    # repel = TRUE, 
                    col.ind = sample_groups, 
                    addEllipses = TRUE, 
                    axes.linetype=NA,  
                    mean.point=F
) +
  
  coord_fixed(ratio = 1)+ 
  xlab(paste0("PC1 (",round(percentVar[1,'variance.percent'],1),"%)")) +
  ylab(paste0("PC2 (",round(percentVar[2,'variance.percent'],1),"%)")) + 
  theme_classic()+
  # theme(legend.position = "none")+  # "none" REMOVE legend
  # ggforce::geom_mark_ellipse(aes(fill = group_list,color = group_list)) +  
  scale_color_lancet()+
  theme(
    legend.key = element_rect(fill = NA, color = NA),  
    legend.position = "right",     
    legend.background = element_blank(),  
    legend.box.background = element_blank()  
  )


ggsave(pca, filename = output_pca ,width = 5,height =4)
dev.off()


