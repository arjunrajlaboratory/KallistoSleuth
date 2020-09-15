### functions_kallisto_sleuth.R

# Packages
if (!requireNamespace("BiocManager", quietly = TRUE)){install.packages("BiocManager")}
library("BiocManager")
if (!requireNamespace("biomaRt", quietly = TRUE)){BiocManager::install("biomaRt")}
library("biomaRt")
if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
library("devtools")
if (!requireNamespace("rhdf5", quietly = TRUE)){BiocManager::install("rhdf5")}
library("rhdf5")
if (!requireNamespace("sleuth", quietly = TRUE)){BiocManager::install("pachterlab/sleuth")}
library("sleuth")
if (!requireNamespace("ggplot2", quietly = TRUE)){install.packages("ggplot2")}
library("ggplot2")
if (!requireNamespace("dplyr", quietly = TRUE)){devtools::install_github("hadley/dplyr")}
library("dplyr")
if (!requireNamespace("ggdendro", quietly = TRUE)){install.packages("ggdendro")}
library("ggdendro")
if (!requireNamespace("reshape2", quietly = TRUE)){install.packages("reshape2")}
library("reshape2")

# List of functions called by kallisto_sleuth.R workook.

# generates a heatmap that shows the gene names on y-axis.
# uses clustering on genes and samples to order the plot.
generate.heatmap = function(sleuthObject = so, sleuth_table = sleuth_table_tx, topNgenes=25){
  sleuth_table = sleuth_table[1:topNgenes,]
  merge.df= data.frame(merge(sleuthObject$obs_raw, sleuth_table, by = "target_id"))
  
  clust.mtx = acast(merge.df, ext_gene~sample, fun.aggregate = sum,value.var="tpm")
  clust.genes.dendro <- as.dendrogram(hclust(d = dist(x = clust.mtx)))
  clust.samples.dendro <- as.dendrogram(hclust(d = dist(x = t(clust.mtx))))
  
  genes.order = order.dendrogram(clust.genes.dendro)
  sample.order = order.dendrogram(clust.samples.dendro)
  
  merge.df$genes.ord = factor(x = merge.df$ext_gene,
                              levels = rownames(clust.mtx)[genes.order])
  merge.df$samples.ord = factor(x = merge.df$sample,
                                levels = colnames(clust.mtx)[sample.order])
  
  print(
    ggplot2::ggplot(merge.df, aes(samples.ord,genes.ord))+
      geom_tile(aes(fill=log2(tpm+0.5)))+
      scale_fill_gradient(low="purple",high = "orange")+
      ggtitle(paste0("Top ",topNgens," DE genes"))+
      theme_bw()+
      theme(
        axis.title= element_blank(),
        axis.text.x=element_text(angle=-90)))
  
}

generate.volcano = function(sleuthObject, qvalueCutoff= 0.01, pvalueCutoff=10**-10, betaCutoff=2, topNgenes = 25){
  sleuthObject <- dplyr::filter(sleuthObject, qval <= qvalueCutoff)
  
  select.genes = sleuthObject[1:topNgenes,]
  
  print(
    ggplot()+
      geom_point(data = sleuthObject,
                 aes(b, -1*log10(pval),label = ext_gene),
                 col="grey")+
      geom_point(data = sleuthObject[(abs(sleuthObject$b) >= betaCutoff)&
                                            (sleuthObject$pval <= pvalueCutoff),],
                 aes(b, -1*log10(pval)), col = "red")+
      geom_vline(xintercept = c(-1*betaCutoff,betaCutoff), linetype =2)+
      geom_hline(yintercept = -1*log10(pvalueCutoff), linetype =2)+
      labs(x = "Beta", y = "-log10(p-value)")+ 
      geom_text_repel(data = select.genes,    
                      box.padding = unit(0.35, "lines"),
                      point.padding = unit(0.5, "lines"),
                      aes(b, -1*log10(pval),label = ext_gene))+ 
      theme_bw())
  
}

# Functions 