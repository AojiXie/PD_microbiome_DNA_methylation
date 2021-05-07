## plot tree for Fig S1 C
library(data.tree)
library(treemap)

#######################PLOT IN FIG S1#######################################################



#write.table(microbiom, "/Users/Aoji.Xie/Desktop/PD_microbiome/microbiome_all_cov_analysis/all_cov/microbiom_diagnosis_results_tree.txt", quote=F)

library(ggtree)
library(ggplot2)
##https://www.bioconductor.org/packages/release/bioc/vignettes/ggtree/inst/doc/treeAnnotation.html


correlation = read.table("/Users/Aoji.Xie/Desktop/microbiome_analysis/microbiome/microbiome_butyrate_cor/cor_microbiom_butyric_acid_all_microbe_table_vs_butyric_in_All_group_Genus.txt", header = T, sep = "\t", stringsAsFactors = F)
head(correlation)
dim(correlation)
correlation = subset(correlation, correlation$Genus %in% microbiom_diagnosis_results$Genus)

correlation$pathString <- paste("microbiom", correlation$Kingdom, correlation$Phylum, correlation$Class, correlation$Order,correlation$Family, correlation$Genus, sep = "/")
microbiom <- as.Node(correlation)
microbiom_string <- as.phylo.Node(microbiom)

#write.table(microbiom, "/Users/Aoji.Xie/Desktop/PD_microbiome/microbiome_all_cov_analysis/all_cov/microbiom_diagnosis_results_tree_Genus.txt", quote=F)


p <- ggtree(microbiom_string)  

#p <- p %<+% dd + geom_nodelab(aes(label = node))
# generate some random values for each tip label in the data
d1 <- data.frame(id= correlation$Genus, val= (-1) * sign(correlation$Estimate) * log10(correlation$fdr))
# Make a second plot with the original, naming the new plot "dot", 
# using the data you just created, with a point geom.
p2 <- facet_plot(p, panel = "dot", data=d1, geom=geom_point, aes(x=val),  color='red')




d <- data.frame(x = -log10(0.05), .panel = 'dot')
p3 <- p2 + geom_vline(data = d, aes(xintercept= x),color = "black",linetype='dashed',size=0.5)
d2 <- data.frame(x = log10(0.05), .panel = 'dot')
p4 <- p3 + geom_vline(data = d2, aes(xintercept= x),color = "black",linetype='dashed',size=0.5)
p4 + theme_tree2()


pdf("/Users/Aoji.Xie/Desktop/microbiome_analysis/microbiome/microbiome_butyrate_cor/microbiom_cor_butyrate_in_All_group_Genus.pdf", height = 4, width = 4)

p4+ theme_tree2()







microbiom_GDS_results  <- read.table("/Users/Aoji.Xie/Desktop/microbiome_analysis/microbiome/mentagenomeSeq_1/result/microbiome_GDS_Genus_filtered.txt", header=T, sep="\t", stringsAsFactors = F)
head(microbiom_GDS_results)
dim(microbiom_GDS_results)

microbiom_GDS_results$pathString <- paste("microbiom", microbiom_GDS_results$Kingdom, microbiom_GDS_results$Phylum, microbiom_GDS_results$Class, microbiom_GDS_results$Order,microbiom_GDS_results$Family, microbiom_GDS_results$Genus, sep = "/")
microbiom <- as.Node(microbiom_GDS_results)
microbiom_string <- as.phylo.Node(microbiom)

#write.table(microbiom, "/Users/Aoji.Xie/Desktop/microbiome_analysis/microbiome/mentagenomeSeq/microbiom_GDS_results_tree.txt", quote=F)

library(ggtree)
library(ggplot2)
##https://www.bioconductor.org/packages/release/bioc/vignettes/ggtree/inst/doc/treeAnnotation.html



p <- ggtree(microbiom_string)  

#p <- p %<+% dd + geom_nodelab(aes(label = node))
# generate some random values for each tip label in the data
d1 <- data.frame(id= microbiom_GDS_results$Genus, val= (-1) * sign(microbiom_GDS_results$GDS) * log10(microbiom_GDS_results$adjPvalues))
# Make a second plot with the original, naming the new plot "dot", 
# using the data you just created, with a point geom.
p2 <- facet_plot(p, panel = "dot", data=d1, geom=geom_point, aes(x=val),  color='#56B4E9')




d <- data.frame(x = -log10(0.05), .panel = 'dot')
p3 <- p2 + geom_vline(data = d, aes(xintercept= x),color = "black",linetype='dashed',size=0.5)
d2 <- data.frame(x = log10(0.05), .panel = 'dot')
p4 <- p3 + geom_vline(data = d2, aes(xintercept= x),color = "black",linetype='dashed',size=0.5)

pdf("/Users/Aoji.Xie/Desktop/microbiome_analysis/microbiome/mentagenomeSeq_1/microbiom_GDS_results_tree_test_Genus_modified.pdf", height = 4, width = 4)

p4+ theme_tree2()

dev.off()