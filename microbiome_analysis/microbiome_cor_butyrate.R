## correlation between butyrate and microbiome
setwd("/Users/Aoji.Xie/Desktop/microbiome_analysis/microbiome/microbiome_butyrate_cor")


microbiom <- read.table("/Users/Aoji.Xie/Desktop/microbiome_analysis/microbiome/mentagenomeSeq/order/microbiom_order_GenusLevel.txt", header=T, sep = "\t")
head(microbiom)
dim(microbiom)
microbiom <- microbiom[,7:ncol(microbiom)]

rownames(microbiom) <-  microbiom$Genus
microbiom$Genus = NULL

################################IN pd

butyric_acid <- read.table("butyric_acid_normalized_df.txt", sep = "\t", header = T)
butyric_acid_PD <- subset(butyric_acid, butyric_acid$Diagnosis == "P")
rownames(butyric_acid_PD) <- butyric_acid_PD$Subject_ID
butyric_acid_PD <- na.omit(butyric_acid_PD)


col = intersect(butyric_acid_PD$Subject_ID, colnames(microbiom))

microbiom_intrested <- microbiom[, col]
butyric_acid_PD <- butyric_acid_PD[col,]
microbiom_intrested <- microbiom_intrested[,rownames(butyric_acid_PD)]


estimates <- numeric(nrow(microbiom_intrested))
pvalues <- numeric(nrow(microbiom_intrested))


for (i in 1:nrow(microbiom_intrested)){
  test <- cor.test(butyric_acid_PD$butyric_acid_normalized, t(microbiom_intrested[i,1:ncol(microbiom_intrested)]))
  estimates[i] <- test$estimate
  pvalues[i] <- test$p.value
}

cor_otu_butyric_acid <- data.frame(microbiom = rownames(microbiom_intrested), P.value = "", Estimate = "")

cor_otu_butyric_acid$P.value <- pvalues
cor_otu_butyric_acid$Estimate <- estimates
cor_otu_butyric_acid$fdr <- p.adjust(cor_otu_butyric_acid$P.value, method="fdr")

tst <-  subset(cor_otu_butyric_acid, cor_otu_butyric_acid$fdr <= 0.05)

taxa <- read.table("/Users/Aoji.Xie/Desktop/microbiome_analysis/microbiome/mentagenomeSeq/order/microbiom_order_GenusLevel.txt", header = T, sep = "\t")


taxa <-  taxa[,2:7]
colnames(cor_otu_butyric_acid)[1] <- "Genus"

cor_otu_butyric_acid <- merge(taxa, cor_otu_butyric_acid, by = "Genus")


write.table(cor_otu_butyric_acid,"cor_microbiom_butyric_acid_all_microbe_table_vs_butyric_in_PD_Genus.txt",sep = "\t", row.names = F, quote = F)




############ALL GROUP ##################################


butyric_acid <- read.table("butyric_acid_normalized_df.txt", sep = "\t", header = T)



rownames(butyric_acid) <- butyric_acid$Subject_ID
butyric_acid <- na.omit(butyric_acid)
dim(butyric_acid_PD)
col <- intersect(butyric_acid$Subject_ID, colnames(microbiom))

microbiom_intrested <-microbiom[, col]
butyric_acid <- butyric_acid[col,]

microbiom_intrested <- microbiom_intrested[,rownames(butyric_acid)]


estimates <- numeric(nrow(microbiom_intrested))
pvalues <- numeric(nrow(microbiom_intrested))


for (i in 1:nrow(microbiom_intrested)){
  test <- cor.test(butyric_acid$butyric_acid_normalized, t(microbiom_intrested[i,1:ncol(microbiom_intrested)]))
  estimates[i] <- test$estimate
  pvalues[i] <- test$p.value
}

cor_otu_butyric_acid <- data.frame(microbiom = rownames(microbiom_intrested), P.value = "", Estimate = "")

#cor_food_butyric_acid$Food =colnames(food_info)
cor_otu_butyric_acid$P.value <- pvalues
cor_otu_butyric_acid$Estimate <- estimates
cor_otu_butyric_acid$fdr <- p.adjust(cor_otu_butyric_acid$P.value, method="fdr")


colnames(cor_otu_butyric_acid)[1] <-  "Genus"

cor_otu_butyric_acid <-  merge(taxa, cor_otu_butyric_acid, by = "Genus")

write.table(cor_otu_butyric_acid,"cor_microbiom_butyric_acid_all_microbe_table_vs_butyric_in_All_group_Genus.txt",sep = "\t", row.names = F, quote = F)