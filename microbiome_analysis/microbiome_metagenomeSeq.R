library("metagenomeSeq")

setwd("/Users/Aoji.Xie/Desktop/microbiome_analysis/microbiome/mentagenomeSeq_1")

taxa <- read.delim("taxa.txt", sep = "\t", stringsAsFactors=FALSE,row.names=1)
taxa <- taxa[,1:6]

butyrate <- read.table("/Users/Aoji.Xie/Desktop/microbiome_analysis/clinical/butyric_acid_normalized_df.txt", header = T, sep = "\t")
butyrate <- na.omit(butyrate)
dim(butyrate)

sample_info_all_cov <- read.table("sample_info_all_cov.txt", header = T, sep = "\t")
rownames(sample_info_all_cov) = sample_info_all_cov$Subject_ID
colnames(sample_info_all_cov)
sample_info_all_cov = subset(sample_info_all_cov, sample_info_all_cov$Subject_ID %in% butyrate$Subject_ID)
dim(sample_info_all_cov)



sample_info_all_cov_pd <- subset(sample_info_all_cov,sample_info_all_cov$Diagnosis == "P") 

count_table <- read.csv("scheperjans_16S_otu_counts.csv", header = TRUE)
count_table <- count_table[,-1]

count_table <- count_table[,as.character(sample_info_all_cov$Subject_ID)]
count_table_pd <- count_table[,as.character(sample_info_all_cov_pd$Subject_ID)]
#dim(count_table_all)
#######################################  in ALL group ###########################
sample_data_all <- AnnotatedDataFrame(sample_info_all_cov) 
OTU_data <- AnnotatedDataFrame(taxa) 
sample_data_pd <- AnnotatedDataFrame(sample_info_all_cov_pd) 


obj_all <- newMRexperiment(count_table,phenoData=sample_data_all,featureData=OTU_data)
obj_pd <-  newMRexperiment(count_table_pd,phenoData=sample_data_pd,featureData=OTU_data)

# 5% of samples
rareFeatures <- which(rowSums(MRcounts(obj_all)> 0) < 6)
#rareFeatures
obj_all <- obj_all[-rareFeatures,]
#dim(obj_all)

#obj_all = filterData(obj_all, present = 10, depth = 1)

count_table_after_filter <- MRcounts(obj_all)
class(count_table_after_filter)



Age <- pData(obj_all)$Age
Sex <- pData(obj_all)$Sex
Diagnosis <- pData(obj_all)$Diagnosis
BMI <-  pData(obj_all)$BMI
tobacco_100_in_life <- pData(obj_all)$tobacco_100_in_life


##################################
run_metagenom_seq_1 <- function(exp_obj, maxit){
  p = cumNormStat(exp_obj)
  exp_obj = cumNorm(exp_obj, p = p)
  
  normFactor = normFactors(exp_obj)
  normFactor = log2(normFactor/median(normFactor) + 1)
  mod = model.matrix( ~  Diagnosis + BMI + Age +  Sex +tobacco_100_in_life + normFactor)
  
  settings = zigControl(maxit = 10, verbose = TRUE)
  fit = fitZig(obj = exp_obj, mod = mod , control = settings)
  
  mod_r = model.matrix( ~  BMI + Age +  Sex +tobacco_100_in_life +  normFactor)
  fit_r = fitZig(obj = exp_obj, mod = mod_r , control = settings)
  
  #residuals_ <- residuals.MArrayLM(fit_r,  count_table_after_filter)
  
  
  sig = calculateEffectiveSamples(fit)
  effective = as.data.frame(sig)
  colnames(effective)[1] = "ave"
  effective_otu = subset(effective, effective$ave >  mean(effective$ave))
  
  
  
  topTable = MRcoefs(fit,coef = 2 , number = 1000)
  topTable_filtered = merge(effective_otu, topTable, by = 0)
  rownames(topTable_filtered) = topTable_filtered$Row.names
  
  #residuals_filtered =  merge(effective_otu, residuals_, by = 0)
  #rownames(residuals_filtered) = residuals_filtered$Row.names
  
  topTable_fdr0.05 = subset(topTable, topTable$adjPvalues <= 0.05)
  topTable_effective_fdr0.05 =  merge(effective_otu, topTable_fdr0.05, by = 0)
  rownames(topTable_effective_fdr0.05) = topTable_effective_fdr0.05$Row.names
  return(list(topTable,topTable_filtered, topTable_fdr0.05,topTable_effective_fdr0.05,fit_r))
  
}

output_result_otu = function (results, outputfile){
  
  result_1 <-  merge(taxa, results[1],by = 0)
  outfile_1 = paste(outputfile, ".txt", sep="")
  write.table(result_1, outfile_1, row.names = F, sep = "\t", quote = F)
  
  result_2 <-  merge(taxa, results[2],by = 0)
  outfile_2 <- paste(outputfile, "_filtered.txt", sep="")
  write.table(result_2, outfile_2, row.names = F, sep = "\t", quote = F)
  
  
  result_3 <-  merge(taxa, results[3],by = 0)
  outfile_3 <- paste(outputfile, "_fdr0.05.txt", sep="")
  write.table(result_3, outfile_3, row.names = F, sep = "\t", quote = F)
  
  
  
  result_4 <-  merge(taxa, results[4],by = 0)
  outfile_4 <- paste(outputfile, "_filtered_fdr0.05.txt", sep="")
  write.table(result_4, outfile_4, row.names = F, sep = "\t", quote = F)
  
}



output_result_levels <- function (results, outputfile,level_file, level){
  
  result_1 <- as.data.frame(results[1])
  result_1[4] <- rownames(result_1)
  colnames(result_1)[4] <- level
  result_1 <- result_1[,c(4,1:3)]
  result_1 <- merge(level_file, result_1, by = level)
  
  outfile_1 <- paste(outputfile, ".txt", sep="")
  write.table(result_1, outfile_1, row.names = F, sep = "\t", quote = F)
  
  
  
  result_2 <- as.data.frame(results[2])
  result_2$Row.names = NULL
  result_2[5] <- rownames(result_2)
  colnames(result_2)[5] <- level
  
  result_2 =result_2[,c(5,1:4)]
  result_2 = merge(level_file, result_2, by = level)
  
  outfile_2 <- paste(outputfile, "_filtered.txt", sep="")
  write.table(result_2, outfile_2, row.names = F, sep = "\t", quote = F)
  
  
  result_3 <- as.data.frame(results[3])
  result_3[4] <- rownames(result_3)
  colnames(result_3)[4] <- level
  result_3 <- result_3[,c(4,1:3)]
  result_3 <- merge(level_file, result_3, by = level)
  
  outfile_3 <- paste(outputfile, "_fdr0.05.txt", sep="")
  write.table(result_3, outfile_3, row.names = F, sep = "\t", quote = F)
  
  
  result_4 <- as.data.frame(results[4])
  result_4$Row.names = NULL
  result_4[5] <- rownames(result_4)
  colnames(result_4)[5] <- level
  result_4 <- result_4[,c(5,1:4)]
  result_4 <-  merge(level_file, result_4, by = level)
  
  outfile_4 <-  paste(outputfile, "_filtered_fdr0.05.txt", sep="")
  write.table(result_4, outfile_4, row.names = F, sep = "\t", quote = F)
  
}
obj_Genus <-aggTax(obj_all, lvl = "Genus") 
obj_Family <- aggTax(obj_all, lvl = "Family") 
obj_Order <-  aggTax(obj_all, lvl = "Order") 
obj_Phylum <-  aggTax(obj_all, lvl = "Phylum") 
obj_Class <- aggTax(obj_all, lvl = "Class") 

Class_level <- read.table("/Users/Aoji.Xie/Desktop/microbiome_analysis/microbiome/mentagenomeSeq/order/Class_level_file.txt", sep = "\t", header = T)
Genus_level <- read.table("/Users/Aoji.Xie/Desktop/microbiome_analysis/microbiome/mentagenomeSeq/order/Genus_level_file.txt", sep = "\t", header = T)
Family_level <- read.table("/Users/Aoji.Xie/Desktop/microbiome_analysis/microbiome/mentagenomeSeq/order/Family_level_file.txt", sep = "\t", header = T)
Order_level <- read.table("/Users/Aoji.Xie/Desktop/microbiome_analysis/microbiome/mentagenomeSeq/order/Order_level_file.txt", sep = "\t", header = T)
Phylum_level <- read.table("/Users/Aoji.Xie/Desktop/microbiome_analysis/microbiome/mentagenomeSeq/order/Phylum_level_file.txt", sep = "\t", header = T)


topTable_OTU_1 <- run_metagenom_seq_1(obj_all, 10)
output_result_otu(topTable_OTU_1, outputfile = "result/microbiome_diagnosis") 
Genus_1 <-  run_metagenom_seq_1(obj_Genus,10)
output_result_levels(Genus_1, outputfile = "result/microbiome_diagnosis_Genus", level = "Genus", level_file = Genus_level)
Family_1 <- run_metagenom_seq_1(obj_Family,10)
output_result_levels(Family_1, outputfile = "result/microbiome_diagnosis_Family", level = "Family", level_file = Family_level)
Order_1 <- run_metagenom_seq_1(obj_Order,10)
output_result_levels(Order_1, outputfile = "result/microbiome_diagnosis_Order", level = "Order", level_file = Order_level)
Class_1 <- run_metagenom_seq_1(obj_Class,10)
output_result_levels(Class_1, outputfile = "result/microbiome_diagnosis_Class", level = "Class", level_file = Class_level)
Phylum_1 <- run_metagenom_seq_1(obj_Phylum,10)
output_result_levels(Phylum_1, outputfile = "result/microbiome_diagnosis_Phylum", level = "Phylum", level_file = Phylum_level)

### save residuals
fit_r <- topTable_OTU_1[5]
fit_r[[1]]@fit
residuals_ <- residuals.MArrayLM(fit_r[[1]]@fit, count_table_after_filter)
residual_table <- as.data.frame(residuals_)

topTable_OTU_1[2]
residual_table_filtered <- merge(residual_table,topTable_OTU_1[2], by = 0)


residual_table_filtered <- residual_table_filtered[,1:129]

write.table (residual_table_filtered, "result/microbiome_diagnosis_residual_table_all_s_marker.txt", sep = "\t", row.names = F, quote = F)

###############################In PD#######################

## 5% of sample
rareFeatures <- which(rowSums(MRcounts(obj_pd)> 0) < 3)
#rareFeatures
obj_pd = obj_pd[-rareFeatures,]
dim(obj_pd)



count_table_after_filter = MRcounts(obj_pd)
class(count_table_after_filter)


GDS <- pData(obj_pd)$GDS
Age <- pData(obj_pd)$Age
Sex <- pData(obj_pd)$Sex
BMI <- pData(obj_pd)$BMI
tobacco_100_in_life <- pData(obj_pd)$tobacco_100_in_life



run_metagenom_seq_2 <- function(exp_obj, maxit){
  p <- cumNormStat(exp_obj)
  exp_obj <- cumNorm(exp_obj, p = p)
  
  normFactor <- normFactors(exp_obj)
  normFactor <- log2(normFactor/median(normFactor) + 1)
  mod <- model.matrix( ~  GDS + BMI + Age +  Sex +tobacco_100_in_life + normFactor)
  
  settings <- zigControl(maxit = 10, verbose = TRUE)
  fit <- fitZig(obj = exp_obj, mod = mod , control = settings)
  
  mod_r <- model.matrix( ~  BMI + Age +  Sex +tobacco_100_in_life  + 
                          normFactor)
  fit_r <- fitZig(obj = exp_obj, mod = mod_r , control = settings)
  
  #residuals_ <- residuals.MArrayLM(fit_r,  count_table_after_filter)
  
  
  sig <- calculateEffectiveSamples(fit)
  effective <- as.data.frame(sig)
  colnames(effective)[1] <- "ave"
  effective_otu <- subset(effective, effective$ave >  mean(effective$ave))
  
  
  
  topTable <- MRcoefs(fit,coef = 2 , number = 1000)
  topTable_filtered <- merge(effective_otu, topTable, by = 0)
  rownames(topTable_filtered) <- topTable_filtered$Row.names
  
  #residuals_filtered =  merge(effective_otu, residuals_, by = 0)
  #rownames(residuals_filtered) = residuals_filtered$Row.names
  
  topTable_fdr0.05 <- subset(topTable, topTable$adjPvalues <= 0.05)
  topTable_effective_fdr0.05 <- merge(effective_otu, topTable_fdr0.05, by = 0)
  rownames(topTable_effective_fdr0.05) <- topTable_effective_fdr0.05$Row.names
  return(list(topTable,topTable_filtered, topTable_fdr0.05,topTable_effective_fdr0.05,fit_r))
  
}

obj_Genus_pd <-  aggTax(obj_pd, lvl = "Genus") 
obj_Family_pd <- aggTax(obj_pd, lvl = "Family") 
obj_Order_pd <- aggTax(obj_pd, lvl = "Order") 
obj_Phylum_pd <- aggTax(obj_pd, lvl = "Phylum") 
obj_Class_pd <- aggTax(obj_pd, lvl = "Class") 
dim(obj_pd)

topTable_OTU_1 <- run_metagenom_seq_2(obj_pd, 10)
output_result_otu(topTable_OTU_1, outputfile = "result/microbiome_GDS") 
Genus_1 <- run_metagenom_seq_2(obj_Genus_pd,10)
output_result_levels(Genus_1, outputfile = "result/microbiome_GDS_Genus", level = "Genus", level_file = Genus_level)
Family_1 <- run_metagenom_seq_2(obj_Family_pd,10)
output_result_levels(Family_1, outputfile = "result/microbiome_GDS_Family", level = "Family", level_file = Family_level)
Order_1 <- run_metagenom_seq_2(obj_Order_pd,10)
output_result_levels(Order_1, outputfile = "result/microbiome_GDS_Order", level = "Order", level_file = Order_level)
Class_1 <- run_metagenom_seq_2(obj_Class_pd,10)
output_result_levels(Class_1, outputfile = "result/microbiome_GDS_Class", level = "Class", level_file = Class_level)
Phylum_1 <- run_metagenom_seq_2(obj_Phylum_pd,10)
output_result_levels(Phylum_1, outputfile = "result/microbiome_GDS_Phylum", level = "Phylum", level_file = Phylum_level)
