library(data.table)
library(limma)

setwd("/Users/Aoji.Xie/Desktop/DNA_methylation_analysis")
probe.features = data.frame(fread("/Users/Aoji.Xie/Desktop/DNA_Methylation/test_M/probe_features.txt"),row.names = 1)
head(probe.features)
p_data = read.table("pheno_data.txt", sep = "\t", header = T)

#### use the normalized data
methy_data_array_slide_combat_for_model = data.frame(fread("methy_data_array_slide_combat_for_model.txt"), row.names = 1)


run_linear_model <-  function(mod_full, meth_data, design_res, file_path){
  
  
  fit <- lmFit(meth_data,design = mod_full ,method="robust")
  fit <- eBayes(fit)

  results <- cbind(rownames(meth_data),topTable(fit,coef=2,adjust.method="BH",sort="none",n=Inf))
  results_fdr0.05 <- subset(results, results$adj.P.Val <= 0.05)

  
  results <- merge(results,probe.features, by = 0)
  results_fdr0.05 <- merge(results_fdr0.05, probe.features, by = 0)
  
  file_results <- paste(file_path,".txt",sep = "")
  file_results_fdr <- paste(file_path,"_fdr0.05.txt", sep = "")
  
  write.table(results, file = file_results, sep = "\t", quote = F, row.names = F)
  write.table(results_fdr0.05, file = file_results_fdr, sep = "\t", quote = F, row.names = F)
  
  #####################################bed file
  
  if(nrow(results_fdr0.05) > 0) {
    ## create bed file for GREAT analysis
    bed <- results_fdr0.05[, c("CHR","MAPINFO","MAPINFO","Row.names")]
    colnames(bed) <- c("CHR","start","end","probe_ID")
    bed$end <- bed$start+1
    
    bed$CHR <- paste("chr", bed$CHR, sep = "")
    
    
    file_bed <- paste(file_path,"_fdr0.05.bed", sep = "")
    write.table(bed, file=file_bed, quote=F, sep="\t", row.names=F, col.names=F)
  }
  
  
  fit_r <- lmFit(meth_data,design= design_res,method="robust")
  residuals_result <- residuals.MArrayLM(fit_r, meth_data)
  
  residuals_result <- as.data.frame(residuals_result)
  residuals_result$Probe_ID = rownames(residuals_result)
  residuals_result <- residuals_result[,c(ncol(residuals_result),1:ncol(residuals_result) - 1)]
  
  file_res <- paste(file_path,"_resduals.txt",sep = "")
  write.table(residuals_result , file = file_res,sep = "\t", quote = F,row.names = F)
  
  
}




butyric_acid <- read.table("butyric_acid_normalized_df.txt" , header = T, sep = "\t")

p_data_butyric <- merge(p_data, butyric_acid, by ="Subject_ID")
p_data_butyric <- na.omit(p_data_butyric)

p_data_butyric_PD <- subset(p_data_butyric,p_data_butyric$Sample_Group == "P")
p_data_butyric_control <- subset(p_data_butyric,p_data_butyric$Sample_Group == "C")

design_methylation_butyric <- model.matrix(~ butyric_acid_normalized  + Age + Sex + tobacco_100_in_life  + BMI + Mono + Neu + 
                                            CD4T + CD8T + Bcell, data = p_data_butyric_PD)

data_methylation_butyric <-  methy_data_array_slide_combat_for_model[,as.character(p_data_butyric_PD$Subject_ID)]

mod_r_butyric_PD <- model.matrix(~ Age + Sex + tobacco_100_in_life  + BMI  + Mono + Neu + 
                                 CD4T + CD8T + Bcell, data = p_data_butyric_PD)
### run linear model in PD
run_linear_model( mod_full = design_methylation_butyric, meth_data = data_methylation_butyric, mod_r_butyric_PD,
                  file_path = "butyric/methylation_results_PD_butyric" )


design_methylation_butyric <- model.matrix(~ butyric_acid_normalized  + Age + Sex + tobacco_100_in_life + BMI  + Mono + Neu + 
                                            CD4T + CD8T + Bcell, data = p_data_butyric_control)

data_methylation_butyric <- methy_data_array_slide_combat_for_model[,as.character(p_data_butyric_control$Subject_ID)]
dim(data_methylation_butyric)
mod_r_butyric_control <- model.matrix(~ Age + Sex + tobacco_100_in_life + BMI  + Mono + Neu + 
                                      CD4T + CD8T + Bcell, data = p_data_butyric_PD)
## run linear model in Control
run_linear_model( mod_full = design_methylation_butyric, meth_data = data_methylation_butyric, mod_r_butyric_control,
                  file_path = "butyric/methylation_results_control_butyric" )
