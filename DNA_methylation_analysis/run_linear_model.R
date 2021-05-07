## run linear model on different cell type epigenome. related with Fig 4 

library("limma")

p_data <- read.table("/Users/Aoji.Xie/Desktop/DNA_Methylation/test_M/pheno_data.txt", sep = "\t", header = T)
butyric_acid <- read.table("/Users/Aoji.Xie/Desktop/DNA_Methylation/butyric_acid_normalized_df.txt" , header = T, sep = "\t")

p_data_butyric <- merge(p_data, butyric_acid, by ="Subject_ID")
p_data_butyric <- na.omit(p_data_butyric)

p_data_butyric_PD <- subset(p_data_butyric,p_data_butyric$Sample_Group == "P")
p_data_butyric_control <- subset(p_data_butyric,p_data_butyric$Sample_Group == "C")



probe.features <- data.frame(fread("/Users/Aoji.Xie/Desktop/DNA_Methylation/test_M/probe_features.txt"),row.names = 1)








run_linear_model <-  function(mod_full, meth_data, design_res, file_path){
  
  
  fit <- lmFit(meth_data,design = mod_full ,method="robust")
  fit <- eBayes(fit)
  fit
  results <- cbind(rownames(meth_data),toptable(fit,coef=2,adjust.method="BH",sort="none",n=Inf))
  results_fdr0.05 <- subset(results, results$adj.P.Val <= 0.05)
  #dim(methylation_results_diagnosis_comB_fdr0.05) ## 6608
  
  results <-  merge(results,probe.features, by = 0)
  results_fdr0.05 <- merge(results_fdr0.05, probe.features, by = 0)
  
  file_results <- paste(file_path,".txt",sep = "")
  file_results_fdr <- paste(file_path,"_fdr0.05.txt", sep = "")
  
  write.table(results, file = file_results, sep = "\t", quote = F, row.names = F)
  write.table(results_fdr0.05, file = file_results_fdr, sep = "\t", quote = F, row.names = F)
  
  #####################################bed file
  if(nrow(results_fdr0.05) > 0) {
    bed <- results_fdr0.05[, c("CHR","MAPINFO","MAPINFO","Row.names")]
    colnames(bed) <- c("CHR","start","end","probe_ID")
    bed$end <- bed$start+1
    
    bed$CHR <- paste("chr", bed$CHR, sep = "")
    
    
    file_bed <-  paste(file_path,"_fdr0.05.bed", sep = "")
    write.table(bed, file=file_bed, quote=F, sep="\t", row.names=F, col.names=F)
  }
  
  
  fit_r = lmFit(meth_data,design= design_res,method="robust")
  residuals_result <- residuals.MArrayLM(fit_r, meth_data)
  
  residuals_result = as.data.frame(residuals_result)
  residuals_result$Probe_ID = rownames(residuals_result)
  residuals_result = residuals_result[,c(ncol(residuals_result),1:ncol(residuals_result) - 1)]
  
  file_res = paste(file_path,"_resduals.txt", sep = "")
  write.table(residuals_result , file = file_res,sep = "\t", quote = F,row.names = F)
  
  
  
}
###################output from TCA.R###########################
meth_cd8 <- data.frame(fread("TCA_celltypes/type_cd8_pd.txt"), row.names = 1)
meth_cd4 <- data.frame(fread("TCA_celltypes/type_cd4_pd.txt"), row.names = 1)
meth_NK <- data.frame(fread("TCA_celltypes/type_NK_pd.txt"), row.names = 1)
meth_Bcell <- data.frame(fread("TCA_celltypes/type_Bcell_pd.txt"), row.names = 1)
meth_Mono <- data.frame(fread("TCA_celltypes/type_Mono_pd.txt"), row.names = 1)
meth_Neu <- data.frame(fread("TCA_celltypes/type_Neu_pd.txt"),row.names = 1)







design_methylation_diagnosis <- model.matrix(~ butyric_acid_normalized + Age + Sex + tobacco_100_in_life + BMI , data = p_data_butyric_PD)

mod_r = model.matrix(~ Age + Sex + tobacco_100_in_life + BMI , data =p_data_butyric_PD)



run_linear_model( mod_full = design_methylation_diagnosis, meth_data = data_methylation_diagnosis_cd8, mod_r,
                  file_path = "statistical_analysis_celltypes/methylation_butyric_cd8" )
run_linear_model( mod_full = design_methylation_diagnosis, meth_data = data_methylation_diagnosis_cd4, mod_r,
                  file_path = "statistical_analysis_celltypes/methylation_butyric_cd4" )
run_linear_model( mod_full = design_methylation_diagnosis, meth_data = data_methylation_diagnosis_NK, mod_r,
                  file_path = "statistical_analysis_celltypes/methylation_butyric_NK" )
run_linear_model( mod_full = design_methylation_diagnosis, meth_data = data_methylation_diagnosis_Bcell, mod_r,
                  file_path = "statistical_analysis_celltypes/methylation_butyric_Bcell" )
run_linear_model( mod_full = design_methylation_diagnosis, meth_data = data_methylation_diagnosis_Mono, mod_r,
                  file_path = "statistical_analysis_celltypes/methylation_butyric_Mono_result" )
run_linear_model( mod_full = design_methylation_diagnosis, meth_data = data_methylation_diagnosis_Neu, mod_r,
                  file_path = "statistical_analysis_celltypes/methylation_butyric_Neu" )
