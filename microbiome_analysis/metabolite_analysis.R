library("limma")
library("sva")
library("GGally")
setwd("/Users/Aoji.Xie/Desktop/microbiome_analysis/clinical")

df_clinical <-  read.csv("HelsinkiBlood_Master_OKtoUse_data.csv", header = TRUE)
df_clinical <- subset(df_clinical, df_clinical$OK_to_use == "yes")  # filter the OK_to_use samples

dim(df_clinical)


sample_info <- data.frame(df_clinical$X, df_clinical$gender, df_clinical$Group, df_clinical$age_at_stool_collection, df_clinical$tobacco_100_in_life,df_clinical$Wexner_total,df_clinical$Rome_III_IBS_criteria_fulfilled,df_clinical$LED, df_clinical$BMI, df_clinical$GDS_15, df_clinical$meds_SSRI, df_clinical$meds_Tricyclics, df_clinical$history_appendectomy)

colnames(sample_info) <- c("Subject_ID", "Sex","Diagnosis","Age","tobacco_100_in_life","Wexner_total","Rome_III_IBS_criteria_fulfilled","LED","BMI","GDS","meds_SSRI","meds_Tricyclics", "history_appendectomy")


sample_info$meds <- sample_info$meds_SSRI + sample_info$meds_Tricyclics
## write the sample information data
sample_info <- na.omit(sample_info)
dim(sample_info)
write.table(sample_info, "sample_info_raw.txt",row.names = F, sep="\t",quote = F)

sample_info_pd  subset(sample_info, grepl("P", sample_info$Diagnosis))
sample_info_all_pd <- na.omit(sample_info_pd)
dim(sample_info_all_pd )
write.table(sample_info_all_pd, "sample_info_pd.txt", row.names = F, sep = "\t", quote = F)

#################################processed metabolite ##################

metabolite_info <- data.frame(df_clinical$X, df_clinical$acetic_acid_normalized, df_clinical$propionic_acid_normalized, df_clinical$isobutyric_acid_normalized, df_clinical$butyric_acid_normalized, df_clinical$isovaleric_acid_normalized, df_clinical$valeric_acid_normalized)

colnames(metabolite_info) <- c("Subject_ID", "acetic_acid_normalized","propionic_acid_normalized","isobutyric_acid_normalized","butyric_acid_normalized","isovaleric_acid_normalized","valeric_acid_normalized")
metabolite_info <- as.data.frame(t(metabolite_info))
#class(metabolite_info[1,])
colnames(metabolite_info) <- as.character(unlist(metabolite_info[1,]))
metabolite_info <- metabolite_info[-1,]
#head(metabolite_info)

indx <- sapply(metabolite_info, is.character)
metabolite_info[indx] <- lapply(metabolite_info[indx], function(x) as.numeric(as.character(x)))


metabolite_order_all <- metabolite_info[, c(as.character(sample_info$Subject_ID))]
metabolite_order_all_pd <- metabolite_info[, c(as.character(sample_info_all_pd$Subject_ID))]



##################################################linear model #######################################


design_metabolite_diagnosis_all <- model.matrix(~ Diagnosis + Age + Sex + tobacco_100_in_life + BMI, data = sample_info)
dim(design_metabolite_diagnosis_all)

# data matrix
data_metabolite_diagnosis_all <- metabolite_order_all

fit <- lmFit(data_metabolite_diagnosis_all,design= design_metabolite_diagnosis_all,method="robust")
fit <- eBayes(fit)


metabolite_results_diagnosis_all = cbind(rownames(data_metabolite_diagnosis_all),topTable(fit,coef=2,adjust.method="BH",sort="none",n=Inf))
metabolite_results_diagnosis_all_fdr0.05 = subset(metabolite_results_diagnosis_all, metabolite_results_diagnosis_all$adj.P.Val <= 0.05)
metabolite_results_diagnosis_all[,2:6]

write.table(metabolite_results_diagnosis_all, "metabolite/metabolite_results_diagnosis.txt",sep = "\t", quote = F, row.names = F)
write.table(metabolite_results_diagnosis_all_fdr0.05, "metabolite/metabolite_results_diagnosis_fdr0.05.txt",sep = "\t", quote = F,row.names = F)



sample_info_all_pd <- subset(sample_info, sample_info$Diagnosis == "P")
design_metabolite_GDS_all_pd = model.matrix(~ GDS + Age + Sex + tobacco_100_in_life + BMI , data = sample_info_all_pd)
dim(design_metabolite_GDS_all_pd)
data_metabolite_GDS_all_pd = metabolite_order_all_pd
dim(data_metabolite_GDS_all_pd)
fit <- lmFit(data_metabolite_GDS_all_pd,design = design_metabolite_GDS_all_pd,method="robust")
fit <- eBayes(fit)
metabolite_results_GDS_all_pd = cbind(rownames(data_metabolite_GDS_all_pd),topTable(fit,coef=2,adjust.method="BH",sort="none",n=Inf))
metabolite_results_GDS_all_pd_fdr0.05 = subset(metabolite_results_GDS_all_pd, metabolite_results_GDS_all_pd$adj.P.Val <= 0.05)
metabolite_results_GDS_all_pd[,2:6]
write.table(metabolite_results_GDS_all_pd, "metabolite/metabolite_results_GDS_pd.txt" , sep = "\t", quote = F, row.names = F)
write.table(metabolite_results_GDS_all_pd_fdr0.05, "metabolite/metabolite_results_GDS_pd_fdr0.05.txt" , sep = "\t", quote = F, row.names = F)