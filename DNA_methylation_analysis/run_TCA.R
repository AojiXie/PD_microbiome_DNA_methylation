## use TCA to estimate epigenome of different cell types 
## this is used in Fig 4 
library("data.table")
library("TCA")
library(doParallel)
detectCores()

X1 <- data.frame(fread(file = "methy_data_array_slide_combat_for_model.txt"),row.names = 1)
W1 <- data.frame(fread(file = "/Users/Aoji.Xie/Desktop/test_TCA/combined_result.txt"), row.names=1)

dim(X1)
head(X1)[1:10]
dim(W1)
head(W1)


X1 <- X1[,65:128]
colnames(X1)
dim(X1)

W1["C0005",] = (W1["C0005",] + W1["C0005rep",])/2
W1["C0148",] = (W1["C0148",] + W1["C0148rep",])/2
W1["P0072",] = (W1["P0072",] + W1["P0072rep",])/2
W1["P0042",] = (W1["P0042",] + W1["P0042rep",])/2


dim(W1)

row.names.remove <- c("C0005rep","C0148rep","P0072rep","P0042rep")
W1 = W1[!(row.names(W1) %in% row.names.remove),]





W1 = W1[colnames(X1),]




tca.mdl <- tca(X = X1, W = W1, 
               parallel = TRUE, num_cores = 12, log_file = "GSE42861.tca.log")

mdl.tca.sub_test = tcasub(tca.mdl, features = rownames(X1), log_file = "GSE42861.tcasub.log")
Z_hat_test <- tensor(X = X1, mdl.tca.sub_test, log_file = "GSE42861.tensor.log")


meth_cd8 <- data.frame(Z_hat_test[[1]])
meth_cd4 <- data.frame(Z_hat_test[[2]])
meth_NK <- data.frame(Z_hat_test[[3]])
meth_Bcell <- data.frame(Z_hat_test[[4]])
meth_Mono <- data.frame(Z_hat_test[[5]])
meth_Neu <- data.frame(Z_hat_test[[6]])
head(meth_cd8)

write.table(meth_cd8, "TCA_celltypes/type_cd8_pd.txt",sep = "\t", quote = FALSE)
write.table(meth_cd4, "TCA_celltypes/type_cd4_pd.txt",sep = "\t", quote = FALSE)
write.table(meth_NK, "TCA_celltypes/type_NK_pd.txt",sep = "\t", quote = FALSE)
write.table(meth_Bcell, "TCA_celltypes/type_Bcell_pd.txt",sep = "\t", quote = FALSE)
write.table(meth_Mono, "TCA_celltypes/type_Mono_pd.txt",sep = "\t", quote = FALSE)
write.table(meth_Neu, "TCA_celltypes/type_Neu_pd.txt",sep = "\t", quote = FALSE)
