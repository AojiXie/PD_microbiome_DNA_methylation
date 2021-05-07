## process the result of find nearest file
## constract Butyrate related methyDNA information for LDSC analysis
## we used https://github.com/bulik/ldsc  for this analysis

library(dplyr)
library(data.table)

setwd("/Users/Aoji.Xie/Desktop/LDSC_analysis/common_snp")
all_chr <- list()

for (m in 1:22) {
  name <- paste("chr", m, sep = "_")
  all_chr[[m]] <- read.csv(name, sep = "\t", header = T)
  
}

v1 <- all_chr[[1]]

ncol(v1)
ncol(all_chr[[1]])

for(n in 2:22){
  v1 = rbind(v1, all_chr[[n]])
}



prob_snp <- subset(v1, v1$SNP != "no")


write.table(prob_snp,"prob_snp.txt", sep = "\t", row.names = FALSE, quote = FALSE)






## all snp information
all_snp <- read.table("/Users/Aoji.Xie/Desktop/LDSC/w_hm3.snplist", header = T, sep = "\t")
head(all_snp)

prob_snp <- merge(prob_snp,all_snp, by = "SNP")
head(prob_snp)

##select SNP within 5000bp of each probe
prob_snp <- subset(prob_snp,prob_snp$distance < 5000)

##pick the smallest p value if a probe is related with multple SNPs 

table <- as.data.table(prob_snp)
setkey(table, SNP)
result <- table[, list(pval = min(P.Value), probe=Row.names[which.min(P.Value)], B=B[which.min(P.Value)], A1 = A1[which.min(P.Value)], A2 = A2[which.min(P.Value)]) ,by = SNP]
# convert result into normal data.frame 
result <- as.data.frame(result)



write.table(result,"test_file.txt", sep = "\t", row.names = FALSE, quote = FALSE)
