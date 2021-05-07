## this is for pathway analysis and the result is shown on Fig 5 B
## This file is to find the intersect genome location between butyrate and PD and on of other 
## diseases(BIP, SCZ, MDD, AD, RA, UC, CD)
## we used https://www.gsea-msigdb.org/gsea/index.jsp for pathway analysis (reactome)

setwd("/Users/Aoji.Xie/Desktop/LDSC_analysis")
library(data.table)
library(stringr)
bip = data.frame(fread("/Users/Aoji.Xie/Desktop/LDSC/LD-bip/pgc.bip.full.2012-04.txt"))

colnames(bip)[1] <- 'SNP'
colnames(bip)[2] <- "CHR"
colnames(bip)[3] <- "BP"
colnames(bip)[4] <- "A1"
colnames(bip)[5] <- "A2"
colnames(bip)[6] <- "or"
colnames(bip)[7] <- "se"
colnames(bip)[8] <- "pval"

bip$q = p.adjust(bip$pval)
bip_sig = subset(bip, bip$pval < 0.05)
bip_bed = bip_sig[,c("CHR","BP","BP","SNP")]
bip_bed$CHR = paste("chr", bip_bed$CHR, sep = "")
  



mdd <- data.frame(fread("/Users/Aoji.Xie/Desktop/LDSC/LD-depression/pgc.mdd.full.2012-04.txt"))

colnames(mdd)[1] <- 'SNP'
colnames(mdd)[2] <- "CHR"
colnames(mdd)[3] <- "BP"
colnames(mdd)[4] <- "A1"
colnames(mdd)[5] <- "A2"
colnames(mdd)[6] <- "or"
colnames(mdd)[7] <- "se"
colnames(mdd)[8] <- "pval"


mdd_sig <- subset(mdd, mdd$pval < 0.05)
mdd_bed <- mdd_sig[,c("CHR","BP","BP","SNP")]
mdd_bed$CHR = paste("chr", mdd_bed$CHR, sep = "")

########### ####################

SCZ <- data.frame(fread("/Users/Aoji.Xie/Desktop/LDSC/LD-scz/pgc.scz.full.2012-04.txt"))
head(SCZ)
colnames(SCZ)[1] <- "SNP"
colnames(SCZ)[2] <- "CHR"
colnames(SCZ)[3] <- "BP"
colnames(SCZ)[4] <- "A1"
colnames(SCZ)[5] <- "A2"
colnames(SCZ)[6] <- "or"
colnames(SCZ)[7] <- "se"
colnames(SCZ)[8] <- "pval"

SCZ_sig <- subset(SCZ, SCZ$pval < 0.05)
SCZ_bed <- SCZ_sig[,c("CHR","BP","BP","SNP")]
SCZ_bed$CHR <- paste("chr", SCZ_bed$CHR, sep = "")


AD <- data.frame(fread("/Users/Aoji.Xie/Desktop/LDSC/IGAP_summary_statistics/IGAP_stage_1_2_combined.txt"))
AD_sig <- subset(AD, AD$Pvalue < 0.05)
AD_bed <- AD_sig[, c("Chromosome","Position","Position","MarkerName")]
colnames(AD_bed)[4] <- "SNP"
colnames(AD_bed)[1] <- "CHR"
AD_bed$CHR = paste("chr", AD_bed$CHR, sep = "")
#########################################################

pd <- data.frame(fread("/Users/Aoji.Xie/Desktop/LDSC/LD_PD_2019/nallsEtAl2019_excluding23andMe_allVariants.tab"))
pd_sig <- subset(pd, pd$p <0.05)
pd_1 <- data.frame(fread("/Users/Aoji.Xie/Desktop/LDSC/LD_PD_2019/pd_trait_2019.txt"))
annotation <- data.frame(fread("/Users/Aoji.Xie/Desktop/LDSC/Chang_2017_PDWBS-5.0/gwas_anno-5.0/all_snp_info-5.0.txt"))
annotation$SNP <- paste(annotation$scaffold, annotation$position, sep = ":")
annotation <- annotation[,c("assay.name", "SNP")]



pd_sig <- merge(annotation, pd_sig, by ="SNP")
df <- data.frame(str_split_fixed(pd_sig$SNP, ":", 2))
colnames(df)[1] <- "CHR"
colnames(df)[2] <- "BP"

pd_sig <- cbind(df, pd_sig)
pd_bed <- pd_sig[,c("CHR","BP","BP","assay.name")]
colnames(pd_bed)[4] <- "SNP"


UC <- data.frame(fread("/Users/Aoji.Xie/Desktop/LDSC/LD-UC/UC-trait.txt"))
head(UC)
UC_sig <- subset(UC,UC$pval<0.05)
UC_bed <- UC_sig[,c("CHR","POS","POS","SNP")]
UC_bed$CHR <-  paste("chr", UC_bed$CHR, sep = "")




CD = data.frame(fread("/Users/Aoji.Xie/Desktop/LDSC/LD-IBD/CD/EUR.CD.gwas_info03_filtered.assoc"))
CD_sig = subset(CD, CD$P < 0.05)
CD_bed = CD_sig[,c("CHR","BP","BP","SNP")]
CD_bed$CHR = paste("chr", CD_bed$CHR, sep = "")



ra <- data.frame(fread("/Users/Aoji.Xie/Desktop/LDSC/LD-ra/RA_GWASmeta_European_v2.txt"))
ra_sig <- subset(ra, ra$pval < 0.05)
ra_bed <- ra_sig[,c("Chr","Position.hg19.","Position.hg19.","SNPID")]
colnames(ra_bed) <- c("CHR","BP","BP","SNP")
ra_bed$CHR <- paste("chr", ra_bed$CHR, sep = "")



PD_IBD <- subset(pd_bed, pd_bed$SNP %in% IBD_bed$SNP)
PD_CD <- subset(pd_bed, pd_bed$SNP %in% CD_bed$SNP)
PD_mdd <- subset(pd_bed, pd_bed$SNP %in% mdd_bed$SNP)
PD_SCZ <- subset(pd_bed, pd_bed$SNP %in% SCZ_bed$SNP)
PD_bip <- subset(pd_bed, pd_bed$SNP %in% bip_bed$SNP)
PD_UC <- subset(pd_bed, pd_bed$SNP %in% UC_bed$SNP)
PD_AD <- subset(pd_bed, pd_bed$SNP %in% AD_bed$SNP)
PD_ra <- subset(pd_bed, pd_bed$SNP %in% ra_bed$SNP)





PD_CD_butyric <- subset(PD_CD, PD_CD$SNP%in%butyric_sig$SNP)
write.table(PD_CD_butyric, "bed_file/PD_CD_butyric.bed", quote = F, sep = "\t", row.names = F, col.names = F)


PD_SCZ_butyric <-  subset(PD_SCZ, PD_SCZ$SNP%in%butyric_sig$SNP)
write.table(PD_SCZ_butyric, "bed_file/PD_SCZ_butyric.bed", quote = F, sep = "\t", row.names = F, col.names = F)

PD_mdd_butyric <- subset(PD_mdd, PD_mdd$SNP%in%butyric_sig$SNP)
write.table(PD_mdd_butyric, "bed_file/PD_mdd_butyric.bed", quote = F, sep = "\t", row.names = F, col.names = F)

PD_bip_butyric <- subset(PD_bip, PD_bip$SNP%in%butyric_sig$SNP)
write.table(PD_bip_butyric, "bed_file/PD_bip_butyric.bed", quote = F, sep = "\t", row.names = F, col.names = F)

PD_UC_butyric <- subset(PD_UC, PD_UC$SNP%in%butyric_sig$SNP)
write.table(PD_UC_butyric, "bed_file/PD_UC_butyric.bed", quote = F, sep = "\t", row.names = F, col.names = F)

PD_AD_butyric <- subset(PD_AD, PD_AD$SNP%in%butyric_sig$SNP)
write.table(PD_AD_butyric, "bed_file/PD_AD_butyric.bed", quote = F, sep = "\t", row.names = F, col.names = F)

PD_ra_butyric <- subset(PD_ra, PD_ra$SNP%in%butyric_sig$SNP)
write.table(PD_ra_butyric, "bed_file/PD_ra_butyric.bed", quote = F, sep = "\t", row.names = F, col.names = F)