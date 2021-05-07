## plot the resuts for Fig 5A
setwd("/Users/Aoji.Xie/Desktop/LDSC_analysis/common_snp")

library(data.table)
##########################pd##################

my_data = read.table("test.sumstats", header = T, sep = "\t") 
colnames(my_data) <- c("SNP",  "Z_ba","A1_ba", "A2_ba",  "N_ba")


pd_snp = read.table("LD_PD_2019/trait.sumstats", header = T, sep = "\t")
colnames(pd_snp) <- c("SNP",  "A1_pd", "A2_pd","Z_pd", "N_pd")
all_ld_score = data.frame(fread("/Users/Aoji.Xie/Desktop/LDSC/all_LD_score"))



keys <- c("SNP")

all_ld_score_table <- data.table(all_ld_score, key=keys)
my_data_table <-  data.table(my_data, key=keys)
pd_snp_table <-  data.table(pd_snp , key=keys)



overlap <- all_ld_score_table[my_data_table, nomatch=0]
overlap <- overlap[pd_snp_table, nomatch=0]

overlap <- data.frame(overlap)
dim(overlap)  ## 298666 

write.table(overlap, "overlap.txt", quote=FALSE, sep="\t", row.names=FALSE)

overlap <- read.table("overlap.txt", header=T, sep="\t")
overlap <- na.omit(overlap)


model <- lm(Z_pd~ Z_ba+ L2, data= overlap)

summary(model)


summary_pd = data.frame(disease = "Parkinson's disease",
                        p_val = 1.16e-06,
                        class = "Neurodegenerative")

'''
Call:
lm(formula = Z_pd ~ Z_ba + L2, data = overlap)

Residuals:
    Min      1Q  Median      3Q     Max 
-1.3788 -0.5190 -0.1412  0.3720  9.9963 

Coefficients:
             Estimate Std. Error t value Pr(>|t|)    
(Intercept) 0.8014857  0.0026329 304.411  < 2e-16 ***
Z_ba        0.0069337  0.0014262   4.862 1.16e-06 ***
L2          0.0021688  0.0000631  34.371  < 2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.6703 on 296862 degrees of freedom
Multiple R-squared:  0.00403,	Adjusted R-squared:  0.004024 
F-statistic: 600.7 on 2 and 296862 DF,  p-value: < 2.2e-16

'''

########################SCZ#########################################
scz_snp = read.table("LD-scz/trait.sumstats", header = T, sep = "\t")
colnames(scz_snp) <- c("SNP",  "A1_scz", "A2_scz", "Z_scz","N_scz")
head(scz_snp)
all_ld_score = data.frame(fread("/Users/Aoji.Xie/Desktop/LDSC/all_LD_score"))



keys <- c("SNP")

all_ld_score_table <- data.table(all_ld_score, key=keys)
my_data_table <-  data.table(my_data, key=keys)
scz_snp_table <-  data.table(scz_snp , key=keys)



overlap <- all_ld_score_table[my_data_table, nomatch=0]
overlap <- overlap[scz_snp_table, nomatch=0]

overlap <- data.frame(overlap)

write.table(overlap, "overlap.txt", quote=FALSE, sep="\t", row.names=FALSE)

overlap <- read.table("overlap.txt", header=T, sep="\t")
overlap <- na.omit(overlap) 


model <- lm(Z_scz~ Z_ba+ L2, data= overlap)

summary(model)
'''
Call:
lm(formula = Z_pd ~ Z_ba + L2, data = overlap)

Residuals:
    Min      1Q  Median      3Q     Max 
-1.3788 -0.5190 -0.1412  0.3720  9.9963 

Coefficients:
             Estimate Std. Error t value Pr(>|t|)    
(Intercept) 0.8014857  0.0026329 304.411  < 2e-16 ***
Z_ba        0.0069337  0.0014262   4.862 1.16e-06 ***
L2          0.0021688  0.0000631  34.371  < 2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.6703 on 296862 degrees of freedom
Multiple R-squared:  0.00403,	Adjusted R-squared:  0.004024 
F-statistic: 600.7 on 2 and 296862 DF,  p-value: < 2.2e-16

'''
cor.test(overlap$Z_ba, overlap$Z_scz)

summary_scz = data.frame(disease = "Schizophrenia",
                        p_val = 0.144 ,
                        class = "Psychiatric")


###########################depression####################
depression_snp = read.table("LD-depression/mdd.sumstats", header = T, sep = "\t")
colnames(depression_snp) <- c("SNP",  "A1_depression", "A2_depression", "Z_depression","N_depression")
head(depression_snp)
all_ld_score = data.frame(fread("/Users/Aoji.Xie/Desktop/LDSC/all_LD_score"))



keys <- c("SNP")

all_ld_score_table <- data.table(all_ld_score, key=keys)
my_data_table <-  data.table(my_data, key=keys)
depression_snp_table <-  data.table(depression_snp , key=keys)



overlap <- all_ld_score_table[my_data_table, nomatch=0]
overlap <- overlap[depression_snp_table, nomatch=0]
overlap <- data.frame(overlap)
write.table(overlap, "overlap.txt", quote=FALSE, sep="\t", row.names=FALSE)

overlap <- read.table("overlap.txt", header=T, sep="\t")
overlap <- na.omit(overlap)


model <- lm(Z_depression~ Z_ba+ L2, data= overlap)

summary(model)

'''

Call:
lm(formula = Z_depression ~ Z_ba + L2, data = overlap)

Residuals:
    Min      1Q  Median      3Q     Max 
-1.0804 -0.4990 -0.1282  0.3684  4.0374 

Coefficients:
              Estimate Std. Error t value Pr(>|t|)    
(Intercept)  8.055e-01  2.988e-03 269.587   <2e-16 ***
Z_ba        -6.644e-04  1.601e-03  -0.415    0.678    
L2           1.078e-03  6.705e-05  16.080   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.6289 on 212212 degrees of freedom
Multiple R-squared:  0.001217,	Adjusted R-squared:  0.001208 
F-statistic: 129.3 on 2 and 212212 DF,  p-value: < 2.2e-16
'''

summary_mdd = data.frame(disease = "Major depression",
                        p_val = 0.678,
                        class = "Psychiatric")

 
#####################BIP#############################

bip_snp = read.table("LD-bip/trait.sumstats", header = T, sep = "\t")
colnames(bip_snp) <- c("SNP",  "A1_bip", "A2_bip", "Z_bip","N_bip")
head(bip_snp)
all_ld_score = data.frame(fread("/Users/Aoji.Xie/Desktop/LDSC/all_LD_score"))



keys <- c("SNP")

all_ld_score_table <- data.table(all_ld_score, key=keys)
my_data_table <-  data.table(my_data, key=keys)
bip_snp_table <-  data.table(bip_snp , key=keys)

overlap <- all_ld_score_table[my_data_table, nomatch=0]
overlap <- overlap[bip_snp_table, nomatch=0]
overlap <- data.frame(overlap)


write.table(overlap, "overlap.txt", quote=FALSE, sep="\t", row.names=FALSE)

overlap <- read.table("overlap.txt", header=T, sep="\t")
overlap <- na.omit(overlap)  

model <- lm(Z_bip~ Z_ba+ L2, data= overlap)

summary(model)
'''
Call:
lm(formula = Z_bip ~ Z_ba + L2, data = overlap)

Residuals:
    Min      1Q  Median      3Q     Max 
-1.3187 -0.5295 -0.1384  0.3841  5.0310 

Coefficients:
             Estimate Std. Error t value Pr(>|t|)    
(Intercept) 8.107e-01  3.476e-03 233.237  < 2e-16 ***
Z_ba        6.123e-03  1.863e-03   3.286  0.00102 ** 
L2          2.360e-03  7.701e-05  30.648  < 2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.6753 on 183656 degrees of freedom
Multiple R-squared:  0.005159,	Adjusted R-squared:  0.005148 
F-statistic: 476.2 on 2 and 183656 DF,  p-value: < 2.2e-16

'''

summary_bip = data.frame(disease = "Bipolar disorder",
                         p_val =  0.00102,
                         class = "Psychiatric")




################################UC###############################
UC_snp = read.table("LD-UC/trait.sumstats", header = T, sep = "\t")
colnames(UC_snp) <- c("SNP",  "Z_UC","A1_UC", "A2_UC", "N_UC")
head(UC_snp)
all_ld_score = data.frame(fread("/Users/Aoji.Xie/Desktop/LDSC/all_LD_score"))



keys <- c("SNP")

all_ld_score_table <- data.table(all_ld_score, key=keys)
my_data_table <-  data.table(my_data, key=keys)
UC_snp_table <-  data.table(UC_snp , key=keys)



overlap <- all_ld_score_table[my_data_table, nomatch=0]
overlap <- overlap[UC_snp_table, nomatch=0]
overlap <- data.frame(overlap)

write.table(overlap, "overlap.txt", quote=FALSE, sep="\t", row.names=FALSE)

overlap <- read.table("overlap.txt", header=T, sep="\t")
overlap <- na.omit(overlap)  #206056


model <- lm(Z_UC~ Z_ba+ L2, data= overlap)

summary(model)

'''
Call:
lm(formula = Z_UC ~ Z_ba + L2, data = overlap)

Residuals:
    Min      1Q  Median      3Q     Max 
-1.3219 -0.5632 -0.1575  0.3884 12.7340 

Coefficients:
             Estimate Std. Error t value Pr(>|t|)    
(Intercept) 8.600e-01  3.517e-03 244.520  < 2e-16 ***
Z_ba        9.409e-03  1.915e-03   4.913 8.97e-07 ***
L2          2.230e-03  8.237e-05  27.077  < 2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.7452 on 204975 degrees of freedom
Multiple R-squared:  0.003673,	Adjusted R-squared:  0.003664 
F-statistic: 377.9 on 2 and 204975 DF,  p-value: < 2.2e-16
'''

summary_UC = data.frame(disease = "Ulcerative colitis",
                         p_val =   8.97e-07,
                         class = "Gastrointestinal")


#################################AD################################
AD_snp = read.table("IGAP_summary_statistics/trait.sumstats", header = T, sep = "\t")
colnames(AD_snp) <- c("SNP", "A1_AD", "A2_AD", "Z_AD", "N_AD")
head(AD_snp)
all_ld_score = data.frame(fread("/Users/Aoji.Xie/Desktop/LDSC/all_LD_score"))



keys <- c("SNP")

all_ld_score_table <- data.table(all_ld_score, key=keys)
my_data_table <-  data.table(my_data, key=keys)
AD_snp_table <-  data.table(AD_snp , key=keys)



overlap <- all_ld_score_table[my_data_table, nomatch=0]
overlap <- overlap[AD_snp_table, nomatch=0]
overlap <- data.frame(overlap)


write.table(overlap, "overlap.txt", quote=FALSE, sep="\t", row.names=FALSE)

overlap <- read.table("overlap.txt", header=T, sep="\t")
overlap <- na.omit(overlap)

model <- lm(Z_AD~ Z_ba+ L2, data= overlap)

summary(model)
'''
Call:
lm(formula = Z_AD ~ Z_ba + L2, data = overlap)

Residuals:
    Min      1Q  Median      3Q     Max 
-3.5039 -0.7492 -0.2115  0.5657  7.0897 

Coefficients:
             Estimate Std. Error t value Pr(>|t|)    
(Intercept)  3.735951   0.093814  39.823   <2e-16 ***
Z_ba        -0.082923   0.048940  -1.694   0.0906 .  
L2           0.001850   0.001634   1.133   0.2577    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 1.256 on 832 degrees of freedom
Multiple R-squared:  0.004634,	Adjusted R-squared:  0.002242 
F-statistic: 1.937 on 2 and 832 DF,  p-value: 0.1448
'''



summary_AD = data.frame(disease = "Alzheimer's disease",
                        p_val = 0.0906,
                        class = "Neurodegenerative")


###########################CD###########################

CD_snp = read.table("LD_CD/trait.sumstats", header = T, sep = "\t")
colnames(CD_snp) <- c("SNP", "A1_CD", "A2_CD", "Z_CD", "N_CD")
all_ld_score = data.frame(fread("/Users/Aoji.Xie/Desktop/LDSC/all_LD_score"))



keys <- c("SNP")

all_ld_score_table <- data.table(all_ld_score, key=keys)
my_data_table <-  data.table(my_data, key=keys)
CD_snp_table <-  data.table(CD_snp , key=keys)



overlap <- all_ld_score_table[my_data_table, nomatch=0]
overlap <- overlap[CD_snp_table, nomatch=0]
overlap <- data.frame(overlap)

write.table(overlap, "overlap.txt", quote=FALSE, sep="\t", row.names=FALSE)

overlap <- read.table("overlap.txt", header=T, sep="\t")
overlap <- na.omit(overlap)


model <- lm(Z_CD ~ Z_ba+ L2, data= overlap)

summary(model)
'''
Call:
lm(formula = Z_AD ~ Z_ba + L2, data = overlap)

Residuals:
    Min      1Q  Median      3Q     Max 
-3.5039 -0.7492 -0.2115  0.5657  7.0897 

Coefficients:
             Estimate Std. Error t value Pr(>|t|)    
(Intercept)  3.735951   0.093814  39.823   <2e-16 ***
Z_ba        -0.082923   0.048940  -1.694   0.0906 .  
L2           0.001850   0.001634   1.133   0.2577    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 1.256 on 832 degrees of freedom
Multiple R-squared:  0.004634,	Adjusted R-squared:  0.002242 
F-statistic: 1.937 on 2 and 832 DF,  p-value: 0.1448
'''


summary_CD = data.frame(disease = "Crohn's disease",
                        p_val =  2.83e-05,
                        class = "Gastrointestinal")



###################################ra#############################

ra_snp = read.table("LD_RA/trait_ra.sumstats", header = T, sep = "\t")
colnames(ra_snp) <- c("SNP",  "A1_ra", "A2_ra","Z_ra", "N_ra")
all_ld_score = data.frame(fread("/Users/Aoji.Xie/Desktop/LDSC/all_LD_score"))



keys <- c("SNP")

all_ld_score_table <- data.table(all_ld_score, key=keys)
my_data_table <-  data.table(my_data, key=keys)
ra_snp_table <-  data.table(ra_snp , key=keys)



overlap <- all_ld_score_table[my_data_table, nomatch=0]
overlap <- overlap[ra_snp_table, nomatch=0]
overlap <- data.frame(overlap)


write.table(overlap, "overlap.txt", quote=FALSE, sep="\t", row.names=FALSE)

overlap <- read.table("overlap.txt", header=T, sep="\t")
overlap <- na.omit(overlap)


model <- lm(Z_ra~ Z_ba+ L2, data= overlap)

summary(model)
'''
Call:
lm(formula = Z_AD ~ Z_ba + L2, data = overlap)

Residuals:
    Min      1Q  Median      3Q     Max 
-3.5039 -0.7492 -0.2115  0.5657  7.0897 

Coefficients:
             Estimate Std. Error t value Pr(>|t|)    
(Intercept)  3.735951   0.093814  39.823   <2e-16 ***
Z_ba        -0.082923   0.048940  -1.694   0.0906 .  
L2           0.001850   0.001634   1.133   0.2577    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 1.256 on 832 degrees of freedom
Multiple R-squared:  0.004634,	Adjusted R-squared:  0.002242 
F-statistic: 1.937 on 2 and 832 DF,  p-value: 0.1448
'''


summary_ra = data.frame(disease = "Rheumatoid arthritis",
                        p_val = 0.000367,
                        class = "Autoimmune")


############################
summary_table = data.frame(disease = character(),
                           p_val = numeric(),
                           class = character())

summary_table = rbind(summary_AD,summary_bip, summary_CD, summary_ra, summary_mdd, summary_pd, summary_scz, summary_UC)

summary_table$logP = log10(summary_table$p_val) * -1

summary_table <- summary_table[order(summary_table[,c("logP")]),]
summary_table[,1] <- factor(summary_table[,1], levels=c(as.character(summary_table[, 1])))
library(ggplot2)
pdf("plot.pdf", width = 7.5, height = 4)
p10 <- ggplot(summary_table, aes(x = disease, y = logP, fill=class)) + geom_bar(stat="identity") +
    #scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
    scale_y_continuous(name = "-logP")  +
    scale_x_discrete() +
    
    theme(axis.line = element_line(colour = "black"), panel.background = element_blank(), plot.title = element_text(size = 20, family="Helvetica", face = "plain"), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text = element_text(size = 15),
          axis.title = element_text(size = 15, colour = "black", family="Helvetica", face = "plain"),
          axis.text.x=element_text(size = 15, colour = "black", family="Helvetica", face = "plain"),
          axis.text.y=element_text(size =15, colour = "black", family="Helvetica", face = "plain"),
          legend.title = element_blank()) +  
    scale_fill_brewer(palette = "Accent") +
    geom_hline(yintercept=1.3, linetype="dashed", 
               color = "red", size=1.5) +
     labs(title="LD score regression") 
p10 + coord_flip()
dev.off()


