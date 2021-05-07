## plot Fig 5B


CD_result <- read.table("/Users/Aoji.Xie/Desktop/LDSC_analysis/bed_file/CD_2019_.tsv", header = T, sep = "\t")
#IBD_result = read.table("/Users/Aoji.Xie/Desktop/LDSC_analysis/bed_file/IBD_2019_.tsv", header = T, sep = "\t")
MDD_result <- read.table("/Users/Aoji.Xie/Desktop/LDSC_analysis/bed_file/MDD_2019_.tsv", header = T, sep = "\t", stringsAsFactors = FALSE)
AD_result <- read.table("/Users/Aoji.Xie/Desktop/LDSC_analysis/bed_file/AD_2019_.tsv", header = T, sep = "\t", stringsAsFactors = FALSE)
BIP_result <- read.table("/Users/Aoji.Xie/Desktop/LDSC_analysis/bed_file/BIP_2019_.tsv", header = T, sep = "\t", stringsAsFactors = FALSE)
SCZ_result <- read.table("/Users/Aoji.Xie/Desktop/LDSC_analysis/bed_file/SCZ_2019_.tsv", header = T, sep = "\t", stringsAsFactors = FALSE)
UC_result <- read.table("/Users/Aoji.Xie/Desktop/LDSC_analysis/bed_file/UC_2019_.tsv", header = T, sep = "\t", stringsAsFactors = FALSE)
RA_result <- read.table("/Users/Aoji.Xie/Desktop/LDSC_analysis/bed_file/ra_2019_.tsv", header = T, sep = "\t", stringsAsFactors = FALSE)
CD_result  <- na.omit(CD_result)
MDD_result <- na.omit(MDD_result)

CD_result <- subset(CD_result, CD_result$Description == "Innate Immune System" |CD_result$Description == "Signaling by GPCR" | CD_result$Description == "Adaptive Immune System" | CD_result$Description == "Neutrophil degranulation"| CD_result$Description =="Metabolism of lipids" | CD_result$Description == "Vesicle-mediated transport")

#IBD_result = subset(IBD_result, IBD_result$Description == "Innate Immune System" |IBD_result$Description == "Signaling by GPCR" | IBD_result$Description == "Adaptive Immune System" | IBD_result$Description == "Neutrophil degranulation"| IBD_result$Description =="Metabolism of lipids" | IBD_result$Description == "Vesicle-mediated transport")

MDD_result <- subset(MDD_result, MDD_result$Description == "Innate Immune System" |MDD_result$Description == "Signaling by GPCR" | MDD_result$Description == "Adaptive Immune System" | MDD_result$Description == "Neutrophil degranulation"| MDD_result$Description =="Metabolism of lipids" | MDD_result$Description == "Vesicle-mediated transport")

BIP_result <- subset(BIP_result, BIP_result$Description == "Innate Immune System" |BIP_result$Description == "Signaling by GPCR" | BIP_result$Description == "Adaptive Immune System" | BIP_result$Description == "Neutrophil degranulation"| BIP_result$Description =="Metabolism of lipids" | BIP_result$Description == "Vesicle-mediated transport")

SCZ_result <- subset(SCZ_result, SCZ_result$Description == "Innate Immune System" |SCZ_result$Description == "Signaling by GPCR" | SCZ_result$Description == "Adaptive Immune System" | SCZ_result$Description == "Neutrophil degranulation"| SCZ_result$Description =="Metabolism of lipids" | SCZ_result$Description == "Vesicle-mediated transport")

AD_result <- subset(AD_result, AD_result$Description == "Innate Immune System" |AD_result$Description == "Signaling by GPCR" | AD_result$Description == "Adaptive Immune System" | AD_result$Description == "Neutrophil degranulation"| AD_result$Description =="Metabolism of lipids" | AD_result$Description == "Vesicle-mediated transport")

UC_result <- subset(UC_result, UC_result$Description == "Innate Immune System" |UC_result$Description == "Signaling by GPCR" | UC_result$Description == "Adaptive Immune System" | UC_result$Description == "Neutrophil degranulation"| UC_result$Description =="Metabolism of lipids" | UC_result$Description == "Vesicle-mediated transport")

RA_result <- subset(RA_result, RA_result$Description == "Innate Immune System" |RA_result$Description == "Signaling by GPCR" | RA_result$Description == "Adaptive Immune System" | RA_result$Description == "Neutrophil degranulation"| RA_result$Description =="Metabolism of lipids" | RA_result$Description == "Vesicle-mediated transport")


CD_result$disease <- "Crohn's disease"
#IBD_result$disease <- "Inflammatory Bowel Disease"
MDD_result$disease <- "Major depression"
BIP_result$disease <- "Bipolar disorder"
SCZ_result$disease <- "Schizophrenia"
UC_result$disease <- "Ulcerative colitis"
AD_result$disease <- "Alzheimer's disease"
RA_result$disease <- "Rheumatoid arthritis"

result_table <- rbind(UC_result,CD_result,RA_result,BIP_result,AD_result,SCZ_result,MDD_result, RA_result)
result_table$logfdr <- -log(result_table$FDR.q.value)

level_order <- c('Ulcerative colitis',  'Crohn\'s disease' , 'Rheumatoid arthritis','Bipolar disorder','Schizophrenia','Alzheimer\'s disease','Major depression' )
#scale_fill_gradientn(colours=c("darkgreen", "green", "greenyellow", "yellow", "red", "darkred"), breaks=seq(0, 1, 0.1), "CM", limit=c(min_cm, max_cm))
library(ggplot2)
mine.heatmap <- ggplot(data = result_table, mapping = aes(x = factor(disease, levels = level_order),
                                                       y = Description,
                                                       fill = logfdr)) +
  geom_tile() +
  xlab(label = "") + 
  ylab(label = "Pathway")+
  scale_fill_gradient2(
                      low = "yellow",
                      mid = "red",
                      high = "darkred",
                      midpoint = 10) +
  theme(axis.line = element_line(colour = "black"), panel.background = element_blank(), plot.title = element_text(size = 5, family="Helvetica", face = "plain"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15),
        axis.title = element_text(size = 15, colour = "black", family="Helvetica", face = "plain"),
        axis.text.x=element_text(size = 15, colour = "black", family="Helvetica", face = "plain", angle = 45,hjust = 1 ),
        axis.text.y=element_text(size =15, colour = "black", family="Helvetica", face = "plain")
        
  ) 

pdf("/Users/Aoji.Xie/Desktop/LDSC_analysis/bed_file/plot/heatmap_1.pdf", height = 5, width = 7)

mine.heatmap
dev.off()



