## plot the dot plot in Fig 2
library(ggplot2)

setwd("/Users/Aoji.Xie/Desktop/microbiome_analysis/microbiome/microbiome_butyrate_cor")
microbiom <- read.table("/Users/Aoji.Xie/Desktop/microbiome_analysis/microbiome/mentagenomeSeq/order/microbiom_order_GenusLevel.txt", header=T, sep = "\t")
head(microbiom)
dim(microbiom)
microbiom <- microbiom[,7:ncol(microbiom)]

rownames(microbiom) <-  microbiom$Genus
microbiom$Genus <- NULL



butyric_acid <- read.table("butyric_acid_normalized_df.txt", sep = "\t", header = T)

butyric_acid_PD <- subset(butyric_acid, butyric_acid$Diagnosis == "P")
rownames(butyric_acid_PD) <- butyric_acid_PD$Subject_ID
butyric_acid_PD = na.omit(butyric_acid_PD)



col <- intersect(butyric_acid_PD$Subject_ID, colnames(microbiom))

microbiom_intrested <- microbiom[, col]
butyric_acid_PD <- butyric_acid_PD[col,]

microbiom_intrested <- microbiom_intrested[,rownames(butyric_acid_PD)]


estimates <- numeric(nrow(microbiom_intrested))
pvalues <- numeric(nrow(microbiom_intrested))



m <- t(microbiom_intrested[,1:ncol(microbiom_intrested)])
m <- data.frame(m)
M = data.frame(m$Roseburia)



plot_table = cbind(butyric_acid_PD, M)
pdf("Roseburia_butyric_plot.pdf", width = 3, height = 4)
p<-ggplot(plot_table, aes(x= log10(butyric_acid_normalized), y= log(m.Roseburia))) + geom_point(color = "gray", size = 3) + geom_smooth(method=lm, se=FALSE,color = "#E69F00") +
  scale_y_continuous(name = "log10(Roseburia Counts)")+
  scale_x_continuous(name = "log10(Butyrate)") + 
  theme(axis.line = element_line(colour = "black"), panel.background = element_blank(), plot.title = element_text(size = 20, family="Helvetica", face = "plain"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 20),
        axis.title = element_text(),
        axis.text.x=element_text(size = 20, colour = "black", family="Helvetica", face = "plain"),
        axis.text.y=element_text(size = 20, colour = "black", family="Helvetica", face = "plain"),
        legend.position = "none") + labs(title="Roseburia
") 
p
dev.off()


M = data.frame(m$Romboutsia)



plot_table = cbind(butyric_acid_PD, M)
pdf("Romboutsia_butyric_plot.pdf", width = 3, height = 4)
p<-ggplot(plot_table, aes(x= log10(butyric_acid_normalized), y= log(m.Rombousia)))) + geom_point(color = "gray", size = 3) + geom_smooth(method=lm, se=FALSE,color = "#E69F00") +
  scale_y_continuous(name = "log10(Romboutsia Counts)")+
  scale_x_continuous(name = "log10(Butyrate)") + 
  theme(axis.line = element_line(colour = "black"), panel.background = element_blank(), plot.title = element_text(size = 20, family="Helvetica", face = "plain"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 20),
        axis.title = element_text(),
        axis.text.x=element_text(size = 20, colour = "black", family="Helvetica", face = "plain"),
        axis.text.y=element_text(size = 20, colour = "black", family="Helvetica", face = "plain"),
        legend.position = "none") + labs(title="Romboutsia
") 
p
dev.off()