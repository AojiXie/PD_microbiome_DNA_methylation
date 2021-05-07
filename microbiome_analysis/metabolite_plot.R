


library("ggplot2")
setwd("/Users/Aoji.Xie/Desktop/microbiome_analysis/clinical/metabolite")

diagnosis_metabolite <- read.table("metabolite_results_diagnosis.txt", header = T, sep = "\t")



diagnosis_metabolite$value <- (-1) * sign(diagnosis_metabolite$logFC) * log10(diagnosis_metabolite$P.Value)
diagnosis_metabolite_order <- diagnosis_metabolite[order(diagnosis_metabolite[,c("value")]),]
diagnosis_metabolite_order$sign <- as.character( (-1) * sign(diagnosis_metabolite$logFC))

diagnosis_metabolite_order[,1] <- factor(diagnosis_metabolite_order[,1], levels=c(as.character(diagnosis_metabolite_order[, 1])))

library(stringr)
infor <- data.frame(str_split_fixed(diagnosis_metabolite_order[, 1], "\\_", 3))
mylabels <- paste(infor[, 1], infor[, 2], sep=" ")
mylabels[1] <- "butyrate"
mylabels[2] <- "propionate"
mylabels[3] <- "acetate"
mylabels[4] <- "isobutyrate"
mylabels[5] <- "valerate"
mylabels[6] <- "isovalerate"


fdr_sig <-  log10(1/6 * 0.05)









pdf("diagnosis_metabolite.pdf", width = 5, height = 3)

p10 <- ggplot(diagnosis_metabolite_order, aes(x = rownames.data_metabolite_diagnosis_all., y = value, fill=sign)) + geom_bar(stat="identity", color="black") +
  scale_y_continuous(name = "SLP") + geom_hline(yintercept= c(0), linetype=c( "solid"), color = c("black")) +
  scale_x_discrete(name = " ",  labels = mylabels) +
  theme(axis.line = element_line(colour = "black"), panel.background = element_blank(), plot.title = element_text(size = 20, family="Helvetica", face = "plain"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 20),
        axis.title = element_text(),
        axis.text.x=element_text(size = 20, colour = "black", family="Helvetica", face = "plain"),
        axis.text.y=element_text(size =20, colour = "black", family="Helvetica", face = "plain"),
        legend.position = "none") +  #scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
  scale_fill_brewer(palette = "Accent") + geom_hline(yintercept= fdr_sig  , linetype="dashed", 
                                                     color = "red", size=1) + labs(title="PD vs Control") 
p10 + coord_flip()
dev.off()





GDS_metabolite = read.table("metabolite_results_GDS_pd.txt", header = T, sep = "\t")



GDS_metabolite$value <- (-1) * sign(GDS_metabolite$logFC) * log10(GDS_metabolite$P.Value)
GDS_metabolite_order <- GDS_metabolite[order(GDS_metabolite[,c("value")]),]
GDS_metabolite_order$sign <- as.character( (-1) * sign(GDS_metabolite$logFC))

GDS_metabolite_order[,1] <- factor(GDS_metabolite_order[,1], levels=c(as.character(GDS_metabolite_order[, 1])))

library(stringr)
infor <- data.frame(str_split_fixed(GDS_metabolite_order[, 1], "\\_", 3))
mylabels <- paste(infor[, 1], infor[, 2], sep=" ")

mylabels[1] <- "butyrate"
mylabels[2] <-  "isobutyrate"
mylabels[3] <-  "propionate" 
mylabels[4] <- "acetate"
mylabels[5] <- "isovalerate"
mylabels[6] <- "valerate"



fdr_sig =  log10(0.016)





pdf("GDS_metabolite.pdf", width = 5, height = 3)

p10 <- ggplot(GDS_metabolite_order, aes(x = rownames.data_metabolite_GDS_all_pd., y = value, fill=sign)) + geom_bar(stat="identity", color="black") +
  scale_y_continuous(name = "SLP") + geom_hline(yintercept= c(0), linetype=c( "solid"), color = c("black")) +
  scale_x_discrete(name = " ",  labels = mylabels) +
  theme(axis.line = element_line(colour = "black"), panel.background = element_blank(), plot.title = element_text(size = 20, family="Helvetica", face = "plain"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 20),
        axis.title = element_text(),
        axis.text.x=element_text(size = 20, colour = "black", family="Helvetica", face = "plain"),
        axis.text.y=element_text(size = 20, colour = "black", family="Helvetica", face = "plain"),
        legend.position = "none") +  #scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
  scale_fill_brewer(palette = "Accent") + geom_hline(yintercept= fdr_sig  , linetype="dashed", 
                                                     color = "red", size=1) + labs(title="Correlation with GDS") 
       
p10 + coord_flip()
dev.off()




p_data <- read.table("/Users/Aoji.Xie/Desktop/microbiome_analysis/clinical/sample_info_raw.txt", header = T, sep = "\t")

butyric <-  read.table("/Users/Aoji.Xie/Desktop/microbiome_analysis/clinical/butyric_acid_normalized_df.txt", header = T, sep = "\t")

p_data <- merge(p_data, butyric, by = "Subject_ID")

p_data_PD <- subset(p_data, p_data$Diagnosis.x == "P")

pdf("butyric_cor_GDS.pdf", height = 4, width = 3)
p<-ggplot(p_data_PD, aes(x= GDS, y=  log(butyric_acid_normalized))) + geom_point(color = "gray", size = 3) + geom_smooth(method=lm, se=FALSE,color = "deepskyblue") +
  scale_y_continuous(name = "log(butyrate)")+
  theme(axis.line = element_line(colour = "black"), panel.background = element_blank(), plot.title = element_text(size = 25, family="Helvetica", face = "plain"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 20),
        axis.title = element_text(),
        axis.text.x=element_text(size = 20, colour = "black", family="Helvetica", face = "plain"),
        axis.text.y=element_text(size = 20, colour = "black", family="Helvetica", face = "plain"),
        legend.position = "none")
p
dev.off()
