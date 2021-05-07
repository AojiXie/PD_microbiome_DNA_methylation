

### This is for the plot in Fig 3A
library("ggplot2")
library("data.table")
library(magrittr) # needs to be run every time you start R and want to use %>%
library(dplyr)
library(xlsx)
library("latex2exp")



library(xlsx)
GWAS_PD_gene = read.xlsx("/Users/Aoji.Xie/Desktop/DNA_Methylation/Hi-C_data_process/GWAS_PD_gene.xlsx", sheetIndex  = 1)
GWAS_PD_gene = GWAS_PD_gene[,1]
GWAS_depression_gene = read.xlsx("/Users/Aoji.Xie/Desktop/DNA_Methylation/depression_GWAS.xlsx", sheetIndex = 1, startRow = 3)

GWAS_depression_gene = GWAS_depression_gene$geneName
length(GWAS_depression_gene)
is.na(GWAS_depression_gene)


manhattanPlot <- function(fit, qThresh = 0.05, labelGenes = TRUE) {
  
  
  
  fit <- copy(fit)
  # Mark the fit points which have to be labeled
  geneList <- GWAS_PD_gene
  
  
  
  foo <- function(snps, pvals) {
    snps[ which.min(pvals)]
  }
  
  
  
  foo2 <- function(N, SNP) {
    i <- which.max(N)
    return(SNP[i])
  }
  
  
  
  
  getColorGroup <- function(Chr, Q, S) {
    Chr <- as.numeric(Chr)
    Q <- as.numeric(Q)
    S <- as.numeric(S)
    df <- data.frame(Chr, Q, S)
    apply(df, 1, function(I) {
      I <- as.list(I)
      if (I$Chr %% 2 == 0) {
        # Darker colors
        if (I$Q > qThresh) {
          "black"
        } else {
          #ifelse(I$S > 0, "dark positive", "dark negative")
          "black"
        }
      } else {
        # Brighter colors
        if (I$Q > qThresh) {
          "light grey"
        } else {
          #ifelse(I$S > 0, "light positive", "light negative")
          "light grey"
        }
      }
    }) %>% unlist
  }
  
  
  
  if (labelGenes == TRUE) {
    effects <- 
      fit[ toupper(Gene) %in% geneList & adj.P.Val < qThresh, ] %>% 
      # count how many times each gene is negative of positive
      .[, list(.N, SNP = foo(SNP, P.Value), P = min(P.Value)), by=list(Gene, sign(logFC))] %>%
      # select the sign which has more representation
      .[, list(SNP=foo2(N, SNP)), Gene] %>%
      .[, Label := TRUE] %>% 
      merge(fit, ., by=c("Gene", "SNP"), all.x=TRUE)
    effects[is.na(Label), Label := FALSE]
    
    
    
    effects <- effects[Chr != "chrM"]    
  } else {
    effects <- fit[Chr != "chrM"] %>%
      .[, Label := FALSE] %>%
      .[, Gene := ""]
  }
  
  
  
  
  # Format the input data
  effects <- effects %>%
    .[, Diff := logFC ] %>%
    .[, SNP := strsplit(ID, split="_") %>% sapply(., '[[', 2) %>% as.numeric] %>%
    .[, P := P.Value] %>%
    .[, Q := adj.P.Val] %>%
    # Fix chromosome name 
    .[, Chr := gsub("chr", "", Chr) ] %>% 
    # Give chromosomes numbers
    .[, ChrNo :=  Chr %>% gsub("X", "23", .) %>% gsub("Y", "24", .) %>% gsub("M", "25", .) %>%as.numeric]  %>%
    setkey(ChrNo, SNP)
  
  
  
  # Prepare plot data
  effects <- effects[!is.na(Q)]
  pd <- effects %>% 
    # Subsample 
    .[ Q < 0.05 | ID %in% sample(ID, .N*0.3, prob=-log10(P))] %>%  
    # Compute chromosome size
    .[, list(ChrLen = max(SNP)), ChrNo] %>%
    # Calculate cumulative position of each chromosome
    .[, list(ChrNo, Tot = cumsum(ChrLen) - ChrLen)] %>%
    # Add this info to the initial dataset
    merge(effects, ., by=c("ChrNo")) %>%
    # Sort by ChrNo and position 
    .[ order(ChrNo, SNP), ] %>%
    # Compute position for each SNP on x axis
    .[, Loc := SNP + Tot] %>%
    # Get color group
    .[, Color := getColorGroup(ChrNo, Q, Diff) %>% as.factor ] %>%
    setkey(Loc)
  
  
  
  
  # Prepare the x axis
  pdAxis <- pd %>% 
    .[, list(center = ( max(Loc) + min(Loc) ) / 2 ), by=list(Chr, ChrNo)] %>%
    # omit some chromosomes
    .[, Chr := ifelse(ChrNo %in% seq(15, 21, 2), "", Chr)]
  
  
  
  
  # Prepare the core plot
  p1 <- pd %>% 
    ggplot(., aes(x=Loc, y=-1 * sign(Diff) * log10(P))) +
    # Show all points
    geom_point( aes(color=Color), alpha=0.5, size=0.5) + 
    # Plot bigger target points?
    geom_point(data = pd[Q < qThresh], aes(color=Color), size = 1) + 
    # custom X axis:
    scale_x_continuous( label = pdAxis$Chr, breaks= pdAxis$center ) +
    scale_color_manual( values = c("black", "gray")) +
    guides(color = FALSE)
  
  
  
  
  # Add labels
  if (!is.null(pd$Label) & sum(pd$Label) > 0) {
    require("ggrepel")
    # Add highlighted points
    p1 <- p1 + 
      geom_point(data=pd[Label == TRUE], color="orange", size=1) +
      geom_label_repel(data=pd[Label == TRUE], aes(label=Gene), size=6, segment.colour = "darkgreen", colour = "darkgreen")
  }
  
  
  
  # Label the axes
  p1 <- p1 + 
    xlab("Chromosome") + 
    ylab("SLP")
  
  
  
  
  # Enrichment barplot
  f <- 
    tryCatch({
      t <- fit[, table(
        Hypo=sign(logFC)<0, 
        Significant = adj.P.Val < qThresh)]
      f <- t %>% fisher.test    
      f$estimate <- format(f$estimate, digits=2, scientific=FALSE)
      f$p.value <- sprintf("p = %s", texP(f$p.value))
      f
    }, error = function(e) {
      return(list(
        estimate = 1,
        p.value = NA
      ))
    })
  
  
  
  p2 <- t %>% reshape2::melt() %>% as.data.table %>% 
    .[, list(Perc = value/sum(value), Hypo), by=Significant] %>%
    .[, Significant := factor(Significant, levels=c(FALSE, TRUE), labels=c("Background", "Significant"))] %>% 
    .[, Hypo := factor(Hypo, levels = c(FALSE, TRUE), labels=c("Hyper", "Hypo"))] %>% 
    .[Significant == "Significant"] %>% 
    ggplot(., aes(Hypo, Perc, fill=Hypo)) + 
    geom_bar(stat="identity",width = 0.5) + 
    # scale_y_continuous(labels = scales::percent, limit=c(0, 0.75)) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual("", values=c("#33a02c", "#1f78b4")) + 
    xlab("-methylated") + 
    ylab("mC Proportion, %") + 
    ggtitle("Fisher's test", latex2exp::TeX(f$p.value)) + 
    #ggtitle(glue::glue("Fisher test\np={texP(f$p.value)}")) + 
    
    guides(fill=FALSE)
  
  
  
  
  
  
  
  
  
  # Set theme
  p1 <- p1 +  
    #theme_bw(base_size=8) +
    
    theme(axis.line = element_blank(),#element_line(colour = "black")
          axis.ticks = element_line(size = 1),
          axis.ticks.length = unit(.4, "cm"),
          panel.background = element_blank(), plot.title = element_text(size = 20, family="Helvetica", face = "plain"), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text = element_text(size = 20),
          axis.title = element_text(),
          axis.text.x=element_text(size = 15, colour = "black", family="Helvetica", face = "plain"),
          axis.text.y=element_text(size =20, colour = "black", family="Helvetica", face = "plain"),
          legend.position = "none") +
    geom_hline(yintercept= log10(0.05)  , linetype="dashed", 
               color = "darkred", size=1) +
    geom_hline(yintercept= -log10(0.05)  , linetype="dashed", 
               color = "darkred", size=1)
  
  
  p2 <- p2 + 
    theme_bw(base_size=8) + 
    theme( 
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text=element_text(size=25),
      axis.title=element_text(size=25),
      plot.title = element_text(size = 25),
      plot.subtitle = element_text(size = 23)
      
      
    )
  return(list(p1, p2))
}


input = diagnosis_cd8
process_table = function(input, file_name) {
  test_table = input
  
  test_table$Chr = paste("chr",input$CHR,sep = "")
  
  test_table$Gene = test_table$gene
  #test_table$ID = paste("chr",test_table$Chr, sep = "")
  test_table$ID = test_table$Chr
  test_table$ID = paste(test_table$ID, test_table$MAPINFO, sep = "_")
  test_table$SNP =test_table$MAPINFO
  test_table$Str = ifelse(test_table$Strand == "F", "+","-")
  test_table = test_table[,c("Chr","Gene","ID","SNP","Str","logFC","t","P.Value","adj.P.Val","B")]
  
  test_table = data.table(test_table)
  
  pdf(file_name, height = 5, width = 9)
  manhattanPlot(fit = test_table, qThresh = 0.05, labelGenes = TRUE)
  #dev.off()
}

####################??????




# Takes a p value and formats with Latex in scientific notation
# but only if the value is less than 0.01
texP <- function(p, digits=2, scthreshold = 0.01) {
  p[ p == 0 ] <- 2.2e-16
  scientific <- p < scthreshold
  scientific[is.na(p)] <- FALSE
  output <- character(length(p))
  
  
  
  output[!scientific] <- 
    format(p[!scientific], digits=digits)
  
  
  
  #Transforms the number into scientific notation even if small
  output[scientific] <- 
    format(p[scientific], digits = digits, scientific = TRUE) 
  output[scientific] <- strsplit(output[scientific], split = "e")
  output[scientific] <- 
    sapply(output[scientific], 
           function(X) {
             sprintf("$%s\\,\\mathrm{x}\\,10^{%s}$", 
                     X[1], as.integer(X[2]))
           })
  output[!scientific] <- 
    sapply(output[!scientific],
           function(X) {
             sprintf("$%s$", X)
           })
  output[ is.na(p) ] <- "NA"
  return(unlist(output))
}


####################################RUN FUNCTION########################

setwd("/Users/Aoji.Xie/Desktop/DNA_methylation_final/butyric")

res = data.frame(fread(file = "methylation_results_PD_butyric.txt"))

library(ggplot2)
library("ggthemes")
process_table(input = res, file_name ="methylation_results_PD_butyric_plot_test.pdf")
dev.off()


