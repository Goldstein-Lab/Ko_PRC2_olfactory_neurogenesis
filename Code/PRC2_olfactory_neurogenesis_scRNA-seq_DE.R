library(dplyr)
library(stringr)
library(Seurat)
library(patchwork)
library(ggrepel)
library(RColorBrewer)
library(DESeq2)
library(ggplot2)
library(viridis)
library(ggpubr)
library(rstatix)

EEDKOHET_WT_GBC <-readRDS("./EEDKOHET_WT_GBC.rds")

##### KO tom+ vs tom- #####
Idents(EEDKOHET_WT_GBC)<-"genotype"
EEDKOHET_WT_GBC_Ko <- subset(EEDKOHET_WT_GBC, idents = c("ko"))

#for tom+
Idents(EEDKOHET_WT_GBC_Ko) <- "tomato"
  VlnPlot(EEDKOHET_WT_GBC_Ko, features = c('tdTomato','Gapdh','Actb','Ppia','Ubc','Hprt'))
  GBC_Ko_cluster_markers <- FindMarkers(EEDKOHET_WT_GBC_Ko, ident.1 = "positive",
                                        ident.2 = 'negative', min.pct = 0.25, logfc.threshold = 0.00)
  GBC_Ko_cluster_markers$nameRows <- row.names(GBC_Ko_cluster_markers)
  
  
  #Make Volcano Plot
  #Create categories to color
  GBC_Ko_cluster_markers$Change <- "NO"
  GBC_Ko_cluster_markers $Change[GBC_Ko_cluster_markers$avg_log2FC > 0.6 & GBC_Ko_cluster_markers$p_val < 0.05] <- "UP"
  GBC_Ko_cluster_markers $Change[GBC_Ko_cluster_markers $avg_log2FC < -0.6 & GBC_Ko_cluster_markers $p_val < 0.05] <- "DOWN"
  mycolors <- c("royalblue4", "firebrick2", "gray")
  names(mycolors) <- c("DOWN", "UP", "NO")
  
  #Now add gene names (stored in rownames) as a columc called label, only labeling genes in Up or Down categories above
  GBC_Ko_cluster_markers$label <- NA
  labels <- GBC_Ko_cluster_markers$...1
  GBC_Ko_cluster_markers$label <- labels
  GBC_Ko_cluster_markers$label[GBC_Ko_cluster_markers$Change == "NO"] <- NA
  
  write.csv(GBC_Ko_cluster_markers,"./GBC_Ko_cluster_markers.csv")
  
  GBC_Ko_cluster_markers<- GBC_Ko_cluster_markers %>% 
    filter(`...1` != "tdTomato")
  
  GBCEedKo_tomPosVsTomNeg_DEplot <- ggplot(data=GBC_Ko_cluster_markers, 
                                           aes(x=avg_log2FC, y=-log10(p_val),
                                               col=Change, label=label)) +
    geom_point(alpha = 1/1.3) + 
    theme_classic() + geom_vline(xintercept=c(-0.6, 0.6), linetype='dashed') +
    geom_hline(yintercept=-log10(0.05), linetype='dashed') + 
    scale_colour_manual(values = mycolors) + geom_text_repel(color="black") + 
    labs(title = "EED KO GBCs") +
    theme(plot.title = element_text(hjust = 0.5, size=18, face='bold'), 
          text = element_text(size = 12, face="bold"),
          legend.position = "none")
  
  
  GBCEedKo_tomPosVsTomNeg_DEplot

    ##### KO tom+ vs Het tom + #####
    Idents(EEDKOHET_WT_GBC)<-"tomato"
      EEDKOHET_GBC_Tom <- subset(EEDKOHET_WT_GBC, idents = c("positive"))
      Idents(EEDKOHET_GBC_Tom)<-"genotype"
      VlnPlot(EEDKOHET_GBC_Tom, features = c('tdTomato','Gapdh','Actb','Ppia','Ubc','Hprt'))
      GBC_HetKo_cluster_markers <- FindMarkers(EEDKOHET_GBC_Tom, ident.1 = "ko",
                                               ident.2 = 'het', min.pct = 0.25, logfc.threshold = 0.00)
      GBC_HetKo_cluster_markers$nameRows <- row.names(GBC_HetKo_cluster_markers)
      GBC_HetKo_cluster_markers$p_val_bh <- p.adjust(GBC_HetKo_cluster_markers$p_val, method = "fdr")
      #Make Volcano Plot
      #Create categories to color
      GBC_HetKo_cluster_markers$Change <- "NO"
      GBC_HetKo_cluster_markers $Change[GBC_HetKo_cluster_markers$avg_log2FC > 0.6 & GBC_HetKo_cluster_markers$p_val < 0.05] <- "UP"
      GBC_HetKo_cluster_markers $Change[GBC_HetKo_cluster_markers $avg_log2FC < -0.6 & GBC_HetKo_cluster_markers $p_val < 0.05] <- "DOWN"
      mycolors <- c("royalblue4", "firebrick2", "gray")
      names(mycolors) <- c("DOWN", "UP", "NO")
      
      #Now add gene names (stored in rownames) as a columc called label, only labeling genes in Up or Down categories above
      GBC_HetKo_cluster_markers$label <- NA
      labels <- rownames(GBC_HetKo_cluster_markers)
      GBC_HetKo_cluster_markers$label <- labels
      GBC_HetKo_cluster_markers$label[GBC_HetKo_cluster_markers$Change == "NO"] <- NA
      
      
      
      write.csv(GBC_HetKo_cluster_markers,"./GBC_HetKo_cluster_markers.csv")
      
      GBC_HetKo_cluster_markers.clean <- GBC_HetKo_cluster_markers %>% 
        filter(!grepl("Rpl", nameRows)) %>% 
        filter(!grepl("Mrpl", nameRows)) %>% 
        filter(!grepl("Mrps", nameRows)) %>% 
        filter(!grepl("^Rps", nameRows)) %>% 
        filter(!grepl("^Hba", nameRows)) %>%
        filter(!grepl("^Hbb", nameRows)) %>% 
        filter(!grepl("^mt", nameRows)) %>%
        filter(!grepl("Uty", nameRows)) %>% 
        filter(!grepl("Eif2s3y",nameRows)) %>% 
        filter(!grepl("Shroom2",nameRows))
      
      
      
      write.csv(GBC_HetKo_cluster_markers_clean,"./GBC_HetKo_cluster_markers_clean.csv")
      
      GBCEedKoVsHet_DEplot <- ggplot(data=GBC_HetKo_cluster_markers_clean, aes(x=avg_log2FC, y=-log10(p_val_bh), col=Change, label=label)) +
        geom_point(alpha = 1/1.3) + 
        theme_classic() + geom_vline(xintercept=c(-0.6, 0.6), linetype='dashed') +
        geom_hline(yintercept=-log10(0.05), linetype='dashed') + 
        scale_colour_manual(values = mycolors) + geom_text_repel(color="black") + 
        labs(title = "tdTom+ GBC- EEDKO vs Het") +
        theme(plot.title = element_text(hjust = 0.5, size=20, face='bold'), text = element_text(size = 12, face="bold"))
      
      GBCEedKoVsHet_DEplot

          