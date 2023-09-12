#This script uses the output of bedtools intersected replicates for peak annotation
#in Galaxy: individual BAMs -> MACS peak call -> bedtools intersect
##### Annotation of Peaks ####
#libraries
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
#set working directory
setwd("/hpc/group/goldsteinlab/tbk13/PRC2-iv-ChIP/in-vivo-intersected/")

#H3K27me3 peaks
me3broadPeaks <- readPeakFile(peakfile = "./Galaxy129-[H3K27me3_consensus-bedtools_Intersect_intervals_on_data_125_and_data_119].bed")
levels(me3broadPeaks@seqnames@values) <- c("chr1", "chr10","chr11",
                                           "chr12","chr13","chr14",
                                           "chr15","chr16","chr17",
                                           "chr18","chr19","chr2",
                                           "chr3","chr4","chr5",
                                           "chr6","chr7","chr8",
                                           "chr9","chrX","chrY")
peakAnno_me3 <- annotatePeak(me3broadPeaks, 
                             tssRegion=c(-5000, 5000), 
                             TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, 
                             annoDb="org.Mm.eg.db")
me3_annoPeaks <- as.data.frame(peakAnno_me3@anno)

#H3K27ac peaks
acbroadPeaks <- readPeakFile(peakfile = "./Galaxy128-[H3K27ac_consensus-bedtools_Intersect_intervals_on_data_121_and_data_123].bed")
levels(acbroadPeaks@seqnames@values) <- c("chr1", "chr10","chr11",
                                          "chr12","chr13","chr14",
                                          "chr15","chr16","chr17",
                                          "chr18","chr19","chr2",
                                          "chr3","chr4","chr5",
                                          "chr6","chr7","chr8",
                                          "chr9","chrX","chrY")
peakAnno_ac <- annotatePeak(acbroadPeaks, 
                            tssRegion=c(-5000, 5000), 
                            TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, 
                            annoDb="org.Mm.eg.db")
ac_annoPeaks <- as.data.frame(peakAnno_ac@anno)

write.csv(me3_annoPeaks, "~/me3_annopeaks.csv")
write.csv(ac_annoPeaks, "~/ac_annopeaks.csv")

##### Non-duplicated gene lists #####
genesMe3 <- me3_annoPeaks %>%
  select(geneId, ENSEMBL, SYMBOL, GENENAME) %>% 
  filter(!duplicated(ENSEMBL))
genesAc <- ac_annoPeaks%>%
  select(geneId, ENSEMBL, SYMBOL, GENENAME) %>% 
  filter(!duplicated(ENSEMBL))

write.csv(genesMe3, file = "./genesMe3_intersected.csv")
write.csv(genesAc, file = "./genesAc_intersected.csv")

#H3K27Ac and H3K27me3 marked genes
genesAcMe3_bivalent <- semi_join(genesAc, genesMe3, by = c("ENSEMBL"))
genesMe3_only <- anti_join(genesMe3, genesAc, by = c("ENSEMBL"))
genesAc_only <- anti_join(genesAc, genesMe3, by = c("ENSEMBL"))

write.csv(genesAcMe3_bivalent, "./genesAcMe3_bivalent_intersected.csv")
write.csv(genesAc_only, "./genesAc_only_intersected.csv")
write.csv(genesMe3_only, "./genesMe3_only_intersected.csv")

##### Toppgene graphs #####
library(tidyr)
library(forcats)
#load Toppgene tables- generated from Toppgene on non-duplicated gene lists
Toppgene_genesMe3_only <- read_delim("./Toppgene_genesMe3_only_intersected.txt", 
                                     delim = "\t", escape_double = FALSE, 
                                     trim_ws = TRUE)
Toppgene_genesAc_only <- read_delim("./Toppgene_genesAc_only_intersected.txt", 
                                    delim = "\t", escape_double = FALSE, 
                                    trim_ws = TRUE)
Toppgene_genes_AcMe3_bivalent <- read_delim("./Toppgene_genes_AcMe3_bivalent_intersected.txt", 
                                            delim = "\t", escape_double = FALSE, 
                                            trim_ws = TRUE)

#make the (GO:#####) column
Toppgene_genesMe3_only<-Toppgene_genesMe3_only %>% 
  select(!`Hit in Query List`) %>% 
  mutate(ID = paste0("(",ID,")"))
Toppgene_genesAc_only<-Toppgene_genesAc_only %>% 
  select(!`Hit in Query List`) %>% 
  mutate(ID = paste0("(",ID,")"))
Toppgene_genes_AcMe3_bivalent<-Toppgene_genes_AcMe3_bivalent %>% 
  select(!`Hit in Query List`) %>% 
  mutate(ID = paste0("(",ID,")"))

Toppgene_genesMe3_only<-Toppgene_genesMe3_only %>% 
  mutate(Name_ID = paste(Name, ID, sep = " "))
Toppgene_genesAc_only<-Toppgene_genesAc_only %>% 
  mutate(Name_ID = paste(Name, ID, sep = " "))
Toppgene_genes_AcMe3_bivalent<-Toppgene_genes_AcMe3_bivalent%>% 
  mutate(Name_ID = paste(Name, ID, sep = " "))

#make the gene proportion column (x count in query list/ y count in genome)
Toppgene_genesMe3_only <- Toppgene_genesMe3_only %>% 
  mutate(geneProportion = paste0("(",`Hit Count in Query List`, "/", `Hit Count in Genome`, ")"))
Toppgene_genesAc_only <- Toppgene_genesAc_only %>% 
  mutate(geneProportion = paste0("(",`Hit Count in Query List`, "/", `Hit Count in Genome`, ")"))
Toppgene_genes_AcMe3_bivalent <- Toppgene_genes_AcMe3_bivalent %>% 
  mutate(geneProportion = paste0("(",`Hit Count in Query List`, "/", `Hit Count in Genome`, ")"))

#data visualization
Top15_tg_me3_only <- c("synaptic signaling",
                       "neuropeptide signaling pathway",
                       "anterograde trans-synaptic signaling",
                       'chemical synaptic transmission',
                       'regulation of membrane potential',
                       'trans-synaptic signaling',
                       'regulation of postsynaptic membrane potential',
                       'nervous system process',
                       'cell fate commitment',
                       'central nervous system development',
                       'neuron fate commitment',
                       'chemical synaptic transmission, postsynaptic',
                       'excitatory postsynaptic potential',
                       'action potential',
                       'neuron fate specification'
)
plotTgMe3only <- Toppgene_genesMe3_only %>% 
  filter(Name %in% Top15_tg_me3_only) %>% 
  mutate(Name_ID = fct_reorder(Name_ID, -log10(`q-value Bonferroni`))) %>% 
  ggplot(aes(x = Name_ID,y= -log10(`q-value Bonferroni`)))+geom_col(fill="#D41159")+
  geom_text(aes(label = geneProportion),nudge_y = 0.1, size = 5, hjust = 1.1)+
  ylab("-log(q value)")+
  theme_classic()+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 22),
        axis.text.x = element_text(size = 18),
        plot.margin = margin(r = 15, b = 10, t = 10, l = 5))+
  coord_flip()
plotTgMe3only #1130x541

Top15_tg_ac_only <- c('macromolecule catabolic process',
                      'regulation of organelle organization',
                      'protein transport',
                      'DNA metabolic process',
                      'chromosome organization',
                      'mitotic cell cycle',
                      'regulation of cell cycle',
                      'cell cycle phase transition',
                      'vesicle organization',
                      'microtubule cytoskeleton organization',
                      'cell division',
                      'DNA replication',
                      'chromatin remodeling',
                      'histone modification',
                      'epigenetic regulation of gene expression')
plotTgAcOnly <-Toppgene_genesAc_only %>% 
  filter(Name %in% Top15_tg_ac_only) %>% 
  mutate(Name_ID = fct_reorder(Name_ID, -log10(`q-value Bonferroni`))) %>% 
  ggplot(aes(x = Name_ID,y= -log10(`q-value Bonferroni`)))+geom_col(fill = "#1A85FF")+
  geom_text(aes(label = geneProportion), nudge_y = 0.1, size = 5, hjust = 1.1)+
  ylab("-log(P value)")+
  theme_classic()+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 22),
        axis.text.x = element_text(size = 18),
        plot.margin = margin(r = 25,b = 10, t = 10, l = 5))+
  coord_flip()
plotTgAcOnly

Top15_tg_acMe3 <- c('neuron development',
                    'cell morphogenesis',
                    'sensory organ development',
                    'cell adhesion',
                    'central nervous system development',
                    'epithelium development',
                    'axon development',
                    'embryonic morphogenesis',
                    'synapse organization',
                    'response to growth factor',
                    'actin cytoskeleton organization',
                    'regulation of cell development',
                    'regulation of nervous system development',
                    'axon guidance',
                    'morphogenesis of an epithelium')
plotTgAcMe3 <-Toppgene_genes_AcMe3_bivalent %>% 
  filter(Name %in% Top15_tg_acMe3) %>% 
  mutate(Name_ID = fct_reorder(Name_ID, -log10(`q-value Bonferroni`))) %>% 
  ggplot(aes(x = Name_ID,y= -log10(`q-value Bonferroni`)))+geom_col(fill = "#A73FF3")+
  geom_text(aes(label = geneProportion), nudge_y = 0.1, size = 5, hjust = 1.1)+
  ylab("-log(P value)")+
  theme_classic()+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 22),
        axis.text.x = element_text(size = 18),
        plot.margin = margin(r = 25,b = 10, t = 10, l = 5))+
  coord_flip()
plotTgAcMe3


##### Transcription factor lists ####
library(readxl)
tfList <- read_excel("in-vivo-intersected/41467_2017_BFncomms15089_MOESM3449_ESM.xlsx", 
                     sheet = "allTFs")
tfNamesOnly <- tfList %>% 
  mutate(SYMBOL = GeneSymbol) %>% 
  select(SYMBOL) 

genesAc_only_intersected <- read_csv("in-vivo-intersected/genesAc_only_intersected.csv")
genesAcMe3_bivalent_intersected <- read_csv("in-vivo-intersected/genesAcMe3_bivalent_intersected.csv")
genesMe3_only_intersected <- read_csv("in-vivo-intersected/genesMe3_only_intersected.csv")

acTFs <- semi_join(genesAc_only_intersected, tfNamesOnly, by = c("SYMBOL"))
acMe3TFs <- semi_join(genesAcMe3_bivalent_intersected, tfNamesOnly, by = c("SYMBOL"))
me3TFs <- semi_join(genesMe3_only_intersected, tfNamesOnly, by = c('SYMBOL'))

write.csv(me3TFs, "me3TFs.csv")
write.csv(acMe3TFs, "acMe3TFs.csv")
write.csv(acTFs, "acTFs.csv")
