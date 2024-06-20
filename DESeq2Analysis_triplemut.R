setwd("/binf-isilon/PBgrp/xpj980/datascratch/triplemutants/04.STAR_mappings/old_style_mapping/test/")

library(tidyverse)
library(DESeq2)
library(RColorBrewer)
#install.packages("pheatmap")
library(pheatmap)
library(ggpubr)
library(ggrepel)
library(rtracklayer) 
library(VennDiagram)
library(png)
library(systemPipeR)
library(eulerr)
library(GeneOverlap)
library(svglite)
library(reshape2)
library(wesanderson)
library(GenomicRanges)
library(GenomicAlignments)
library(readxl)
library(magrittr)
library(forcats)
library(ComplexUpset)
library(openxlsx)

#### miRNA target list ### 

miRNAtargets_rawtable <- readxl::read_xlsx("/binf-isilon/PBgrp/xpj980/TAIR/miRNAtargetlist_Maria_degradome_done.xlsx", col_names = T)

miRNAtargets <- miRNAtargets_rawtable %>% 
  select(miRNA, 
         gene = target_id, 
         target_name, 
         miRNA_ID)

miRNAtargets <- filter(miRNAtargets, !duplicated(miRNAtargets$gene))

##import my featureCounts table of all libraries  
raw_featureCount_table <- read.table("/binf-isilon/PBgrp/xpj980/datascratch/triplemutants/04.STAR_mappings/old_style_mapping/featureCountsCarlottagtfmodifiedcountall", header = TRUE, row.names = 1)
glimpse(raw_featureCount_table)

## make the table DESeq2-usable by removing columns (chr, start, end, strand and length)
DESeq_table <- raw_featureCount_table %>%
  rownames_to_column() %>%
  select(1,7:50) %>% 
  column_to_rownames(var = "rowname") %>%
  glimpse()


new_names <- names(DESeq_table) %>% 
  str_remove(pattern = ".fastq.gz.trimmed.fq.gzAligned.out.sam")

names(DESeq_table) <- new_names
glimpse(DESeq_table)

## reorder columns 
colorder <- c("wt_R1", "wt_R2", "wt_R3", "ski2.2_R1", "ski2.2_R2", "ski2.2_R3", "ski2.5_R1", "ski2.5_R2", "ski2.5_R3", "rrp4_R1", "rrp4_R2", "rrp4_R3", "urt1.1_R1", "urt1.1_R2", "urt1.1_R3", "urt1.2_R1", "urt1.2_R2", "urt1.2_R3", "heso1.3_R1", "heso1.3_R2", "heso1.3_R3", "GKheso1_R1", "GKheso1_R2", "GKheso1_R3", "heso1urt1_R1", "heso1urt1_R2", "heso1urt1_R3", "rrp4urt1heso1_R1", "rrp4urt1heso1_R2", "rrp4urt1heso1_R3", "rrp4urt1dcl2.A1_R1", "rrp4urt1dcl2.A1_R2", "rrp4urt1dcl2.A1_R3", "rrp4urt1dcl2.B9_R1", "rrp4urt1dcl2.B9_R2", "rrp4urt1dcl2.B9_R3", "rrp4urt1dcl2het_R1", "rrp4urt1dcl2het_R2", "rrp4urt1dcl2het_R3", "ski2urt1dcl2.107_R1", "ski2urt1dcl2.107_R2", "ski2urt1dcl2.123_R1", "ski2urt1dcl2.123_R2", "ski2urt1dcl2.123_R3")
reordered_deseq_table <- DESeq_table[,colorder]

## make metadata matrix
introduce_genotype_to_table  <- data.frame(row.names=colnames(reordered_deseq_table),
                                           "genotype"=c(rep("wt", 3),
                                                        rep("ski2.2", 3), 
                                                        rep("ski2.5", 3),
                                                        rep("rrp4", 3), 
                                                        rep("urt1.1", 3), 
                                                        rep("urt1.2", 3), 
                                                        rep("heso1.3", 3), 
                                                        rep("GKheso1", 3), 
                                                        rep("heso1urt1", 3), 
                                                        rep("rrp4urt1heso1", 3), 
                                                        rep("rrp4urt1dcl2.A1", 3), 
                                                        rep("rrp4urt1dcl2.B9", 3), 
                                                        rep("rrp4urt1dcl2het", 3), 
                                                        rep("ski2urt1dcl2.107", 2), 
                                                        rep("ski2urt1dcl2.123", 3)))


str(introduce_genotype_to_table)
levels(introduce_genotype_to_table)
introduce_genotype_to_table$genotype  <- as.factor(introduce_genotype_to_table$genotype)
introduce_genotype_to_table$genotype <- relevel(introduce_genotype_to_table$genotype, ref="wt")
levels(introduce_genotype_to_table$genotype)

design  <- model.matrix(~genotype, data=introduce_genotype_to_table)

## convert featurecount table into a matrix
DESeq_matrix <- as.matrix(reordered_deseq_table)
head(DESeq_matrix)

## check that sample order is the same in my metadata matrix and in my featurecount matrix
all(rownames(design) == colnames(DESeq_matrix))

## make the DESeq2 object 
dds <- DESeqDataSetFromMatrix(countData=DESeq_matrix, colData=introduce_genotype_to_table, design = ~ genotype)
head(dds)

##########QA before DESeq analysis ########
## normalize counts & extract them
dds2 <- estimateSizeFactors(dds)
sizeFactors(dds2)
normalized_counts <- counts(dds2, normalized = TRUE)

## log_transformation with vst function 
log_transformeddds2 <- vst(dds2, blind = TRUE)

## distance matrix ###

sampleDists <- dist(t(assay(log_transformeddds2)))

sampleDistMatrix <- as.matrix(sampleDists)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
test <- pheatmap(sampleDistMatrix,
                 clustering_distance_rows=sampleDists,
                 clustering_distance_cols=sampleDists,
                 col=colors)

#ggsave(test, file = "/binf-isilon/PBgrp/xpj980/figures/triplemutants/distanceMatrix20210518.svg", units = "in", width = 15, height = 15, dpi = 300)

## for the TUTase paper: Figure 2A) A PCA plot with WT, ski2, rrp4, urt1, heso1, dbmut ###

DESeq_table <- raw_featureCount_table %>%
  rownames_to_column() %>%
  select(1,7:9,13:18,34:36,42:44,48:50) %>% 
  column_to_rownames(var = "rowname") %>%
  glimpse()


new_names <- names(DESeq_table) %>% 
  str_remove(pattern = ".fastq.gz.trimmed.fq.gzAligned.out.sam")

names(DESeq_table) <- new_names
glimpse(DESeq_table)

## reorder columns 
colorder <- c("wt_R1", "wt_R2", "wt_R3", "ski2.5_R1", "ski2.5_R2", "ski2.5_R3", "rrp4_R1", "rrp4_R2", "rrp4_R3", "urt1.1_R1", "urt1.1_R2", "urt1.1_R3", "GKheso1_R1", "GKheso1_R2", "GKheso1_R3", "heso1urt1_R1", "heso1urt1_R2", "heso1urt1_R3")
reordered_deseq_table <- DESeq_table[,colorder]

## make metadata matrix
introduce_genotype_to_table  <- data.frame(row.names=colnames(reordered_deseq_table),
                                           "genotype"=c(rep("wt", 3),
                                                        rep("ski2.5", 3),
                                                        rep("rrp4", 3), 
                                                        rep("urt1.1", 3), 
                                                        rep("GKheso1", 3), 
                                                        rep("heso1urt1", 3)))


str(introduce_genotype_to_table)
levels(introduce_genotype_to_table)
introduce_genotype_to_table$genotype  <- as.factor(introduce_genotype_to_table$genotype)
introduce_genotype_to_table$genotype <- relevel(introduce_genotype_to_table$genotype, ref="wt")
levels(introduce_genotype_to_table$genotype)

design  <- model.matrix(~genotype, data=introduce_genotype_to_table)

## convert featurecount table into a matrix
DESeq_matrix <- as.matrix(reordered_deseq_table)
head(DESeq_matrix)

## check that sample order is the same in my metadata matrix and in my featurecount matrix
all(rownames(design) == colnames(DESeq_matrix))

## make the DESeq2 object 
dds <- DESeqDataSetFromMatrix(countData=DESeq_matrix, colData=introduce_genotype_to_table, design = ~ genotype)
head(dds)

## normalize counts & extract them
dds2 <- estimateSizeFactors(dds)
sizeFactors(dds2)
normalized_counts <- counts(dds2, normalized = TRUE)

## log_transformation with vst function 
log_transformeddds2 <- vst(dds2, blind = TRUE)

####### PCA plot Figure 2A #######
PCAplot_genotype <- plotPCA(log_transformeddds2, intgroup = "genotype", returnData = T)

percentVar <- round(100 * attr(PCAplot_genotype, "percentVar"))

to_plot <- separate(PCAplot_genotype, name, into = c("genotype", "replicate"), sep = "_")
#save(to_plot,file = "/binf-isilon/PBgrp/xpj980/R_scripts_collected/triplemutants_PCAplot.RData")
#write.xlsx(x = to_plot, file = "/binf-isilon/PBgrp/xpj980/R_scripts_collected/triplemutants_PCA.xlsx")

colourCount = length(unique(to_plot$genotype))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

PCAplot <- ggplot(to_plot, aes(x=PC1, y=PC2, color=factor(genotype), label=genotype, shape=replicate)) +
  geom_point(size=3)+
  labs(color="genotype")+
  #geom_text_repel(size=4, fontface = 'bold', show.legend = FALSE, max.overlaps = Inf) +
  #geom_text(aes(label=sample), size=2, hjust=0.9, vjust=1.5)+
  color_palette(getPalette((colourCount))) +
  ggtitle('PC1 and PC2') +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggpubr::theme_pubr() +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(aspect.ratio = 1) 
coord_fixed() 
ggsave(PCAplot, file = "/binf-isilon/PBgrp/xpj980/figures/triplemutants/PCAplot20231122.svg", units = "in", width = 15, height = 15, dpi = 300)

# distance matrix of samples in Figure 2
sampleDists <- dist(t(assay(log_transformeddds2)))

sampleDistMatrix <- as.matrix(sampleDists)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
test <- pheatmap(sampleDistMatrix,
                 clustering_distance_rows=sampleDists,
                 clustering_distance_cols=sampleDists,
                 col=colors)

#ggsave(test, file = "/binf-isilon/PBgrp/xpj980/figures/triplemutants/distanceMatrix20231122.svg", units = "in", width = 15, height = 15, dpi = 300)

# tidy up samples for DESeq2 analysis 

DESeq_table <- raw_featureCount_table %>%
  rownames_to_column() %>%
  select(1,7:18,20,22:50) %>% 
  column_to_rownames(var = "rowname") %>%
  glimpse()

new_names <- names(DESeq_table) %>% 
  str_remove(pattern = ".fastq.gz.trimmed.fq.gzAligned.out.sam")

names(DESeq_table) <- new_names
glimpse(DESeq_table)

## reorder columns 
colorder <- c("wt_R1", "wt_R2", "wt_R3", "ski2.2_R1", "ski2.2_R2", "ski2.2_R3", "ski2.5_R1", "ski2.5_R2", "ski2.5_R3", "rrp4_R1", "rrp4_R2", "rrp4_R3", "urt1.1_R1", "urt1.1_R2", "urt1.1_R3", "urt1.2_R1", "urt1.2_R2", "urt1.2_R3", "heso1.3_R1", "heso1.3_R2", "heso1.3_R3", "GKheso1_R1", "GKheso1_R2", "GKheso1_R3", "heso1urt1_R1", "heso1urt1_R2", "heso1urt1_R3", "rrp4urt1heso1_R1", "rrp4urt1heso1_R2", "rrp4urt1heso1_R3", "rrp4urt1dcl2.A1_R2", "rrp4urt1dcl2.B9_R1", "rrp4urt1dcl2.B9_R2", "rrp4urt1dcl2.B9_R3", "rrp4urt1dcl2het_R1", "rrp4urt1dcl2het_R2", "rrp4urt1dcl2het_R3", "ski2urt1dcl2.107_R1", "ski2urt1dcl2.107_R2", "ski2urt1dcl2.123_R1", "ski2urt1dcl2.123_R2", "ski2urt1dcl2.123_R3")
reordered_deseq_table <- DESeq_table[,colorder]

## make metadata matrix
introduce_genotype_to_table  <- data.frame(row.names=colnames(reordered_deseq_table),
                                           "genotype"=c(rep("wt", 3),
                                                        rep("ski2.2", 3), 
                                                        rep("ski2.5", 3),
                                                        rep("rrp4", 3), 
                                                        rep("urt1.1", 3), 
                                                        rep("urt1.2", 3), 
                                                        rep("heso1.3", 3), 
                                                        rep("GKheso1", 3), 
                                                        rep("heso1urt1", 3), 
                                                        rep("rrp4urt1heso1", 3), 
                                                        rep("rrp4urt1dcl2.A1", 1), 
                                                        rep("rrp4urt1dcl2.B9", 3), 
                                                        rep("rrp4urt1dcl2het", 3), 
                                                        rep("ski2urt1dcl2.107", 2), 
                                                        rep("ski2urt1dcl2.123", 3)))


str(introduce_genotype_to_table)
levels(introduce_genotype_to_table)
introduce_genotype_to_table$genotype  <- as.factor(introduce_genotype_to_table$genotype)
introduce_genotype_to_table$genotype <- relevel(introduce_genotype_to_table$genotype, ref="wt")
levels(introduce_genotype_to_table$genotype)

design  <- model.matrix(~genotype, data=introduce_genotype_to_table)

## convert featurecount table into a matrix
DESeq_matrix <- as.matrix(reordered_deseq_table)
head(DESeq_matrix)

## check that sample order is the same in my metadata matrix and in my featurecount matrix
all(rownames(design) == colnames(DESeq_matrix))

## make the DESeq2 object 
dds <- DESeqDataSetFromMatrix(countData=DESeq_matrix, colData=introduce_genotype_to_table, design = ~ genotype)
head(dds)

## normalize counts & extract them
dds2 <- estimateSizeFactors(dds)
sizeFactors(dds2)
normalized_counts <- counts(dds2, normalized = TRUE)

## log_transformation with vst function 
log_transformeddds2 <- vst(dds2, blind = TRUE)

##### PCA plot ###
PCAplot_genotype <- plotPCA(log_transformeddds2, intgroup = "genotype", returnData = T)

percentVar <- round(100 * attr(PCAplot_genotype, "percentVar"))

to_plot <- separate(PCAplot_genotype, name, into = c("genotype", "replicate"), sep = "_")

colourCount = length(unique(to_plot$genotype))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

PCAplot <- ggplot(to_plot, aes(x=PC1, y=PC2, color=factor(genotype), label=genotype, shape=replicate)) +
  geom_point(size=3)+
  labs(color="genotype")+
  #geom_text_repel(size=4, fontface = 'bold', show.legend = FALSE, max.overlaps = Inf) +
  #geom_text(aes(label=sample), size=2, hjust=0.9, vjust=1.5)+
  color_palette(getPalette((colourCount))) +
  ggtitle('PC1 and PC2') +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggpubr::theme_pubr() +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(aspect.ratio = 1) 
coord_fixed() 

## distance matrix ##

sampleDists <- dist(t(assay(log_transformeddds2)))

sampleDistMatrix <- as.matrix(sampleDists)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
test <- pheatmap(sampleDistMatrix,
                 clustering_distance_rows=sampleDists,
                 clustering_distance_cols=sampleDists,
                 col=colors)

#ggsave(test, file = "/binf-isilon/PBgrp/xpj980/figures/triplemutants/distanceMatrix20210803.svg", units = "in", width = 15, height = 15, dpi = 300)

########## DESeq Analysis #########
#deseqqed <- DESeq(dds2)

#resultsNames(deseqqed)

## plot estimated dispersions 
#svg(file = "/binf-isilon/PBgrp/xpj980/figures/triplemutants/dispersion.svg")
#(dispersionEst <- plotDispEsts(deseqqed))
#dev.off()

## extract results 
#DESeqresultsSki2.2 <- DESeq2::results(deseqqed, contrast = c("genotype", "ski2.2", "wt"), alpha = 0.05)
#DESeqresultsSki2.5 <- DESeq2::results(deseqqed, contrast = c("genotype", "ski2.5", "wt"), alpha = 0.05)
#DESeqresultsRrp4 <- DESeq2::results(deseqqed, contrast = c("genotype", "rrp4", "wt"), alpha = 0.05)
#DESeqresultsUrt1.1 <- DESeq2::results(deseqqed, contrast = c("genotype", "urt1.1", "wt"), alpha = 0.05)
#DESeqresultsUrt1.2 <- DESeq2::results(deseqqed, contrast = c("genotype", "urt1.2", "wt"), alpha = 0.05)
#DESeqresultsHeso1.3 <- DESeq2::results(deseqqed, contrast = c("genotype", "heso1.3", "wt"), alpha = 0.05)
#DESeqresultsGKheso <- DESeq2::results(deseqqed, contrast = c("genotype", "GKheso1", "wt"), alpha = 0.05)
#DESeqresultsHeso1Urt1 <- DESeq2::results(deseqqed, contrast = c("genotype", "heso1urt1", "wt"), alpha = 0.05)
#DESeqresultsRrp4Urt1Heso1 <- DESeq2::results(deseqqed, contrast = c("genotype", "rrp4urt1heso1", "wt"), alpha = 0.05)
#DESeqresultsRrp4Urt1Dcl2.A1 <- DESeq2::results(deseqqed, contrast = c("genotype", "rrp4urt1dcl2.A1", "wt"), alpha = 0.05)
#DESeqresultsRrp4Urt1Dcl2.B9 <- DESeq2::results(deseqqed, contrast = c("genotype", "rrp4urt1dcl2.B9", "wt"), alpha = 0.05)
#DESeqresultsRrp4Urt1Dcl2het <- DESeq2::results(deseqqed, contrast = c("genotype", "rrp4urt1dcl2het", "wt"), alpha = 0.05)
#DESeqresultsSki2Urt1Dcl2.107 <- DESeq2::results(deseqqed, contrast = c("genotype", "ski2urt1dcl2.107", "wt"), alpha = 0.05)
#DESeqresultsSki2Urt1Dcl2.123 <- DESeq2::results(deseqqed, contrast = c("genotype", "ski2urt1dcl2.123", "wt"), alpha = 0.05)

###### gathered DESegresults-table #######
#ski2.2_tibble <- DESeqresultsSki2.2 %>% 
# as.data.frame() %>%
#rownames_to_column(var = "gene") %>%
#as_tibble() %>%
#mutate(genotype = rep("ski2.2")) %>%
#glimpse()

#ski2.5_tibble <- DESeqresultsSki2.5 %>% 
# as.data.frame() %>%
#rownames_to_column(var = "gene") %>%
#as_tibble() %>%
#mutate(genotype = rep("ski2.5")) %>%
#glimpse()

#rrp4_tibble <- DESeqresultsRrp4 %>% 
# as.data.frame() %>%
#rownames_to_column(var = "gene") %>%
#as_tibble() %>%
#mutate(genotype = rep("rrp4")) %>%
#glimpse()

#urt1.1_tibble <- DESeqresultsUrt1.1 %>% 
# as.data.frame() %>%
#rownames_to_column(var = "gene") %>%
#as_tibble() %>%
#mutate(genotype = rep("urt1.1")) %>%
#glimpse()

#urt1.2_tibble <- DESeqresultsUrt1.2 %>% 
# as.data.frame() %>%
#rownames_to_column(var = "gene") %>%
#as_tibble() %>%
#mutate(genotype = rep("urt1.2")) %>%
#glimpse()

#heso1.3_tibble <- DESeqresultsHeso1.3 %>% 
# as.data.frame() %>%
#rownames_to_column(var = "gene") %>%
#as_tibble() %>%
#mutate(genotype = rep("heso1.3")) %>%
#glimpse()

#GKheso_tibble <- DESeqresultsGKheso %>% 
#as.data.frame() %>%
#rownames_to_column(var = "gene") %>%
#as_tibble() %>%
#mutate(genotype = rep("GKheso1")) %>%
#glimpse()

#heso1urt1_tibble <- DESeqresultsHeso1Urt1 %>% 
# as.data.frame() %>%
#rownames_to_column(var = "gene") %>%
#as_tibble() %>%
#mutate(genotype = rep("heso1urt1")) %>%
#glimpse()

#rrp4urt1heso1_tibble <- DESeqresultsRrp4Urt1Heso1 %>% 
# as.data.frame() %>%
#rownames_to_column(var = "gene") %>%
#as_tibble() %>%
#mutate(genotype = rep("rrp4urt1heso1")) %>%
#glimpse()

#rrp4urt1dcl2.A1_tibble <- DESeqresultsRrp4Urt1Dcl2.A1 %>% 
# as.data.frame() %>%
#rownames_to_column(var = "gene") %>%
#as_tibble() %>%
#mutate(genotype = rep("rrp4urt1dcl2.A1")) %>%
#glimpse()

#rrp4urt1dcl2.B9_tibble <- DESeqresultsRrp4Urt1Dcl2.B9 %>% 
# as.data.frame() %>%
#rownames_to_column(var = "gene") %>%
#as_tibble() %>%
#mutate(genotype = rep("rrp4urt1dcl2.B9")) %>%
#glimpse()

#rrp4urt1dcl2het_tibble <- DESeqresultsRrp4Urt1Dcl2het %>% 
#  as.data.frame() %>%
# rownames_to_column(var = "gene") %>%
#as_tibble() %>%
#mutate(genotype = rep("rrp4urt1dcl2het")) %>%
#glimpse()

#ski2urt1dcl2.107_tibble <- DESeqresultsSki2Urt1Dcl2.107 %>% 
# as.data.frame() %>%
#  rownames_to_column(var = "gene") %>%
# as_tibble() %>%
#mutate(genotype = rep("ski2urt1dcl2.107")) %>%
#glimpse()

#ski2urt1dcl2.123_tibble <- DESeqresultsSki2Urt1Dcl2.123 %>% 
# as.data.frame() %>%
#rownames_to_column(var = "gene") %>%
#as_tibble() %>%
#mutate(genotype = rep("ski2urt1dcl2.123")) %>%
#glimpse()

#gathered_DESeqresults <- bind_rows(... = ski2.2_tibble, ski2.5_tibble, rrp4_tibble, urt1.1_tibble, urt1.2_tibble, heso1.3_tibble, GKheso_tibble, heso1urt1_tibble, rrp4urt1heso1_tibble,rrp4urt1dcl2.A1_tibble, rrp4urt1dcl2.B9_tibble, rrp4urt1dcl2het_tibble, ski2urt1dcl2.107_tibble, ski2urt1dcl2.123_tibble) 

#save(gathered_DESeqresults, 
#   file = "/binf-isilon/PBgrp/xpj980/R_scripts_collected/triplemutants_gatheredDESeqresults.RData")


#### load DESeqresults #####
#load(file = "/binf-isilon/PBgrp/xpj980/R_scripts_collected/triplemutants_gatheredDESeqresults.RData")
load(file = "/binf-isilon/PBgrp/xpj980/R_scripts_collected/triplemutants_gatheredDESeqresults.RData")

#### excel sheets for paper ####
araport_TE_IDs <- import.gff("/binf-isilon/PBgrp/xpj980/TAIR/Araport11_GFF3_genes_transposons.201606.gff")

TE_IDs <- araport_TE_IDs %>% 
  as_tibble() %>% 
  filter(type %in% c("transposable_element", "pseudogene", "transposable_element_gene")) %>% 
  select(type, ID, locus_type) %>% 
  mutate(locus_type = case_when(is.na(locus_type) ~ "transposable_element", 
                                T ~ locus_type)) %>% 
  rename(ID = "gene")

araport_original <- import.gff("/binf-isilon/PBgrp/xpj980/TAIR/Araport11_GFF3_genes_transposons.201606.gff") %>% 
  as.tibble() %>% 
  select(type, ID, locus_type) %>% 
  filter(type =="gene") %>% 
  rename(ID = "gene")

to_join <- bind_rows(TE_IDs, araport_original)

test <- left_join(gathered_DESeqresults, to_join, by = "gene") %>% 
  mutate(locus_type = case_when(is.na(locus_type) ~ "intergenic", 
                                T ~ locus_type)) %>% 
  filter(padj < 0.05) %>% 
  filter(locus_type !="intergenic") %>% 
  left_join(., miRNAtargets, by = "gene", copy = FALSE) %>% 
  mutate(locus_type = case_when(is.na(miRNA) ~ locus_type, 
                                TRUE ~ "miRNA_target")) %>% 
  filter(genotype != "rrp4urt1dcl2.A1") %>%
  select(gene, baseMean, log2FoldChange, pvalue, padj, locus_type, genotype) 

## for GEO submission 
test <- left_join(gathered_DESeqresults, to_join, by = "gene") %>% 
  mutate(locus_type = case_when(is.na(locus_type) ~ "intergenic", 
                                T ~ locus_type)) %>% 
  filter(padj < 0.05) %>% 
  filter(locus_type !="intergenic") %>% 
  left_join(., miRNAtargets, by = "gene", copy = FALSE) %>% 
  mutate(locus_type = case_when(is.na(miRNA) ~ locus_type, 
                                TRUE ~ "miRNA_target")) %>% 
  filter(genotype %in% c("ski2.5", "rrp4", "urt1.1", "GKheso1", "heso1urt1", "rrp4urt1dcl2.B9", "rrp4urt1heso1")) %>% 
  select(gene, baseMean, log2FoldChange, pvalue, padj, locus_type, genotype) 

#write.xlsx(x = test, file = "/binf-isilon/PBgrp/xpj980/datascratch/GEO_may2024/DESeq2ExperimentBWOintergenic.xlsx")

#### HEATMAPS ####
#### Figure 2B & 6B ####

geneset <- gathered_DESeqresults %>% 
  semi_join(miRNAtargets, by = "gene") %>% 
  filter(genotype %in% c("wt", "ski2.5", "rrp4", "GKheso1", "heso1.3", "urt1.1", "urt1.2", "heso1urt1")) %>%
  filter(padj < 0.05) %>% 
  pull(gene) %>% 
  unique()

### sizes ###
# next time just load this one 
#load("/binf-isilon/PBgrp/xpj980/R_scripts_collected/triplemutantSizesNormalized.RData")
load("/binf-isilon/PBgrp/xpj980/R_scripts_collected/triplemutantSizesNormalized2.RData")

#raw_featureCount_sizes <- read.table("/binf-isilon/PBgrp/xpj980/datascratch/triplemutants/04.STAR_mappings/old_style_mapping/by_size/bam/featureCountsNewGTFSizes", header = TRUE, row.names = 1) %>%
# glimpse()

#raw_featureCount_sizes2 <- read.table("/binf-isilon/PBgrp/xpj980/datascratch/triplemutants/04.STAR_mappings/old_style_mapping/by_size/bam/featureCountstestsize", header = TRUE, row.names = 1) %>%
# glimpse()

#size_table <- raw_featureCount_sizes2 %>% 
# select(contains(".15"), contains(".16"), contains(".17"), contains(".18"), contains(".19"), contains(".20"), contains(".21"), contains(".22"), contains(".23"), contains(".24"), contains(".25"), contains(".26"), contains(".27"), contains(".28"), contains(".29"), contains(".30")) %>%
#glimpse()

#n <- names(size_table) %>% 
# str_remove("fastq.gz.trimmed.fq.gzAligned.out.") %>% 
#str_remove(".sorted.bam") %>% 
#str_replace(pattern = "R1.", replacement = "R1_") %>% 
#str_replace(pattern = "R2.", replacement = "R2_") %>%
#str_replace(pattern = "R3.", replacement = "R3_")

#names(size_table) <- n
#glimpse(size_table)

##### make table plotable ####
#size_table_plotable <- size_table %>% 
# rownames_to_column(var = "gene") %>%
#  as.data.frame() %>% glimpse()

#keycol <- "genotype"
#valcol <- "count"
#gathercols <- c(names(size_table_plotable[,2:704]))

#test <- gather_(size_table_plotable, keycol, valcol, gathercols)

#plotable_split <- separate(test, col = genotype, into = c("genotyp", "rep", "size"), sep = "_")

#normalized_for_errorbars <- plotable_split %>% 
# group_by(genotyp,rep) %>%
#mutate(RPM = count*1000000/sum(count))

#save(normalized_for_errorbars, 
#    file = "/binf-isilon/PBgrp/xpj980/R_scripts_collected/triplemutantSizesNormalized2.RData")
#load("/binf-isilon/PBgrp/xpj980/R_scripts_collected/triplemutantSizesNormalized2.RData")

subsetted_urt1_heso1_ski2_rrp4 <- filter(normalized_for_errorbars, size %in% c(21)) %>% 
  filter(gene %in% geneset) %>% 
  filter(genotyp %in% c("wt", "ski2.5", "rrp4", "GKheso1", "urt1.1", "heso1urt1")) 

subsetted_urt1_heso1_ski2_rrp4_2 <- subsetted_urt1_heso1_ski2_rrp4 %>%
  group_by(gene, genotyp) %>% 
  summarise(meanRPM = mean(RPM)) %>% 
  spread(., genotyp, meanRPM) %>% 
  column_to_rownames(var = "gene") %>%
  filter_all(., any_vars(. > 2))

annot_cols <- subsetted_urt1_heso1_ski2_rrp4_2 %>% 
  rownames_to_column(var = "gene") %>% 
  left_join(., miRNAtargets) %>% 
  select(gene, miRNA) %>% 
  column_to_rownames(var = "gene")

column_order <- c("wt", "ski2.5", "rrp4", "urt1.1", "GKheso1", "heso1urt1")

data_reordered <- subsetted_urt1_heso1_ski2_rrp4_2[, match(column_order, colnames(subsetted_urt1_heso1_ski2_rrp4_2))]
#write.xlsx(x = data_reordered, file = "/binf-isilon/PBgrp/xpj980/R_scripts_collected/triplemutant_heatmapFig2.xlsx")

normalized_heatmap_rowclust <- pheatmap(data_reordered,
                                        annotation_row = annot_cols,
                                        #annotation_colors = ann_colors, 
                                        breaks = c(-2:2),
                                        kmeans_k = NA, 
                                        cluster_rows = 1, 
                                        cluster_cols = 0,
                                        cellwidth = 20, 
                                        cellheight = 7,
                                        fontsize_row = 8,
                                        scale = "row",
                                        cutree_rows = 3,
                                        #cutree_cols = 2,
                                        color = rev(RColorBrewer::brewer.pal(name = "RdYlBu", n=4)),
                                        main = "significant DE miRNA targets RPM")
ggsave(normalized_heatmap_rowclust, file = "/binf-isilon/PBgrp/xpj980/figures/triplemutants/heatmapRPMmiRNAtargetsheso1urt1ANDdbmut21ntRPMMoreContrast.svg", units = "in", width = 15, height = 15, dpi = 300)

## Heatmap fig6
geneset <- gathered_DESeqresults %>% 
  semi_join(miRNAtargets, by = "gene") %>% 
  filter(genotype %in% c("wt", "rrp4", "GKheso1", "urt1.1", "heso1urt1", "rrp4urt1dcl2.B9", "rrp4urt1heso1")) %>%
  filter(padj < 0.05) %>% 
  pull(gene) %>% 
  unique()

subsetted_urt1_heso1_ski2_rrp4 <- filter(normalized_for_errorbars, size %in% c(21)) %>% 
  filter(gene %in% geneset) %>% 
  filter(genotyp %in% c("wt", "rrp4", "GKheso1", "urt1.1", "heso1urt1", "rrp4urt1dcl2.B9", "rrp4urt1heso1")) 

subsetted_urt1_heso1_ski2_rrp4_2 <- subsetted_urt1_heso1_ski2_rrp4 %>%
  group_by(gene, genotyp) %>% 
  summarise(meanRPM = mean(RPM)) %>% 
  spread(., genotyp, meanRPM) %>% 
  column_to_rownames(var = "gene") %>%
  filter_all(., any_vars(. > 2))

annot_cols <- subsetted_urt1_heso1_ski2_rrp4_2 %>% 
  rownames_to_column(var = "gene") %>% 
  left_join(., miRNAtargets) %>% 
  select(gene, miRNA) %>% 
  column_to_rownames(var = "gene")

column_order <- c("wt", "rrp4", "urt1.1",  "GKheso1", "heso1urt1", "rrp4urt1dcl2.B9", "rrp4urt1heso1")

data_reordered <- subsetted_urt1_heso1_ski2_rrp4_2[, match(column_order, colnames(subsetted_urt1_heso1_ski2_rrp4_2))]
write.xlsx(x = data_reordered, file = "/binf-isilon/PBgrp/xpj980/R_scripts_collected/triplemutant_heatmapFig6.xlsx")

normalized_heatmap_rowclust <- pheatmap(data_reordered,
                                        annotation_row = annot_cols,
                                        #annotation_colors = ann_colors, 
                                        breaks = c(-2:2),
                                        kmeans_k = NA, 
                                        cluster_rows = 1, 
                                        cluster_cols = 0,
                                        cellwidth = 20, 
                                        cellheight = 7,
                                        fontsize_row = 8,
                                        scale = "row",
                                        #cutree_rows = 2,
                                        #cutree_cols = 2,
                                        color = rev(RColorBrewer::brewer.pal(name = "RdYlBu", n=4)),
                                        main = "significant DE miRNA targets RPM")
ggsave(normalized_heatmap_rowclust, file = "/binf-isilon/PBgrp/xpj980/figures/triplemutants/heatmapRPMmiRNAtargetTriplemut21ntRPM.svg", units = "in", width = 15, height = 15, dpi = 300)

######### miRNA targets barplot ################
#### Figure 6A & EV1B & EV5 #####

barplot_table <- left_join(x = gathered_DESeqresults, y = miRNAtargets, by = "gene", copy = FALSE) %>%
  filter(padj < 0.05) 

to_plot <- barplot_table %>% 
  mutate(up_or_down = case_when(log2FoldChange > 0 ~ "up", 
                                log2FoldChange < 0 ~ "down"), 
         miRNA = case_when(is.na(miRNA) ~ "other", 
                           TRUE ~ "known\n miRNA target")) %>%
  dplyr::select(gene, genotype, miRNA, up_or_down)

to_plot2 <- to_plot %>% 
  filter(genotype %in% c("ski2.5", "rrp4", "urt1.1", "GKheso1", "heso1urt1"))

to_tally <- to_plot2 %>% 
  group_by(genotype, miRNA,up_or_down) %>%
  add_tally(name = "sum") %>%
  distinct() %>%
  ungroup() %>%
  mutate(sum= case_when(up_or_down == "down" ~ -1*sum,
                        up_or_down == "up" ~ 1*sum), 
         plot = paste0(miRNA, up_or_down))

unique(to_tally$genotype)

to_tally$genotype <- factor(to_tally$genotype, levels = c("ski2.5", "rrp4", "urt1.1", "GKheso1", "heso1urt1"))
to_tally$plot <- factor(to_tally$plot, levels = c("known\n miRNA targetup", "known\n miRNA targetdown", "otherdown", "otherup"))

for2A <- to_tally %>% 
  select(genotype, miRNA, up_or_down, sum, plot) %>% 
  distinct()
#write.xlsx(x = for2A, file = "/binf-isilon/PBgrp/xpj980/R_scripts_collected/triplemutant_barplotEV1.xlsx")

(barplot_across_geno <- ggplot(for2A, aes(x= genotype, y=sum, fill = plot)) +
    geom_col(alpha = 0.8, width = 0.5) +
    geom_text(aes(label = sum), vjust=1.6, color="black", size = 3) +
    scale_fill_brewer("more or less sRNAs", 
                      palette = "RdGy",    
                      labels = c("known miRNA target - more sRNAs", "known miRNA target - less sRNAs", "other gene - less sRNAs", "other gene - more sRNAs")) +
    facet_wrap(~miRNA, nrow = 1, scales = "free") +
    labs(y= "number of genes\n", x = "gene type\n") + 
    cowplot::theme_cowplot() + 
    geom_hline(yintercept = 0) +
    scale_x_discrete(position = "bottom") + 
    ggtitle("Genes with more or less sRNAs compared to WT") + 
    theme(strip.text.x = element_text(size = 12,
                                      face = "bold",
                                      color = "#22292F"),
          strip.background = element_blank(),
          plot.title = element_text(size = 16, 
                                    face = "bold",
                                    color = "#22292F", 
                                    margin = margin(b = 8)), 
          axis.title.x = element_text(size = 10,
                                      face = "bold",
                                      color = "#22292F", 
                                      margin = margin(t = 15)),
          axis.title.y = element_text(size = 10,
                                      face = "bold",
                                      color = "#22292F"), 
          axis.text.x = element_text(size = 8, angle = 90), 
          axis.text.y = element_text(size = 8),
          legend.title = element_text(size = 9,
                                      face = "bold",
                                      color = "#22292F"), 
          legend.text = element_text(size = 8))) + 
  theme(aspect.ratio = 1)

#ggsave(barplot_across_geno, file = "/binf-isilon/PBgrp/xpj980/figures/triplemutants/barplot211223.svg", device = "svg", units = "in", width = 15, height = 15, dpi = 300)

to_plot2 <- to_plot %>% 
  filter(genotype %in% c("rrp4", "urt1.1", "GKheso1", "heso1urt1", "rrp4urt1dcl2.B9", "rrp4urt1heso1"))

to_tally <- to_plot2 %>% 
  group_by(genotype, miRNA,up_or_down) %>%
  add_tally(name = "sum") %>%
  distinct() %>%
  ungroup() %>%
  mutate(sum= case_when(up_or_down == "down" ~ -1*sum,
                        up_or_down == "up" ~ 1*sum), 
         plot = paste0(miRNA, up_or_down))

unique(to_tally$genotype)

to_tally$genotype <- factor(to_tally$genotype, levels = c("ski2.5", "rrp4", "urt1.1", "GKheso1", "heso1urt1", "rrp4urt1dcl2.B9", "rrp4urt1heso1"))
to_tally$plot <- factor(to_tally$plot, levels = c("known\n miRNA targetup", "known\n miRNA targetdown", "otherdown", "otherup"))

for2A <- to_tally %>% 
  select(genotype, miRNA, up_or_down, sum, plot) %>% 
  distinct()
#write.xlsx(x = for2A, file = "/binf-isilon/PBgrp/xpj980/R_scripts_collected/triplemutant_barplotFig6A.xlsx")

(barplot_across_geno <- ggplot(for2A, aes(x= genotype, y=sum, fill = plot)) +
    geom_col(alpha = 0.8, width = 0.5) +
    geom_text(aes(label = sum), vjust=1.6, color="black", size = 3) +
    scale_fill_brewer("more or less sRNAs", 
                      palette = "RdGy",    
                      labels = c("known miRNA target - more sRNAs", "known miRNA target - less sRNAs", "other gene - less sRNAs", "other gene - more sRNAs")) +
    facet_wrap(~miRNA, nrow = 1, scales = "free") +
    labs(y= "number of genes\n", x = "gene type\n") + 
    cowplot::theme_cowplot() + 
    geom_hline(yintercept = 0) +
    scale_x_discrete(position = "bottom") + 
    ggtitle("Genes with more or less sRNAs compared to WT") + 
    theme(strip.text.x = element_text(size = 12,
                                      face = "bold",
                                      color = "#22292F"),
          strip.background = element_blank(),
          plot.title = element_text(size = 16, 
                                    face = "bold",
                                    color = "#22292F", 
                                    margin = margin(b = 8)), 
          axis.title.x = element_text(size = 10,
                                      face = "bold",
                                      color = "#22292F", 
                                      margin = margin(t = 15)),
          axis.title.y = element_text(size = 10,
                                      face = "bold",
                                      color = "#22292F"), 
          axis.text.x = element_text(size = 8, angle = 90), 
          axis.text.y = element_text(size = 8),
          legend.title = element_text(size = 9,
                                      face = "bold",
                                      color = "#22292F"), 
          legend.text = element_text(size = 8))) + 
  theme(aspect.ratio = 1)

ggsave(barplot_across_geno, file = "/binf-isilon/PBgrp/xpj980/figures/triplemutants/barplot211223.svg", device = "svg", units = "in", width = 15, height = 15, dpi = 300)
