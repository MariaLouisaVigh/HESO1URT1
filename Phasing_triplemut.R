####### figure 2F ###### 

setwd("/binf-isilon/PBgrp/xpj980/datascratch/triplemutants/04.STAR_mappings/old_style_mapping/")

library(tidyverse)
library(DESeq2)
library(RColorBrewer)
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


# custum functions
rc_string <- function(string) {
  string %>% Biostrings::DNAStringSet() %>% 
    Biostrings::reverseComplement() %>% 
    as.character()
}

c("ATCG", "AGGT") %>% rc_string()

## load gtf and target list 
mygtf <- import("/binf-isilon/PBgrp/xpj980/TAIR/Arabidopsis_thaliana.TAIR10.44.gtf") 

targets <- read.table("/binf-isilon/PBgrp/xpj980/datascratch/2015_exosome_mutants/03.STAR_mapping/sam_files/annotated_ASRP_miRNA_target_sites.bed.txt", header = T, sep = "\t", quote = "", strip.white = T, stringsAsFactors = F) %>%
  as_tibble()

target_GR <- makeGRangesFromDataFrame(df = targets, keep.extra.columns = T, seqnames.field = "chr", start.field = "start", end.field = "stop", strand.field = "strand")
target_GR <- target_GR[!duplicated(target_GR)]

targets <- read_xlsx("/binf-isilon/PBgrp/xpj980/TAIR/miRNAtargetlist_Maria_degradome_done.xlsx", col_names = T)

### loop ###
target_vector <- c("AT2G39675", "AT1G62910", "AT2G34710", "AT4G18390", "AT1G48410", "AT3G26810", "AT1G63130", "AT1G50055")
genotypes_sub <- c("wt", "ski2.5", "rrp4", "urt1.1", "GKheso1", "heso1urt1")

bam_paths <- list.files("/binf-isilon/PBgrp/xpj980/datascratch/triplemutants/04.STAR_mappings/old_style_mapping/by_size/bam/", full.names = T, pattern = "Aligned.out.21.sorted.bam")
index_paths <- list.files("/binf-isilon/PBgrp/xpj980/datascratch/triplemutants/04.STAR_mappings/old_style_mapping/by_size/bam/bai/", full.names = T, pattern = "Aligned.out.21.sorted.bam")

bam_objects_list <- BamFileList(bam_paths, index = index_paths)

### LOOP
for (target in target_vector) {
  
  #GRange for target of interest
  to_subset <- subset(x = mygtf, 
                      mygtf$gene_id == target & mygtf$type %in% "gene")
  
  #GRange for targetsite of target of interest
  targetsite <- subsetByOverlaps(target_GR, to_subset)
  targetsite <- targetsite[1,]
  
  #find cleavage site of target of interest 
  target_midpoint <- targetsite %>% 
    promoters(., downstream = 11) %>% 
    resize(., width = 1, fix = "end")
  
  target_midpoint_df <- target_midpoint %>% as.data.frame()
  
  ### create GRanges of bam file on target
  what <- c("seq", "cigar")
  which <- GRanges(to_subset)
  param <- ScanBamParam(which=which, what=what)
  
  test <- lapply(bam_objects_list, function(x) {readGAlignments(x, param = param)})
  target_GRanges_list <- lapply(test, function(x) {as(x, "GRanges")})
  
  # rename GRange lists 
  file_names <- list.files("/binf-isilon/PBgrp/xpj980/datascratch/triplemutants/04.STAR_mappings/old_style_mapping/by_size/bam/", full.names = T, pattern = "21.sorted.bam")
  names(target_GRanges_list)
  good_names <- file_names %>%
    str_remove(pattern = ".*by_size/bam//") %>%
    str_remove(pattern = ".fastq.gz.trimmed.fq.gzAligned.out.21.sorted.bam") %>%
    str_replace_all(pattern = "-", replacement = "\\.")
  
  names(target_GRanges_list) <- good_names
  
  ## make a column with genotype and a column with replicate 
  new_list <- lapply(seq_along(target_GRanges_list), function(x) {
    test_name <- names(target_GRanges_list[x])
    test <- target_GRanges_list[[x]]
    test$sample <- test_name
    df <- test %>% as.data.frame() %>% tidyr::separate(sample, c("genotype", "replicate"), sep="_") %>%
      mutate(target = rep(target))
    return(df)
  })
  
  assign(paste0("gathered_df_",target), purrr::reduce(new_list,.f = bind_rows))
}



gathered_gr_AT1G62910 <- as(gathered_df_AT1G62910, "GRanges")
gathered_gr_AT2G34710 <- as(gathered_df_AT2G34710, "GRanges")
gathered_gr_AT2G39675 <- as(gathered_df_AT2G39675, "GRanges")  
gathered_gr_AT4G18390 <- as(gathered_df_AT4G18390, "GRanges")
gathered_gr_AT1G48410 <- as(gathered_df_AT1G48410, "GRanges")
gathered_gr_AT3G26810 <- as(gathered_df_AT3G26810, "GRanges")
gathered_gr_AT1G63130 <- as(gathered_df_AT1G63130, "GRanges")
gathered_gr_AT1G50055 <- as(gathered_df_AT1G50055, "GRanges")

for (target in target_vector) {
  
  to_subset <- subset(x = mygtf, 
                      mygtf$gene_id == target & mygtf$type %in% "gene")
  
  #GRange for targetsite of target of interest
  targetsite <- subsetByOverlaps(target_GR, to_subset)
  targetsite <- targetsite[1,]
  
  #find cleavage site of target of interest 
  assign(paste0("target_midpoint",target),targetsite %>% 
           promoters(., downstream = 11) %>% 
           resize(., width = 1, fix = "end"))
  
  assign(paste0("target_midpoint_df", target), targetsite %>% 
           promoters(., downstream = 11) %>% 
           resize(., width = 1, fix = "end") %>% as_data_frame())
}

to_subset_m <- subset(x = mygtf, 
                      mygtf$gene_id == "AT1G62910" & mygtf$type %in% "gene") 

#GRange for targetsite of target of interest
targetsite <- subsetByOverlaps(target_GR, to_subset_m)
targetsite <- targetsite[1,]

test <- targetsite %>% 
  promoters(., downstream = 11) %>% 
  resize(., width = 1, fix = "end") 

##### make phasing columns ####
gathered_gr_AT1G62910$phase <- 1 + (start(gathered_gr_AT1G62910) - start(test)-1) %% 21
gathered_gr_AT1G62910$side <- case_when(start(gathered_gr_AT1G62910) < start(test) ~ "left", T ~ "right")
if (target_midpoint_dfAT1G62910$strand == "-") {gathered_gr_AT1G62910$CF <- case_when(gathered_gr_AT1G62910$side =="right" ~ "5CF", T ~ "3CF")} else {gathered_gr_AT1G62910$CF <- case_when(gathered_gr_AT1G62910$side =="right" ~ "3CF", T ~ "5CF")}
gathered_AT1G62910 <- as.data.frame(gathered_gr_AT1G62910)

gathered_gr_AT2G34710$phase <- 1 + (start(gathered_gr_AT2G34710) - start(target_midpointAT2G34710)-1) %% 21
gathered_gr_AT2G34710$side <- case_when(start(gathered_gr_AT2G34710) < start(target_midpointAT2G34710) ~ "left", T ~ "right")
if (target_midpoint_dfAT2G34710$strand == "-") {gathered_gr_AT2G34710$CF <- case_when(gathered_gr_AT2G34710$side =="right" ~ "5CF", T ~ "3CF")} else {gathered_gr_AT2G34710$CF <- case_when(gathered_gr_AT2G34710$side =="right" ~ "3CF", T ~ "5CF")}
gathered_AT2G34710 <- as.data.frame(gathered_gr_AT2G34710)

gathered_gr_AT2G39675$phase <- 1 + (start(gathered_gr_AT2G39675) - start(target_midpointAT2G39675)-1) %% 21
gathered_gr_AT2G39675$side <- case_when(start(gathered_gr_AT2G39675) < start(target_midpointAT2G39675) ~ "left", T ~ "right")
if (target_midpoint_dfAT2G39675$strand == "-") {gathered_gr_AT2G39675$CF <- case_when(gathered_gr_AT2G39675$side =="right" ~ "5CF", T ~ "3CF")} else {gathered_gr_AT2G39675$CF <- case_when(gathered_gr_AT2G39675$side =="right" ~ "3CF", T ~ "5CF")}
gathered_AT2G39675 <- data.frame(gathered_gr_AT2G39675)

gathered_gr_AT4G18390$phase <- 1 + (start(gathered_gr_AT4G18390) - start(target_midpointAT4G18390)-1) %% 21
gathered_gr_AT4G18390$side <- case_when(start(gathered_gr_AT4G18390) < start(target_midpointAT4G18390) ~ "left", T ~ "right")
if (target_midpoint_dfAT4G18390$strand == "-") {gathered_gr_AT4G18390$CF <- case_when(gathered_gr_AT4G18390$side =="right" ~ "5CF", T ~ "3CF")} else {gathered_gr_AT4G18390$CF <- case_when(gathered_gr_AT4G18390$side =="right" ~ "3CF", T ~ "5CF")}
gathered_AT4G18390 <- data.frame(gathered_gr_AT4G18390)

gathered_gr_AT3G26810$phase <- 1 + (start(gathered_gr_AT3G26810) - start(target_midpointAT3G26810)-1) %% 21
gathered_gr_AT3G26810$side <- case_when(start(gathered_gr_AT3G26810) < start(target_midpointAT3G26810) ~ "left", T ~ "right")
if (target_midpoint_dfAT3G26810$strand == "-") {gathered_gr_AT3G26810$CF <- case_when(gathered_gr_AT3G26810$side =="right" ~ "5CF", T ~ "3CF")} else {gathered_gr_AT3G26810$CF <- case_when(gathered_gr_AT3G26810$side =="right" ~ "3CF", T ~ "5CF")}
gathered_AT3G26810 <- data.frame(gathered_gr_AT3G26810)

gathered_gr_AT1G63130$phase <- 1 + (start(gathered_gr_AT1G63130) - start(target_midpointAT1G63130)-1) %% 21
gathered_gr_AT1G63130$side <- case_when(start(gathered_gr_AT1G63130) < start(target_midpointAT1G63130) ~ "left", T ~ "right")
if (target_midpoint_dfAT1G63130$strand == "-") {gathered_gr_AT1G63130$CF <- case_when(gathered_gr_AT1G63130$side =="right" ~ "5CF", T ~ "3CF")} else {gathered_gr_AT1G63130$CF <- case_when(gathered_gr_AT1G63130$side =="right" ~ "3CF", T ~ "5CF")}
gathered_AT1G63130 <- data.frame(gathered_gr_AT1G63130)

gathered_gr_AT1G50055$phase <- 1 + (start(gathered_gr_AT1G50055) - start(target_midpointAT1G50055)-1) %% 21
gathered_gr_AT1G50055$side <- case_when(start(gathered_gr_AT1G50055) < start(target_midpointAT1G50055) ~ "left", T ~ "right")
if (target_midpoint_dfAT1G50055$strand == "-") {gathered_gr_AT1G50055$CF <- case_when(gathered_gr_AT1G50055$side =="right" ~ "5CF", T ~ "3CF")} else {gathered_gr_AT1G50055$CF <- case_when(gathered_gr_AT1G50055$side =="right" ~ "3CF", T ~ "5CF")}
gathered_AT1G50055 <- data.frame(gathered_gr_AT1G50055)

#### make a dataframe to plot & count phases ####
dt_plot_AT1G62910 <- gathered_AT1G62910 %>%
  filter(cigar=="21M") %>%
  group_by(phase, strand, side, CF, genotype, replicate, target) %>%
  summarise(count=n()) %>%
  ungroup()

to_plot_AT1G62910 <- if (dt_plot_AT1G62910[which.max(dt_plot_AT1G62910$count),]$CF == "3CF") {filter(dt_plot_AT1G62910, CF %in% "3CF")} else {filter(dt_plot_AT1G62910, CF %in% "5CF")}

summed_to_plot_AT1G62910 <- to_plot_AT1G62910 %>% 
  filter(genotype %in% genotypes_sub) %>%
  group_by(genotype, replicate) %>% 
  summarise(sum = sum(count)) %>%
  #summarise(mean_sum = mean(sum)) %>%
  mutate(target = rep("AT1G62910"))

dt_plot_AT2G34710 <- gathered_AT2G34710 %>%
  filter(cigar=="21M") %>%
  group_by(phase, strand, side, CF, genotype, replicate, target) %>%
  summarise(count=n()) %>%
  ungroup()

to_plot_AT2G34710 <- if (dt_plot_AT2G34710[which.max(dt_plot_AT2G34710$count),]$CF == "3CF") {filter(dt_plot_AT2G34710, CF %in% "3CF")} else {filter(dt_plot_AT2G34710, CF %in% "5CF")}

summed_to_plot_AT2G34710 <- to_plot_AT2G34710 %>% 
  filter(genotype %in% genotypes_sub) %>%
  group_by(genotype, replicate) %>% 
  summarise(sum = sum(count)) %>% 
  #summarise(mean_sum = mean(sum)) %>%
  mutate(target = rep("AT2G34710"))

dt_plot_AT2G39675 <- gathered_AT2G39675 %>%
  filter(cigar=="21M") %>%
  group_by(phase, strand, side, CF, genotype, replicate, target) %>%
  summarise(count=n()) %>%
  ungroup()

to_plot_AT2G39675 <- if (dt_plot_AT2G39675[which.max(dt_plot_AT2G39675$count),]$CF == "3CF") {filter(dt_plot_AT2G39675, CF %in% "3CF")} else {filter(dt_plot_AT2G39675, CF %in% "5CF")}

summed_to_plot_AT2G39675 <- to_plot_AT2G39675 %>% 
  filter(genotype %in% genotypes_sub) %>%
  group_by(genotype, replicate) %>% 
  summarise(sum = sum(count)) %>% 
  #summarise(mean_sum = mean(sum)) %>%
  mutate(target = rep("AT2G39675"))

dt_plot_AT4G18390 <- gathered_AT4G18390 %>%
  filter(cigar=="21M") %>%
  group_by(phase, strand, side, CF, genotype, replicate, target) %>%
  summarise(count=n()) %>%
  ungroup()

to_plot_AT4G18390 <- if (dt_plot_AT4G18390[which.max(dt_plot_AT4G18390$count),]$CF == "3CF") {filter(dt_plot_AT4G18390, CF %in% "3CF")} else {filter(dt_plot_AT4G18390, CF %in% "5CF")}

summed_to_plot_AT4G18390 <- to_plot_AT4G18390 %>% 
  filter(genotype %in% genotypes_sub) %>%
  group_by(genotype, replicate) %>% 
  summarise(sum = sum(count)) %>% 
  #summarise(mean_sum = mean(sum)) %>%
  mutate(target = rep("AT4G18390"))

# AFB2

dt_plot_AT3G26810 <- gathered_AT3G26810 %>%
  filter(cigar=="21M") %>%
  group_by(phase, strand, side, CF, genotype, replicate, target) %>%
  summarise(count=n()) %>%
  ungroup()

to_plot_AT3G26810 <- if (dt_plot_AT3G26810[which.max(dt_plot_AT3G26810$count),]$CF == "3CF") {filter(dt_plot_AT3G26810, CF %in% "3CF")} else {filter(dt_plot_AT3G26810, CF %in% "5CF")}

summed_to_plot_AT3G26810 <- to_plot_AT3G26810 %>% 
  filter(genotype %in% genotypes_sub) %>%
  group_by(genotype, replicate) %>% 
  summarise(sum = sum(count)) %>% 
  #summarise(mean_sum = mean(sum)) %>%
  mutate(target = rep("AT3G26810"))

# RPF6

dt_plot_AT1G63130 <- gathered_AT1G63130 %>%
  filter(cigar=="21M") %>%
  group_by(phase, strand, side, CF, genotype, replicate, target) %>%
  summarise(count=n()) %>%
  ungroup()

to_plot_AT1G63130 <- if (dt_plot_AT1G63130[which.max(dt_plot_AT1G63130$count),]$CF == "3CF") {filter(dt_plot_AT1G63130, CF %in% "3CF")} else {filter(dt_plot_AT1G63130, CF %in% "5CF")}

summed_to_plot_AT1G63130 <- to_plot_AT1G63130 %>% 
  filter(genotype %in% genotypes_sub) %>%
  group_by(genotype, replicate) %>% 
  summarise(sum = sum(count)) %>% 
  #summarise(mean_sum = mean(sum)) %>%
  mutate(target = rep("AT1G63130"))

# TAS1B

dt_plot_AT1G50055 <- gathered_AT1G50055 %>%
  filter(cigar=="21M") %>%
  group_by(phase, strand, side, CF, genotype, replicate, target) %>%
  summarise(count=n()) %>%
  ungroup()

to_plot_AT1G50055 <- if (dt_plot_AT1G50055[which.max(dt_plot_AT1G50055$count),]$CF == "3CF") {filter(dt_plot_AT1G50055, CF %in% "3CF")} else {filter(dt_plot_AT1G50055, CF %in% "5CF")}

summed_to_plot_AT1G50055 <- to_plot_AT1G50055 %>% 
  filter(genotype %in% genotypes_sub) %>%
  group_by(genotype, replicate) %>% 
  summarise(sum = sum(count)) %>% 
  #summarise(mean_sum = mean(sum)) %>%
  mutate(target = rep("AT1G50055"))


all_gathered_df <- bind_rows(to_plot_AT1G62910, to_plot_AT2G34710, to_plot_AT2G39675, to_plot_AT4G18390, to_plot_AT3G26810, to_plot_AT1G63130, to_plot_AT1G50055)
all_gathered_sum_reads <- bind_rows(summed_to_plot_AT1G62910, summed_to_plot_AT2G34710, summed_to_plot_AT2G39675, summed_to_plot_AT4G18390, summed_to_plot_AT3G26810, summed_to_plot_AT1G63130, summed_to_plot_AT1G50055)
all_gathered_df_sub <-subset(all_gathered_df, genotype %in% genotypes_sub) 

all_gathered_df_sub$genotype <- factor(all_gathered_df_sub$genotype, levels = c("wt", "ski2.5", "rrp4", "urt1.1", "GKheso1", "heso1urt1"))
all_gathered_sum_reads$genotype <- factor(all_gathered_sum_reads$genotype, levels = c("wt", "ski2.5", "rrp4", "urt1.1", "GKheso1", "heso1urt1"))
all_gathered_df_sub$target <- factor(all_gathered_df_sub$target, levels = c("AT1G62910","AT2G39675", "AT2G34710", "AT4G18390", "AT3G26810", "AT1G63130", "AT1G50055"))
all_gathered_sum_reads$target <- factor(all_gathered_sum_reads$target, levels = c("AT1G62910","AT2G39675", "AT2G34710", "AT4G18390", "AT3G26810",  "AT1G63130", "AT1G50055")) 


read_intensity_plot <- all_gathered_sum_reads %>% 
  filter(., genotype != "wt") %>%
  ggplot(., aes(x=genotype, y=sum, fill = target, alpha = 0.5)) +
  stat_summary(geom = "bar", fun = "mean", position = "dodge", aes(width=0.5)) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position = "dodge", size = 0.5, aes(width=0.5))+
  #geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~target, scales = "free", nrow = 4) +
  scale_fill_brewer(palette = "Dark2") + 
  cowplot::theme_cowplot()
ggsave(read_intensity_plot, filename = paste0("/isdata/PBgrp/xpj980/figures/readsbarplotURT1phas201223.svg"), units = "cm", width = 18, height = 18)


to_plot_AT1G62910$genotype <- factor(to_plot_AT1G62910$genotype, levels = c("Col.0.WT", "ski2.5", "hen2.5", "rrp4"))
to_plot_AT1G62910$replicate <- factor(to_plot_AT1G62910$replicate, levels = c("R3", "R2", "R1"))
together <- to_plot_AT1G62910 %>% 
  filter(genotype %in% genotypes_sub) %>% 
  group_by(genotype, replicate, target) %>%
  mutate(ncount = count/sum(count)) %>%
  ggplot(aes(x=factor(phase, levels = 1:21), y=ncount)) +
  geom_bar(stat="identity", position= "identity", aes(group=replicate, fill=replicate), col="black", alpha=.8) + 
  coord_polar() +
  facet_grid(target~genotype, switch="both") + 
  #ggtitle(label = my_title) +
  theme(panel.grid = element_line(
    size = 0.25, linetype = 'solid', colour = "grey"),
    strip.background = element_rect(color="black"), 
    axis.title.x = element_text(margin = margin(t = 15))) +
  scale_fill_brewer(palette = rev("Reds")) +
  scale_y_continuous(limits = c(0,1))


# plot with wt, ski2, hen2, rrp4 
plot <- all_gathered_df_sub %>% 
  group_by(genotype, replicate, target) %>%
  mutate(ncount = count/sum(count)) %>%
  ggplot(aes(x=factor(phase, levels = 1:21), y=ncount)) +
  geom_bar(stat="identity", position= "identity", aes(group=target, fill=target), col="black", alpha=.8) + 
  coord_polar() +
  facet_grid(target~genotype, switch="both") + 
  #ggtitle(label = my_title) +
  theme(panel.grid = element_line(
    size = 0.25, linetype = 'solid', colour = "grey"),
    strip.background = element_rect(color="black"), 
    axis.title.x = element_text(margin = margin(t = 15))) +
  scale_fill_brewer(palette = "Dark2") +
  scale_y_continuous(limits = c(0,1))

ggsave(plot, filename = paste0("/isdata/PBgrp/xpj980/figures/phasingURT1AFB2.svg"), units = "cm", width = 18, height = 18)

### AT1G62910 ####

### target site
to_subset_m <- subset(x = mygtf, 
                      mygtf$gene_id == "AT1G62910" & mygtf$type %in% "gene") 

#GRange for targetsite of target of interest
targetsite <- subsetByOverlaps(target_GR, to_subset_m)
targetsite <- targetsite[1,]

test <- targetsite %>% 
  promoters(., downstream = 9) %>% 
  resize(., width = 1, fix = "end") 

### make phasing columns ###
gathered_gr_AT1G62910$phase <- 1 + (start(gathered_gr_AT1G62910) - start(test)-1) %% 21
gathered_gr_AT1G62910$side <- case_when(start(gathered_gr_AT1G62910) < start(test) ~ "left", T ~ "right")
if (target_midpoint_dfAT1G62910$strand == "-") {gathered_gr_AT1G62910$CF <- case_when(gathered_gr_AT1G62910$side =="right" ~ "5CF", T ~ "3CF")} else {gathered_gr_AT1G62910$CF <- case_when(gathered_gr_AT1G62910$side =="right" ~ "3CF", T ~ "5CF")}
gathered_AT1G62910 <- as.data.frame(gathered_gr_AT1G62910)

## make a dataframe to plot & count phases ##
dt_plot_AT1G62910 <- gathered_AT1G62910 %>%
  filter(cigar=="21M") %>%
  group_by(phase, strand, side, CF, genotype, replicate, target) %>%
  summarise(count=n()) %>%
  filter(genotype %in% genotypes_sub) %>% 
  ungroup()

to_plot_AT1G62910 <- if (dt_plot_AT1G62910[which.max(dt_plot_AT1G62910$count),]$CF == "3CF") {filter(dt_plot_AT1G62910, CF %in% "3CF")} else {filter(dt_plot_AT1G62910, CF %in% "5CF")}

summed_to_plot_AT1G62910 <- to_plot_AT1G62910 %>% 
  filter(genotype %in% genotypes_sub) %>%
  group_by(genotype, replicate) %>% 
  summarise(sum = sum(count)) %>%
  #summarise(mean_sum = mean(sum)) %>%
  mutate(target = rep("AT1G62910"))

to_plot_AT1G62910$genotype <- factor(to_plot_AT1G62910$genotype, levels = c("wt", "ski2.5", "rrp4", "GKheso1", "urt1.1"))
to_plot_AT1G62910$replicate <- factor(to_plot_AT1G62910$replicate, levels = c("R3", "R2", "R1"))
together <- to_plot_AT1G62910 %>% 
  group_by(phase, genotype, strand, target) %>%
  summarise(meancount = mean(count)) %>% 
  group_by(genotype, strand, target) %>%
  mutate(ncount = meancount/sum(meancount))

test <- together %>% 
  group_by(genotype, strand) %>% 
  #mutate(avg = sum(ncount)/3)
  summarise(sum = sum(ncount))

Rgene <- together %>%
  ggplot(aes(x=factor(phase, levels = 1:21), y=ncount)) +
  geom_bar(stat="identity", position= "identity", aes(fill=strand), col="black", size=0.25, alpha=.9) +
  geom_hline(yintercept = seq(0, 4, by = 1), colour = "darkgrey", size = 0.1) +
  geom_hline(yintercept = seq(0, 4, by = 0.25), colour = "darkgrey", size = 0.1) +
  geom_hline(yintercept = seq(0, 4, by = 0.50), colour = "darkgrey", size = 0.1) +
  geom_hline(yintercept = seq(0, 4, by = 0.75), colour = "darkgrey", size = 0.1) +
  geom_vline(xintercept = 1:21, colour = "darkgrey", size = 0.1) + 
  coord_polar() +
  facet_grid(target~genotype, switch="both") + 
  theme(panel.grid = element_blank()) +
  #ggtitle(label = my_title) +
  #theme(panel.grid = element_line(
  #size = 0.25, linetype = 'solid', colour = "grey"),
  #strip.background = element_rect(color="black"), 
  #axis.title.x = element_text(margin = margin(t = 15))) +
  scale_fill_manual(values = c("#AF2D2D", "#f05454")) +
  scale_y_continuous(limits = c(0,1))
ggsave(Rgene, filename = paste0("/isdata/PBgrp/xpj980/figures/RgeneavgURT1phasing.svg"), units = "cm", width = 18, height = 18)

summed_to_plot_AT1G62910$genotype <- factor(summed_to_plot_AT1G62910$genotype, levels = c("wt", "ski2.5", "rrp4", "GKheso1", "urt1.1"))

read_intensity_plot <- summed_to_plot_AT1G62910 %>% 
  filter(target %in% "AT1G62910") %>%
  ggplot(., aes(x=genotype, y=sum, fill = genotype)) + 
  stat_summary(geom = "bar", fun = "mean", position = "dodge", aes(width=0.5)) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position = "dodge", size = 0.5, aes(width=0.5))+
  #geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~target, scales = "free", nrow = 4) +
  scale_fill_brewer(palette = "Reds") + 
  cowplot::theme_cowplot()
ggsave(read_intensity_plot, filename = paste0("/isdata/PBgrp/xpj980/figures/RgenereadsavgURT1.svg"), units = "cm", width = 18, height = 18)


#### Tas1C ####

### target site
to_subset_m <- subset(x = mygtf, 
                      mygtf$gene_id == "AT2G39675" & mygtf$type %in% "gene") 

#GRange for targetsite of target of interest
targetsite <- subsetByOverlaps(target_GR, to_subset_m)
targetsite <- targetsite[1,]

test <- targetsite %>% 
  promoters(., downstream = 13) %>% 
  resize(., width = 1, fix = "end") 

gathered_gr_AT2G39675$phase <- 1 + (start(gathered_gr_AT2G39675) - start(test)-1) %% 21
gathered_gr_AT2G39675$side <- case_when(start(gathered_gr_AT2G39675) < start(test) ~ "left", T ~ "right")
if (target_midpoint_dfAT2G39675$strand == "-") {gathered_gr_AT2G39675$CF <- case_when(gathered_gr_AT2G39675$side =="right" ~ "5CF", T ~ "3CF")} else {gathered_gr_AT2G39675$CF <- case_when(gathered_gr_AT2G39675$side =="right" ~ "3CF", T ~ "5CF")}
gathered_AT2G39675 <- data.frame(gathered_gr_AT2G39675)

dt_plot_AT2G39675 <- gathered_AT2G39675 %>%
  filter(cigar=="21M") %>%
  group_by(phase, strand, side, CF, genotype, replicate, target) %>%
  summarise(count=n()) %>%
  filter(genotype %in% genotypes_sub) %>% 
  ungroup()

to_plot_AT2G39675 <- if (dt_plot_AT2G39675[which.max(dt_plot_AT2G39675$count),]$CF == "3CF") {filter(dt_plot_AT2G39675, CF %in% "3CF")} else {filter(dt_plot_AT2G39675, CF %in% "5CF")}

summed_to_plot_AT2G39675 <- to_plot_AT2G39675 %>% 
  filter(genotype %in% genotypes_sub) %>%
  group_by(genotype, replicate) %>% 
  summarise(sum = sum(count)) %>% 
  #summarise(mean_sum = mean(sum)) %>%
  mutate(target = rep("AT2G39675"))

to_plot_AT2G39675$genotype <- factor(to_plot_AT2G39675$genotype, levels = c("wt", "ski2.5", "rrp4", "urt1.1", "GKheso1", "heso1urt1"))
to_plot_AT2G39675$replicate <- factor(to_plot_AT2G39675$replicate, levels = c("R3", "R2", "R1"))

together <- to_plot_AT2G39675 %>% 
  group_by(phase, genotype, target, strand) %>%
  summarise(meancount = mean(count)) %>%
  group_by(genotype, target, strand) %>%
  mutate(ncount = meancount/sum(meancount)) 

tas1 <- together %>% 
  ggplot(aes(x=factor(phase, levels = 1:21), y=ncount)) +
  geom_bar(stat="identity", position= "identity", aes(fill=strand), col="black", size=0.25, alpha=.9) +
  geom_hline(yintercept = seq(0, 4, by = 1), colour = "darkgrey", size = 0.1) +
  geom_hline(yintercept = seq(0, 4, by = 0.25), colour = "darkgrey", size = 0.1) +
  geom_hline(yintercept = seq(0, 4, by = 0.50), colour = "darkgrey", size = 0.1) +
  geom_hline(yintercept = seq(0, 4, by = 0.75), colour = "darkgrey", size = 0.1) +
  geom_vline(xintercept = 1:21, colour = "darkgrey", size = 0.1) + 
  coord_polar() +
  facet_grid(target~genotype, switch="both") + 
  theme(panel.grid = element_blank()) +
  #ggtitle(label = my_title) +
  #theme(panel.grid = element_line(
  #size = 0.25, linetype = 'solid', colour = "grey"),
  #strip.background = element_rect(color="black"), 
  #axis.title.x = element_text(margin = margin(t = 15))) +
  scale_fill_manual(values = c("#7ea04d", "#335d2d")) +
  scale_y_continuous(limits = c(0,1))
ggsave(tas1, filename = paste0("/isdata/PBgrp/xpj980/figures/Tas1geneavgURT1201223.svg"), units = "cm", width = 18, height = 18)

summed_to_plot_AT2G39675$genotype <- factor(summed_to_plot_AT2G39675$genotype, levels = c("wt", "ski2.5", "rrp4", "urt1.1",  "GKheso1", "heso1urt1"))

read_intensity_plot <- summed_to_plot_AT2G39675 %>% 
  filter(target %in% "AT2G39675") %>%
  ggplot(., aes(x=genotype, y=sum, fill = genotype)) + 
  stat_summary(geom = "bar", fun = "mean", position = "dodge", aes(width=0.5)) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position = "dodge", size = 0.5, aes(width=0.5))+
  #geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~target, scales = "free", nrow = 4) +
  scale_fill_brewer(palette = "Greens") + 
  cowplot::theme_cowplot()
ggsave(read_intensity_plot, filename = paste0("/isdata/PBgrp/xpj980/figures/Tas1genereadsavgURT1201223.svg"), units = "cm", width = 18, height = 18)


#### PHB #### 
to_subset_m <- subset(x = mygtf, 
                      mygtf$gene_id == "AT2G34710" & mygtf$type %in% "gene") 

#GRange for targetsite of target of interest
targetsite <- subsetByOverlaps(target_GR, to_subset_m)
targetsite <- targetsite[1,]

test <- targetsite %>% 
  promoters(., downstream = 11) %>% 
  resize(., width = 1, fix = "end") 

gathered_gr_AT2G34710$phase <- 1 + (start(gathered_gr_AT2G34710) - start(test)-1) %% 21
gathered_gr_AT2G34710$side <- case_when(start(gathered_gr_AT2G34710) < start(test) ~ "left", T ~ "right")
if (target_midpoint_dfAT2G34710$strand == "-") {gathered_gr_AT2G34710$CF <- case_when(gathered_gr_AT2G34710$side =="right" ~ "5CF", T ~ "3CF")} else {gathered_gr_AT2G34710$CF <- case_when(gathered_gr_AT2G34710$side =="right" ~ "3CF", T ~ "5CF")}
gathered_AT2G34710 <- as.data.frame(gathered_gr_AT2G34710)

dt_plot_AT2G34710 <- gathered_AT2G34710 %>%
  filter(cigar=="21M") %>%
  group_by(phase, strand, side, CF, genotype, replicate, target) %>%
  summarise(count=n()) %>%
  filter(genotype %in% genotypes_sub) %>%
  ungroup()

to_plot_AT2G34710 <- if (dt_plot_AT2G34710[which.max(dt_plot_AT2G34710$count),]$CF == "3CF") {filter(dt_plot_AT2G34710, CF %in% "3CF")} else {filter(dt_plot_AT2G34710, CF %in% "5CF")}

summed_to_plot_AT2G34710 <- to_plot_AT2G34710 %>% 
  filter(genotype %in% genotypes_sub) %>%
  group_by(genotype, replicate) %>% 
  summarise(sum = sum(count)) %>% 
  #summarise(mean_sum = mean(sum)) %>%
  mutate(target = rep("AT2G34710"))

to_plot_AT2G34710$genotype <- factor(to_plot_AT2G34710$genotype, levels = c("wt", "ski2.5", "rrp4", "GKheso1", "urt1.1"))
to_plot_AT2G34710$replicate <- factor(to_plot_AT2G34710$replicate, levels = c("R3", "R2", "R1"))
together <- to_plot_AT2G34710 %>% 
  group_by(phase, genotype, target, strand) %>%
  summarise(meancount = mean(count)) %>%
  group_by(genotype, target, strand) %>%
  mutate(ncount = meancount/sum(meancount)) 

test <- together %>% 
  group_by(genotype, strand) %>% 
  #mutate(avg = sum(ncount)/3)
  summarise(sum = sum(ncount))

PHB <- together %>% 
  group_by(genotype, target, strand) %>%
  ggplot(aes(x=factor(phase, levels = 1:21), y=ncount)) +
  #stat_summary(geom = "bar", fun = "mean", position = "dodge", aes(fill=strand)) +
  geom_bar(stat="identity", position= "identity", aes(fill=strand), col="black", size=0.25, alpha=.9) +
  geom_hline(yintercept = seq(0, 4, by = 1), colour = "darkgrey", size = 0.1) +
  geom_hline(yintercept = seq(0, 4, by = 0.25), colour = "darkgrey", size = 0.1) +
  geom_hline(yintercept = seq(0, 4, by = 0.50), colour = "darkgrey", size = 0.1) +
  geom_hline(yintercept = seq(0, 4, by = 0.75), colour = "darkgrey", size = 0.1) +
  geom_vline(xintercept = 1:21, colour = "darkgrey", size = 0.1) + 
  coord_polar() +
  facet_grid(target~genotype, switch="both") + 
  theme(panel.grid = element_blank()) +
  #ggtitle(label = my_title) +
  #theme(panel.grid = element_line(
  #size = 0.25, linetype = 'solid', colour = "grey"),
  #strip.background = element_rect(color="black"), 
  #axis.title.x = element_text(margin = margin(t = 15))) +
  scale_fill_manual(values = c("#31112c", "#7d0633")) +
  scale_y_continuous(limits = c(0,1))
ggsave(PHB, filename = paste0("/isdata/PBgrp/xpj980/figures/PHBavgURT1.svg"), units = "cm", width = 18, height = 18)

summed_to_plot_AT2G34710$genotype <- factor(summed_to_plot_AT2G34710$genotype, levels = c("wt", "ski2.5", "rrp4", "GKheso1", "urt1.1"))

read_intensity_plot <- summed_to_plot_AT2G34710 %>% 
  filter(target %in% "AT2G34710") %>%
  ggplot(., aes(x=genotype, y=sum, fill = genotype)) + 
  stat_summary(geom = "bar", fun = "mean", position = "dodge", aes(width=0.5)) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position = "dodge", size = 0.5, aes(width=0.5))+
  #geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~target, scales = "free", nrow = 4) +
  scale_fill_brewer(palette = "Purples") + 
  cowplot::theme_cowplot()
ggsave(read_intensity_plot, filename = paste0("/isdata/PBgrp/xpj980/figures/PHBreadsavgURT1.svg"), units = "cm", width = 18, height = 18)

#### TCP2 ####

to_subset_m <- subset(x = mygtf, 
                      mygtf$gene_id == "AT4G18390" & mygtf$type %in% "gene") 

#GRange for targetsite of target of interest
targetsite <- subsetByOverlaps(target_GR, to_subset_m)
targetsite <- targetsite[1,]

test <- targetsite %>% 
  promoters(., downstream = 11) %>% 
  resize(., width = 1, fix = "end") 

gathered_gr_AT4G18390$phase <- 1 + (start(gathered_gr_AT4G18390) - start(test)-1) %% 21
gathered_gr_AT4G18390$side <- case_when(start(gathered_gr_AT4G18390) < start(test) ~ "left", T ~ "right")
if (target_midpoint_dfAT4G18390$strand == "-") {gathered_gr_AT4G18390$CF <- case_when(gathered_gr_AT4G18390$side =="right" ~ "5CF", T ~ "3CF")} else {gathered_gr_AT4G18390$CF <- case_when(gathered_gr_AT4G18390$side =="right" ~ "3CF", T ~ "5CF")}
gathered_AT4G18390 <- data.frame(gathered_gr_AT4G18390)

dt_plot_AT4G18390 <- gathered_AT4G18390 %>%
  filter(cigar=="21M") %>%
  group_by(phase, strand, side, CF, genotype, replicate, target) %>%
  summarise(count=n()) %>%
  filter(genotype %in% genotypes_sub) %>% 
  ungroup()

to_plot_AT4G18390 <- if (dt_plot_AT4G18390[which.max(dt_plot_AT4G18390$count),]$CF == "3CF") {filter(dt_plot_AT4G18390, CF %in% "3CF")} else {filter(dt_plot_AT4G18390, CF %in% "5CF")}

summed_to_plot_AT4G18390 <- to_plot_AT4G18390 %>% 
  filter(genotype %in% genotypes_sub) %>%
  group_by(genotype, replicate) %>% 
  summarise(sum = sum(count)) %>% 
  #summarise(mean_sum = mean(sum)) %>%
  mutate(target = rep("AT4G18390"))

to_plot_AT4G18390$genotype <- factor(to_plot_AT4G18390$genotype, levels = c("wt", "ski2.5", "rrp4", "GKheso1", "urt1.1"))
to_plot_AT4G18390$replicate <- factor(to_plot_AT4G18390$replicate, levels = c("R3", "R2", "R1"))
together <- to_plot_AT4G18390 %>% 
  group_by(phase, genotype, strand, target) %>%
  summarise(meancount = mean(count)) %>%
  group_by(genotype, strand, target) %>%
  mutate(ncount = meancount/sum(meancount)) 

test <- together %>% 
  group_by(genotype, strand) %>% 
  #mutate(avg = sum(ncount)/3)
  summarise(sum = sum(ncount))  

tcp2 <- together %>% 
  ggplot(aes(x=factor(phase, levels = 1:21), y=ncount)) +
  geom_bar(stat="identity", position= "identity", aes(fill=strand), col="black", size=0.25) +
  geom_hline(yintercept = seq(0, 4, by = 1), colour = "darkgrey", size = 0.1) +
  geom_hline(yintercept = seq(0, 4, by = 0.25), colour = "darkgrey", size = 0.1) +
  geom_hline(yintercept = seq(0, 4, by = 0.50), colour = "darkgrey", size = 0.1) +
  geom_hline(yintercept = seq(0, 4, by = 0.75), colour = "darkgrey", size = 0.1) +
  geom_vline(xintercept = 1:21, colour = "darkgrey", size = 0.1) + 
  coord_polar() +
  facet_grid(target~genotype, switch="both") + 
  theme(panel.grid = element_blank()) +
  #ggtitle(label = my_title) +
  #theme(panel.grid = element_line(
  #size = 0.25, linetype = 'solid', colour = "grey"),
  #strip.background = element_rect(color="black"), 
  #axis.title.x = element_text(margin = margin(t = 15))) +
  scale_fill_manual(values = c("#00334e", "#145374")) +
  scale_y_continuous(limits = c(0,1))
ggsave(tcp2, filename = paste0("/isdata/PBgrp/xpj980/figures/TCP2avgURT1.svg"), units = "cm", width = 18, height = 18)

summed_to_plot_AT4G18390$genotype <- factor(summed_to_plot_AT4G18390$genotype, levels = c("wt", "ski2.5", "rrp4", "GKheso1", "urt1.1"))

read_intensity_plot <- summed_to_plot_AT4G18390 %>% 
  filter(target %in% "AT4G18390") %>%
  ggplot(., aes(x=genotype, y=sum, fill = genotype)) + 
  stat_summary(geom = "bar", fun = "mean", position = "dodge", aes(width=0.5)) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position = "dodge", size = 0.5, aes(width=0.5))+
  #geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~target, scales = "free", nrow = 4) +
  scale_fill_brewer(palette = "Blues") + 
  cowplot::theme_cowplot()
ggsave(read_intensity_plot, filename = paste0("/isdata/PBgrp/xpj980/figures/TCP2readsavgURT1.svg"), units = "cm", width = 18, height = 18)

#### AGO1 ####

### target site
to_subset_m <- subset(x = mygtf, 
                      mygtf$gene_id == "AT1G48410" & mygtf$type %in% "gene") 

#GRange for targetsite of target of interest
targetsite <- subsetByOverlaps(target_GR, to_subset_m)
targetsite <- targetsite[1,]

test <- targetsite %>% 
  promoters(., downstream = 12) %>% 
  resize(., width = 1, fix = "end") 

gathered_gr_AT1G48410$phase <- 1 + (start(gathered_gr_AT1G48410) - start(test)-1) %% 21
gathered_gr_AT1G48410$side <- case_when(start(gathered_gr_AT1G48410) < start(test) ~ "left", T ~ "right")
if (target_midpoint_dfAT1G48410$strand == "-") {gathered_gr_AT1G48410$CF <- case_when(gathered_gr_AT1G48410$side =="right" ~ "5CF", T ~ "3CF")} else {gathered_gr_AT1G48410$CF <- case_when(gathered_gr_AT1G48410$side =="right" ~ "3CF", T ~ "5CF")}
gathered_AT1G48410 <- data.frame(gathered_gr_AT1G48410)

dt_plot_AT1G48410 <- gathered_AT1G48410 %>%
  filter(cigar=="21M") %>%
  group_by(phase, strand, side, CF, genotype, replicate, target) %>%
  summarise(count=n()) %>%
  filter(genotype %in% genotypes_sub) %>% 
  ungroup()

to_plot_AT1G48410 <- if (dt_plot_AT1G48410[which.max(dt_plot_AT1G48410$count),]$CF == "3CF") {filter(dt_plot_AT1G48410, CF %in% "3CF")} else {filter(dt_plot_AT1G48410, CF %in% "5CF")}

summed_to_plot_AT1G48410 <- to_plot_AT1G48410 %>% 
  filter(genotype %in% genotypes_sub) %>%
  group_by(genotype, replicate) %>% 
  summarise(sum = sum(count)) %>% 
  #summarise(mean_sum = mean(sum)) %>%
  mutate(target = rep("AT1G48410"))

to_plot_AT1G48410$genotype <- factor(to_plot_AT1G48410$genotype, levels = c("wt", "ski2.5", "rrp4", "urt1.1", "GKheso1", "heso1urt1"))
to_plot_AT1G48410$replicate <- factor(to_plot_AT1G48410$replicate, levels = c("R3", "R2", "R1"))

together <- to_plot_AT1G48410 %>% 
  group_by(phase, genotype, target, strand) %>%
  summarise(meancount = mean(count)) %>%
  group_by(genotype, target, strand) %>%
  mutate(ncount = meancount/sum(meancount)) 

ago1 <- together %>% 
  ggplot(aes(x=factor(phase, levels = 1:21), y=ncount)) +
  geom_bar(stat="identity", position= "identity", aes(fill=strand), col="black", size=0.25, alpha=.9) +
  geom_hline(yintercept = seq(0, 4, by = 1), colour = "darkgrey", size = 0.1) +
  geom_hline(yintercept = seq(0, 4, by = 0.25), colour = "darkgrey", size = 0.1) +
  geom_hline(yintercept = seq(0, 4, by = 0.50), colour = "darkgrey", size = 0.1) +
  geom_hline(yintercept = seq(0, 4, by = 0.75), colour = "darkgrey", size = 0.1) +
  geom_vline(xintercept = 1:21, colour = "darkgrey", size = 0.1) + 
  coord_polar() +
  facet_grid(target~genotype, switch="both") + 
  theme(panel.grid = element_blank()) +
  #ggtitle(label = my_title) +
  #theme(panel.grid = element_line(
  #size = 0.25, linetype = 'solid', colour = "grey"),
  #strip.background = element_rect(color="black"), 
  #axis.title.x = element_text(margin = margin(t = 15))) +
  scale_fill_manual(values = c("#fdb827", "#ffd66b")) +
  scale_y_continuous(limits = c(0,1))
ggsave(ago1, filename = paste0("/isdata/PBgrp/xpj980/figures/AGO1geneavgURT1201223.svg"), units = "cm", width = 18, height = 18)

summed_to_plot_AT1G48410$genotype <- factor(summed_to_plot_AT1G48410$genotype, levels = c("wt", "ski2.5", "rrp4", "urt1.1", "GKheso1", "heso1urt1"))

read_intensity_plot <- summed_to_plot_AT1G48410 %>% 
  filter(target %in% "AT1G48410") %>%
  ggplot(., aes(x=genotype, y=sum, fill = genotype)) + 
  stat_summary(geom = "bar", fun = "mean", position = "dodge", aes(width=0.5)) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position = "dodge", size = 0.5, aes(width=0.5))+
  #geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~target, scales = "free", nrow = 4) +
  scale_fill_brewer(palette = "YlOrRd") + 
  cowplot::theme_cowplot()
ggsave(read_intensity_plot, filename = paste0("/isdata/PBgrp/xpj980/figures/AGO1genereadsavgURT1201223.svg"), units = "cm", width = 18, height = 18)

#### AFB2 ####

### target site
to_subset_m <- subset(x = mygtf, 
                      mygtf$gene_id == "AT3G26810" & mygtf$type %in% "gene") 

#GRange for targetsite of target of interest
targetsite <- subsetByOverlaps(target_GR, to_subset_m)
targetsite <- targetsite[1,]

test <- targetsite %>% 
  promoters(., downstream = 12) %>% 
  resize(., width = 1, fix = "end") 

gathered_gr_AT3G26810$phase <- 1 + (start(gathered_gr_AT3G26810) - start(test)-1) %% 21
gathered_gr_AT3G26810$side <- case_when(start(gathered_gr_AT3G26810) < start(test) ~ "left", T ~ "right")
if (target_midpoint_dfAT3G26810$strand == "-") {gathered_gr_AT3G26810$CF <- case_when(gathered_gr_AT3G26810$side =="right" ~ "5CF", T ~ "3CF")} else {gathered_gr_AT3G26810$CF <- case_when(gathered_gr_AT3G26810$side =="right" ~ "3CF", T ~ "5CF")}
gathered_AT3G26810 <- data.frame(gathered_gr_AT3G26810)

dt_plot_AT3G26810 <- gathered_AT3G26810 %>%
  filter(cigar=="21M") %>%
  group_by(phase, strand, side, CF, genotype, replicate, target) %>%
  summarise(count=n()) %>%
  filter(genotype %in% genotypes_sub) %>% 
  ungroup()

to_plot_AT3G26810 <- if (dt_plot_AT3G26810[which.max(dt_plot_AT3G26810$count),]$CF == "3CF") {filter(dt_plot_AT3G26810, CF %in% "3CF")} else {filter(dt_plot_AT3G26810, CF %in% "5CF")}

summed_to_plot_AT3G26810 <- to_plot_AT3G26810 %>% 
  filter(genotype %in% genotypes_sub) %>%
  group_by(genotype, replicate) %>% 
  summarise(sum = sum(count)) %>% 
  #summarise(mean_sum = mean(sum)) %>%
  mutate(target = rep("AT3G26810"))

to_plot_AT3G26810$genotype <- factor(to_plot_AT3G26810$genotype, levels = c("wt", "ski2.5", "rrp4", "urt1.1", "GKheso1", "heso1urt1"))
to_plot_AT3G26810$replicate <- factor(to_plot_AT3G26810$replicate, levels = c("R3", "R2", "R1"))

together <- to_plot_AT3G26810 %>% 
  group_by(phase, genotype, target, strand) %>%
  summarise(meancount = mean(count)) %>%
  group_by(genotype, target, strand) %>%
  mutate(ncount = meancount/sum(meancount)) 

afb2 <- together %>% 
  ggplot(aes(x=factor(phase, levels = 1:21), y=ncount)) +
  geom_bar(stat="identity", position= "identity", aes(fill=strand), col="black", size=0.25, alpha=.9) +
  geom_hline(yintercept = seq(0, 4, by = 1), colour = "darkgrey", size = 0.1) +
  geom_hline(yintercept = seq(0, 4, by = 0.25), colour = "darkgrey", size = 0.1) +
  geom_hline(yintercept = seq(0, 4, by = 0.50), colour = "darkgrey", size = 0.1) +
  geom_hline(yintercept = seq(0, 4, by = 0.75), colour = "darkgrey", size = 0.1) +
  geom_vline(xintercept = 1:21, colour = "darkgrey", size = 0.1) + 
  coord_polar() +
  facet_grid(target~genotype, switch="both") + 
  theme(panel.grid = element_blank()) +
  #ggtitle(label = my_title) +
  #theme(panel.grid = element_line(
  #size = 0.25, linetype = 'solid', colour = "grey"),
  #strip.background = element_rect(color="black"), 
  #axis.title.x = element_text(margin = margin(t = 15))) +
  scale_fill_manual(values = c("#DA0C81", "#E95793")) +
  scale_y_continuous(limits = c(0,1))
ggsave(afb2, filename = paste0("/isdata/PBgrp/xpj980/figures/AFB2geneavgURT1201223.svg"), units = "cm", width = 18, height = 18)

summed_to_plot_AT3G26810$genotype <- factor(summed_to_plot_AT3G26810$genotype, levels = c("wt", "ski2.5", "rrp4", "urt1.1", "GKheso1", "heso1urt1"))

read_intensity_plot <- summed_to_plot_AT3G26810 %>% 
  filter(target %in% "AT3G26810") %>%
  ggplot(., aes(x=genotype, y=sum, fill = genotype)) + 
  stat_summary(geom = "bar", fun = "mean", position = "dodge", aes(width=0.5)) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position = "dodge", size = 0.5, aes(width=0.5))+
  #geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~target, scales = "free", nrow = 4) +
  scale_fill_brewer(palette = "YlOrRd") + 
  cowplot::theme_cowplot()
ggsave(read_intensity_plot, filename = paste0("/isdata/PBgrp/xpj980/figures/AFB2genereadsavgURT1201223.svg"), units = "cm", width = 18, height = 18)

#### RPf6 ####

### target site
to_subset_m <- subset(x = mygtf, 
                      mygtf$gene_id == "AT1G63130" & mygtf$type %in% "gene") 

#GRange for targetsite of target of interest
targetsite <- subsetByOverlaps(target_GR, to_subset_m)
targetsite <- targetsite[1,]

test <- targetsite %>% 
  promoters(., downstream = 12) %>% 
  resize(., width = 1, fix = "end") 

gathered_gr_AT1G63130$phase <- 1 + (start(gathered_gr_AT1G63130) - start(test)-1) %% 21
gathered_gr_AT1G63130$side <- case_when(start(gathered_gr_AT1G63130) < start(test) ~ "left", T ~ "right")
if (target_midpoint_dfAT1G63130$strand == "-") {gathered_gr_AT1G63130$CF <- case_when(gathered_gr_AT1G63130$side =="right" ~ "5CF", T ~ "3CF")} else {gathered_gr_AT1G63130$CF <- case_when(gathered_gr_AT1G63130$side =="right" ~ "3CF", T ~ "5CF")}
gathered_AT1G63130 <- data.frame(gathered_gr_AT1G63130)

dt_plot_AT1G63130 <- gathered_AT1G63130 %>%
  filter(cigar=="21M") %>%
  group_by(phase, strand, side, CF, genotype, replicate, target) %>%
  summarise(count=n()) %>%
  filter(genotype %in% genotypes_sub) %>% 
  ungroup()

to_plot_AT1G63130 <- if (dt_plot_AT1G63130[which.max(dt_plot_AT1G63130$count),]$CF == "3CF") {filter(dt_plot_AT1G63130, CF %in% "3CF")} else {filter(dt_plot_AT1G63130, CF %in% "5CF")}

summed_to_plot_AT1G63130 <- to_plot_AT1G63130 %>% 
  filter(genotype %in% genotypes_sub) %>%
  group_by(genotype, replicate) %>% 
  summarise(sum = sum(count)) %>% 
  #summarise(mean_sum = mean(sum)) %>%
  mutate(target = rep("AT1G63130"))

to_plot_AT1G63130$genotype <- factor(to_plot_AT1G63130$genotype, levels = c("wt", "ski2.5", "rrp4", "urt1.1", "GKheso1", "heso1urt1"))
to_plot_AT1G63130$replicate <- factor(to_plot_AT1G63130$replicate, levels = c("R3", "R2", "R1"))

together <- to_plot_AT1G63130 %>% 
  group_by(phase, genotype, target, strand) %>%
  summarise(meancount = mean(count)) %>%
  group_by(genotype, target, strand) %>%
  mutate(ncount = meancount/sum(meancount)) 

rpf6 <- together %>% 
  ggplot(aes(x=factor(phase, levels = 1:21), y=ncount)) +
  geom_bar(stat="identity", position= "identity", aes(fill=strand), col="black", size=0.25, alpha=.9) +
  geom_hline(yintercept = seq(0, 4, by = 1), colour = "darkgrey", size = 0.1) +
  geom_hline(yintercept = seq(0, 4, by = 0.25), colour = "darkgrey", size = 0.1) +
  geom_hline(yintercept = seq(0, 4, by = 0.50), colour = "darkgrey", size = 0.1) +
  geom_hline(yintercept = seq(0, 4, by = 0.75), colour = "darkgrey", size = 0.1) +
  geom_vline(xintercept = 1:21, colour = "darkgrey", size = 0.1) + 
  coord_polar() +
  facet_grid(target~genotype, switch="both") + 
  theme(panel.grid = element_blank()) +
  #ggtitle(label = my_title) +
  #theme(panel.grid = element_line(
  #size = 0.25, linetype = 'solid', colour = "grey"),
  #strip.background = element_rect(color="black"), 
  #axis.title.x = element_text(margin = margin(t = 15))) +
  scale_fill_manual(values = c("#3081D0", "#6DB9EF")) +
  scale_y_continuous(limits = c(0,1))
ggsave(rpf6, filename = paste0("/isdata/PBgrp/xpj980/figures/RPF6geneavgURT1201223.svg"), units = "cm", width = 18, height = 18)

summed_to_plot_AT1G63130$genotype <- factor(summed_to_plot_AT1G63130$genotype, levels = c("wt", "ski2.5", "rrp4", "urt1.1", "GKheso1", "heso1urt1"))

read_intensity_plot <- summed_to_plot_AT1G63130 %>% 
  filter(target %in% "AT1G63130") %>%
  ggplot(., aes(x=genotype, y=sum, fill = genotype)) + 
  stat_summary(geom = "bar", fun = "mean", position = "dodge", aes(width=0.5)) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position = "dodge", size = 0.5, aes(width=0.5))+
  #geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~target, scales = "free", nrow = 4) +
  scale_fill_brewer(palette = "YlOrRd") + 
  cowplot::theme_cowplot()
ggsave(read_intensity_plot, filename = paste0("/isdata/PBgrp/xpj980/figures/rpf6genereadsavgURT1201223.svg"), units = "cm", width = 18, height = 18)

#### TAS1B ####

### target site
to_subset_m <- subset(x = mygtf, 
                      mygtf$gene_id == "AT1G50055" & mygtf$type %in% "gene") 

#GRange for targetsite of target of interest
targetsite <- subsetByOverlaps(target_GR, to_subset_m)
targetsite <- targetsite[1,]

test <- targetsite %>% 
  promoters(., downstream = 12) %>% 
  resize(., width = 1, fix = "end") 

gathered_gr_AT1G50055$phase <- 1 + (start(gathered_gr_AT1G50055) - start(test)-1) %% 21
gathered_gr_AT1G50055$side <- case_when(start(gathered_gr_AT1G50055) < start(test) ~ "left", T ~ "right")
if (target_midpoint_dfAT1G50055$strand == "-") {gathered_gr_AT1G50055$CF <- case_when(gathered_gr_AT1G50055$side =="right" ~ "5CF", T ~ "3CF")} else {gathered_gr_AT1G50055$CF <- case_when(gathered_gr_AT1G50055$side =="right" ~ "3CF", T ~ "5CF")}
gathered_AT1G50055 <- data.frame(gathered_gr_AT1G50055)

dt_plot_AT1G50055 <- gathered_AT1G50055 %>%
  filter(cigar=="21M") %>%
  group_by(phase, strand, side, CF, genotype, replicate, target) %>%
  summarise(count=n()) %>%
  filter(genotype %in% genotypes_sub) %>% 
  ungroup()

to_plot_AT1G50055 <- if (dt_plot_AT1G50055[which.max(dt_plot_AT1G50055$count),]$CF == "3CF") {filter(dt_plot_AT1G50055, CF %in% "3CF")} else {filter(dt_plot_AT1G50055, CF %in% "5CF")}

summed_to_plot_AT1G50055 <- to_plot_AT1G50055 %>% 
  filter(genotype %in% genotypes_sub) %>%
  group_by(genotype, replicate) %>% 
  summarise(sum = sum(count)) %>% 
  #summarise(mean_sum = mean(sum)) %>%
  mutate(target = rep("AT1G50055"))

to_plot_AT1G50055$genotype <- factor(to_plot_AT1G50055$genotype, levels = c("wt", "ski2.5", "rrp4", "urt1.1", "GKheso1", "heso1urt1"))
to_plot_AT1G50055$replicate <- factor(to_plot_AT1G50055$replicate, levels = c("R3", "R2", "R1"))

together <- to_plot_AT1G50055 %>% 
  group_by(phase, genotype, target, strand) %>%
  summarise(meancount = mean(count)) %>%
  group_by(genotype, target, strand) %>%
  mutate(ncount = meancount/sum(meancount)) 

tas1b <- together %>% 
  ggplot(aes(x=factor(phase, levels = 1:21), y=ncount)) +
  geom_bar(stat="identity", position= "identity", aes(fill=strand), col="black", size=0.25, alpha=.9) +
  geom_hline(yintercept = seq(0, 4, by = 1), colour = "darkgrey", size = 0.1) +
  geom_hline(yintercept = seq(0, 4, by = 0.25), colour = "darkgrey", size = 0.1) +
  geom_hline(yintercept = seq(0, 4, by = 0.50), colour = "darkgrey", size = 0.1) +
  geom_hline(yintercept = seq(0, 4, by = 0.75), colour = "darkgrey", size = 0.1) +
  geom_vline(xintercept = 1:21, colour = "darkgrey", size = 0.1) + 
  coord_polar() +
  facet_grid(target~genotype, switch="both") + 
  theme(panel.grid = element_blank()) +
  #ggtitle(label = my_title) +
  #theme(panel.grid = element_line(
  #size = 0.25, linetype = 'solid', colour = "grey"),
  #strip.background = element_rect(color="black"), 
  #axis.title.x = element_text(margin = margin(t = 15))) +
  scale_fill_manual(values = c("#304D30", "#B6C4B6")) +
  scale_y_continuous(limits = c(0,1))
ggsave(tas1b, filename = paste0("/isdata/PBgrp/xpj980/figures/tas1bgeneavgURT1201223.svg"), units = "cm", width = 18, height = 18)

summed_to_plot_AT1G50055$genotype <- factor(summed_to_plot_AT1G50055$genotype, levels = c("wt", "ski2.5", "rrp4", "urt1.1", "GKheso1", "heso1urt1"))

read_intensity_plot <- summed_to_plot_AT1G50055 %>% 
  filter(target %in% "AT1G50055") %>%
  ggplot(., aes(x=genotype, y=sum, fill = genotype)) + 
  stat_summary(geom = "bar", fun = "mean", position = "dodge", aes(width=0.5)) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position = "dodge", size = 0.5, aes(width=0.5))+
  #geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~target, scales = "free", nrow = 4) +
  scale_fill_brewer(palette = "YlOrRd") + 
  cowplot::theme_cowplot()
ggsave(read_intensity_plot, filename = paste0("/isdata/PBgrp/xpj980/figures/tas1bgenereadsavgURT1201223.svg"), units = "cm", width = 18, height = 18)
