#### environment ####

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

#### 1Dminus ####
wt_R1 <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1C/wt_R1_seq.txt", header = F, sep = "\t")
wt_R2 <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1C/wt_R2_seq.txt", header = F, sep = "\t")
wt_R3 <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1C/wt_R3_seq.txt", header = F, sep = "\t")
GKheso1_R1 <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1C/GKheso1_R1_seq.txt", header = F, sep = "\t")
GKheso1_R2 <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1C/GKheso1_R2_seq.txt", header = F, sep = "\t")
GKheso1_R3 <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1C/GKheso1_R3_seq.txt", header = F, sep = "\t")
heso1urt1_R1 <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1C/heso1urt1_R1_seq.txt", header = F, sep = "\t")
heso1urt1_R2 <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1C/heso1urt1_R2_seq.txt", header = F, sep = "\t")
heso1urt1_R3 <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1C/heso1urt1_R3_seq.txt", header = F, sep = "\t")
urt1.1_R1 <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1C/urt1_R1_seq.txt", header = F, sep = "\t")
urt1.1_R2 <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1C/urt1_R2_seq.txt", header = F, sep = "\t")
urt1.1_R3 <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1C/urt1_R3_seq.txt", header = F, sep = "\t")

wt_R1_table <- wt_R1 %>%
  mutate(genotype = rep("wt")) %>%
  mutate(rep = rep("R1")) 

wt_R2_table <- wt_R2 %>%
  mutate(genotype = rep("wt")) %>%
  mutate(rep = rep("R2")) 

wt_R3_table <- wt_R3 %>%
  mutate(genotype = rep("wt")) %>%
  mutate(rep = rep("R3")) 

GKHeso1_R1_table <- GKheso1_R1 %>%
  mutate(genotype = rep("GKheso1")) %>%
  mutate(rep = rep("R1"))

GKHeso1_R2_table <- GKheso1_R2 %>%
  mutate(genotype = rep("GKheso1")) %>%
  mutate(rep = rep("R2")) 

GKHeso1_R3_table <- GKheso1_R3 %>%
  mutate(genotype = rep("GKheso1")) %>%
  mutate(rep = rep("R3"))

heso1urt1_R1_table <-heso1urt1_R1 %>%
  mutate(genotype = rep("heso1urt1")) %>%
  mutate(rep = rep("R1"))

heso1urt1_R2_table <-heso1urt1_R2 %>%
  mutate(genotype = rep("heso1urt1")) %>%
  mutate(rep = rep("R2")) 

heso1urt1_R3_table <-heso1urt1_R3 %>%
  mutate(genotype = rep("heso1urt1")) %>%
  mutate(rep = rep("R3")) 

urt1.1_R1_table <-urt1.1_R1 %>%
  mutate(genotype = rep("urt1.1")) %>%
  mutate(rep = rep("R1")) 

urt1.1_R2_table <-urt1.1_R2 %>%
  mutate(genotype = rep("urt1.1")) %>%
  mutate(rep = rep("R2")) 

urt1.1_R3_table <-urt1.1_R3 %>%
  mutate(genotype = rep("urt1.1")) %>%
  mutate(rep = rep("R3")) 

all_dfs <- ls(pattern = "_table$")

# Initialize an empty data frame to store the merged data
merged_df <- data.frame()

# Loop through each data frame and merge it using rbind
for (df_name in all_dfs) {
  df <- get(df_name)  # Get the data frame by name
  merged_df <- rbind(merged_df, df)  # Merge with rbind
}

## add read type 
merged_df <- merged_df %>%
  mutate(read = rep("D1"))

twoU <- merged_df %>% 
  #filter(V1 %in% c("AATGGTCTATTCGCTTGTA", "AATGGTCTATTCGCTTGTAT", "AATGGTCTATTCGCTTGTATT")) %>% 
  mutate(twoU = case_when(V1 == "AATGGTCTATTCGCTTGTA" ~ "0", 
                          V1 == "AATGGTCTATTCGCTTGTAT" ~ "1", 
                          V1 == "AATGGTCTATTCGCTTGTATT" ~ "2", 
                          T ~ "other" )) 

twoU$genotype = factor(twoU$genotype, levels = c("wt", "urt1.1", "GKheso1", "heso1urt1"))

prettyplot <- twoU %>% 
  group_by(genotype, rep, twoU) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(genotype, rep) %>%
  mutate(percentage = count/sum(count) * 100) %>%
  ggplot(., aes(x = twoU, y = percentage, fill = twoU)) + 
  stat_summary(geom = "bar", fun = "mean", position = "dodge") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position = "dodge", size = 0.2) +
  facet_wrap(~genotype, nrow = 1) + 
  scale_fill_brewer(palette = "RdPu") + 
  cowplot::theme_cowplot()
#ggsave(prettyplot, filename = "/isdata/PBgrp/xpj980/figures/triplemutants/TAS1C3D1_twoU.svg", device = "svg", units = "cm", width = 15, height = 10, dpi = 300)

### length dist. 

length_tail <- merged_df %>%
  mutate(length = nchar(V1)) %>%
  group_by(genotype, rep, length) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(genotype, rep) %>%
  mutate(percentage = count/sum(count) * 100)

bars <- length_tail %>%
  ggplot(aes(x = length, y = percentage)) + 
  geom_bar(stat = "identity") + 
  facet_grid(genotype~rep)
#ggsave(bars, filename = "/isdata/PBgrp/xpj980/figures/triplemutants/length_TAS1C3D1.svg", device = "svg")

length_tail <- merged_df %>%
  mutate(length = nchar(V1)) %>%
  filter(length < 25) %>%
  group_by(genotype, rep, length) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(genotype, rep) %>%
  mutate(percentage = count/sum(count) * 100)

length_tail$genotype = factor(length_tail$genotype, levels = c("wt", "urt1.1", "GKheso1", "heso1urt1"))

prettyplot <- length_tail %>%
  ggplot(., aes(x = length, y = percentage, fill = genotype)) + 
  stat_summary(geom = "bar", fun = "mean", position = "dodge") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position = "dodge", size = 0.2) +
  #facet_wrap(~genotype, nrow = 4) + 
  scale_fill_brewer(palette = "RdPu") + 
  cowplot::theme_cowplot()
#ggsave(prettyplot, filename = "/isdata/PBgrp/xpj980/figures/triplemutants/length_TAS1C3D1_prettyplot.svg", device = "svg", units = "cm", width = 15, height = 10, dpi = 300)

### incorporate RP10M ###

featureCounts_sum <- read.table("/binf-isilon/PBgrp/xpj980/datascratch/triplemutants/04.STAR_mappings/old_style_mapping/featureCountsCarlottagtfmodifiedcountall.summary", header = T, row.names = 1)  %>%
  as.data.frame()
new_names <- names(featureCounts_sum) %>% str_remove(".fastq.gz.trimmed.fq.gzAligned.out.sam")

featureCounts_renamed <- featureCounts_sum
names(featureCounts_renamed) <- new_names

total_alignments <- featureCounts_renamed[1,] %>%
  gather(key = "genotype", value = "totalread") %>% 
  separate(genotype, c("genotype", "replicate"), sep="_") %>%
  filter(genotype %in% c("wt", "urt1.1", "GKheso1", "heso1urt1")) %>%
  rename(replicate = "rep")

left_joined <- merge(merged_df, total_alignments, by = c("genotype", "rep"))

length_tail <- left_joined %>%
  mutate(length = nchar(V1)) %>%
  filter(length < 25) %>%
  group_by(genotype, rep, length, totalread) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(genotype, rep) %>%
  mutate(RP10M = ifelse(count == 0, 0, 10000000 * count / totalread)) %>%
  ungroup() 

length_tail$genotype = factor(length_tail$genotype, levels = c("wt", "urt1.1", "GKheso1", "heso1urt1"))

prettyplot <- length_tail %>%
  ggplot(., aes(x = length, y = RP10M, fill = genotype)) + 
  stat_summary(geom = "bar", fun = "mean", position = "dodge") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position = "dodge", size = 0.2) +
  #facet_wrap(~genotype, nrow = 4) + 
  scale_fill_brewer(palette = "RdPu") + 
  cowplot::theme_cowplot()
#ggsave(prettyplot, filename = "/isdata/PBgrp/xpj980/figures/triplemutants/length_TAS1C3D1_prettyplot_RP10M.svg", device = "svg", units = "cm", width = 15, height = 10, dpi = 300)

### size compared to 2D minus and plus ###
#### 2D plus ####

wt_R1_plus <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1C/2D_plus/wt_R1_TAS1C_2Dplus_seq.txt", header = F, sep = "\t")
wt_R2_plus <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1C/2D_plus/wt_R2_TAS1C_2Dplus_seq.txt", header = F, sep = "\t")
wt_R3_plus <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1C/2D_plus/wt_R3_TAS1C_2Dplus_seq.txt", header = F, sep = "\t")
GKheso_R1_plus <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1C/2D_plus/GKheso1_R1_TAS1C_2Dplus_seq.txt", header = F, sep = "\t")
GKheso_R2_plus <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1C/2D_plus/GKheso1_R2_TAS1C_2Dplus_seq.txt", header = F, sep = "\t")
GKheso_R3_plus <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1C/2D_plus/GKheso1_R3_TAS1C_2Dplus_seq.txt", header = F, sep = "\t")
urt1_R1_plus <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1C/2D_plus/urt1_R1_TAS1C_2Dplus_seq.txt", header = F, sep = "\t")
urt1_R2_plus <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1C/2D_plus/urt1_R2_TAS1C_2Dplus_seq.txt", header = F, sep = "\t")
urt1_R3_plus <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1C/2D_plus/urt1_R3_TAS1C_2Dplus_seq.txt", header = F, sep = "\t")
heso1urt1_R1_plus <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1C/2D_plus/heso1urt1_R1_TAS1C_2Dplus_seq.txt", header = F, sep = "\t")
heso1urt1_R2_plus <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1C/2D_plus/heso1urt1_R2_TAS1C_2Dplus_seq.txt", header = F, sep = "\t")
heso1urt1_R3_plus <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1C/2D_plus/heso1urt1_R3_TAS1C_2Dplus_seq.txt", header = F, sep = "\t")

wt_R1_table_plus <- wt_R1_plus %>%
  mutate(genotype = rep("wt")) %>%
  mutate(rep = rep("R1"))

wt_R2_table_plus <- wt_R2_plus %>%
  mutate(genotype = rep("wt")) %>%
  mutate(rep = rep("R2"))

wt_R3_table_plus <- wt_R3_plus %>%
  mutate(genotype = rep("wt")) %>%
  mutate(rep = rep("R3")) 

GKHeso1_R1_table_plus <- GKheso_R1_plus %>%
  mutate(genotype = rep("GKheso1")) %>%
  mutate(rep = rep("R1"))

GKHeso1_R2_table_plus <- GKheso_R2_plus %>%
  mutate(genotype = rep("GKheso1")) %>%
  mutate(rep = rep("R2")) 

GKHeso1_R3_table_plus <- GKheso_R3_plus %>%
  mutate(genotype = rep("GKheso1")) %>%
  mutate(rep = rep("R3"))

heso1urt1_R1_table_plus <-heso1urt1_R1_plus %>%
  mutate(genotype = rep("heso1urt1")) %>%
  mutate(rep = rep("R1"))

heso1urt1_R2_table_plus <-heso1urt1_R2_plus %>%
  mutate(genotype = rep("heso1urt1")) %>%
  mutate(rep = rep("R2")) 

heso1urt1_R3_table_plus <-heso1urt1_R3_plus %>%
  mutate(genotype = rep("heso1urt1")) %>%
  mutate(rep = rep("R3")) 

urt1.1_R1_table_plus <-urt1_R1_plus %>%
  mutate(genotype = rep("urt1.1")) %>%
  mutate(rep = rep("R1")) 

urt1.1_R2_table_plus <-urt1_R2_plus %>%
  mutate(genotype = rep("urt1.1")) %>%
  mutate(rep = rep("R2")) 

urt1.1_R3_table_plus <-urt1_R3_plus %>%
  mutate(genotype = rep("urt1.1")) %>%
  mutate(rep = rep("R3")) 

all_dfs <- ls(pattern = "_table_plus$")

# Initialize an empty data frame to store the merged data
merged_df_plus <- data.frame()

# Loop through each data frame and merge it using rbind
for (df_name in all_dfs) {
  df <- get(df_name)  # Get the data frame by name
  merged_df_plus <- rbind(merged_df_plus, df)  # Merge with rbind
}

## add read type 
merged_df_plus <- merged_df_plus %>%
  mutate(read = rep("D2plus"))

length_tail <- merged_df_plus %>%
  mutate(length = nchar(V1)) %>%
  filter(length < 25) %>%
  group_by(genotype, rep, length) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(genotype, rep) %>%
  mutate(percentage = count/sum(count) * 100)

length_tail$genotype = factor(length_tail$genotype, levels = c("wt", "urt1.1", "GKheso1", "heso1urt1"))

prettyplot <- length_tail %>%
  ggplot(., aes(x = length, y = percentage, fill = genotype)) + 
  stat_summary(geom = "bar", fun = "mean", position = "dodge") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position = "dodge", size = 0.2) +
  #facet_wrap(~genotype, nrow = 4) + 
  scale_fill_brewer(palette = "RdPu") + 
  cowplot::theme_cowplot()

together <- rbind(merged_df, merged_df_plus)

length_tail <- together %>%
  mutate(length = nchar(V1)) %>%
  filter(length < 25) %>%
  group_by(genotype, rep, length, read) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(genotype, rep, read) %>%
  mutate(percentage = count/sum(count) * 100)

length_tail$genotype = factor(length_tail$genotype, levels = c("wt", "urt1.1", "GKheso1", "heso1urt1"))

prettyplot <- length_tail %>%
  ggplot(., aes(x = length, y = percentage, fill = genotype)) + 
  stat_summary(geom = "bar", fun = "mean", position = "dodge") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position = "dodge", size = 0.2) +
  facet_wrap(~read) +
  scale_x_continuous(breaks = c(19:24)) +
  scale_fill_brewer(palette = "RdPu") + 
  cowplot::theme_cowplot()

# RP10M 

left_joined <- merge(together, total_alignments, by = c("genotype", "rep"))

length_tail <- left_joined %>%
  mutate(length = nchar(V1)) %>%
  filter(length < 25) %>%
  group_by(genotype, rep, length, read, totalread) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(genotype, rep, read) %>%
  mutate(RP10M = ifelse(count == 0, 0, 10000000 * count / totalread)) %>%
  ungroup() 

length_tail$genotype = factor(length_tail$genotype, levels = c("wt", "urt1.1", "GKheso1", "heso1urt1"))

prettyplot <- length_tail %>%
  ggplot(., aes(x = length, y = RP10M, fill = genotype)) + 
  stat_summary(geom = "bar", fun = "mean", position = "dodge") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position = "dodge", size = 0.2) +
  facet_wrap(~read, scales = "free") + 
  scale_x_continuous(breaks = c(19:24)) +
  scale_fill_brewer(palette = "RdPu") + 
  cowplot::theme_cowplot()

#### 2D minus #### 

wt_R1_minus <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1C/2D_minus/wt_R1_TAS1C_2Dminus_seq.txt", header = F, sep = "\t")
wt_R2_minus <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1C/2D_minus/wt_R2_TAS1C_2Dminus_seq.txt", header = F, sep = "\t")
wt_R3_minus <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1C/2D_minus/wt_R3_TAS1C_2Dminus_seq.txt", header = F, sep = "\t")
GKheso_R1_minus <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1C/2D_minus/GKheso1_R1_TAS1C_2Dminus_seq.txt", header = F, sep = "\t")
GKheso_R2_minus <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1C/2D_minus/GKheso1_R2_TAS1C_2Dminus_seq.txt", header = F, sep = "\t")
GKheso_R3_minus <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1C/2D_minus/GKheso1_R3_TAS1C_2Dminus_seq.txt", header = F, sep = "\t")
urt1_R1_minus <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1C/2D_minus/urt1_R1_TAS1C_2Dminus_seq.txt", header = F, sep = "\t")
urt1_R2_minus <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1C/2D_minus/urt1_R2_TAS1C_2Dminus_seq.txt", header = F, sep = "\t")
urt1_R3_minus <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1C/2D_minus/urt1_R3_TAS1C_2Dminus_seq.txt", header = F, sep = "\t")
heso1urt1_R1_minus <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1C/2D_minus/heso1urt1_R1_TAS1C_2Dminus_seq.txt", header = F, sep = "\t")
heso1urt1_R2_minus <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1C/2D_minus/heso1urt1_R2_TAS1C_2Dminus_seq.txt", header = F, sep = "\t")
heso1urt1_R3_minus <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1C/2D_minus/heso1urt1_R3_TAS1C_2Dminus_seq.txt", header = F, sep = "\t")

wt_R1_table_minus <- wt_R1_minus %>%
  mutate(genotype = rep("wt")) %>%
  mutate(rep = rep("R1"))

wt_R2_table_minus <- wt_R2_minus %>%
  mutate(genotype = rep("wt")) %>%
  mutate(rep = rep("R2"))

wt_R3_table_minus <- wt_R3_minus %>%
  mutate(genotype = rep("wt")) %>%
  mutate(rep = rep("R3")) 

GKHeso1_R1_table_minus <- GKheso_R1_minus %>%
  mutate(genotype = rep("GKheso1")) %>%
  mutate(rep = rep("R1"))

GKHeso1_R2_table_minus <- GKheso_R2_minus %>%
  mutate(genotype = rep("GKheso1")) %>%
  mutate(rep = rep("R2")) 

GKHeso1_R3_table_minus <- GKheso_R3_minus %>%
  mutate(genotype = rep("GKheso1")) %>%
  mutate(rep = rep("R3"))

heso1urt1_R1_table_minus <-heso1urt1_R1_minus %>%
  mutate(genotype = rep("heso1urt1")) %>%
  mutate(rep = rep("R1"))

heso1urt1_R2_table_minus <-heso1urt1_R2_minus %>%
  mutate(genotype = rep("heso1urt1")) %>%
  mutate(rep = rep("R2")) 

heso1urt1_R3_table_minus <-heso1urt1_R3_minus %>%
  mutate(genotype = rep("heso1urt1")) %>%
  mutate(rep = rep("R3")) 

urt1.1_R1_table_minus <-urt1_R1_minus %>%
  mutate(genotype = rep("urt1.1")) %>%
  mutate(rep = rep("R1")) 

urt1.1_R2_table_minus <-urt1_R2_minus %>%
  mutate(genotype = rep("urt1.1")) %>%
  mutate(rep = rep("R2")) 

urt1.1_R3_table_minus <-urt1_R3_minus %>%
  mutate(genotype = rep("urt1.1")) %>%
  mutate(rep = rep("R3")) 

all_dfs <- ls(pattern = "_table_minus$")

# Initialize an empty data frame to store the merged data
merged_df_minus <- data.frame()

# Loop through each data frame and merge it using rbind
for (df_name in all_dfs) {
  df <- get(df_name)  # Get the data frame by name
  merged_df_minus <- rbind(merged_df_minus, df)  # Merge with rbind
}

## add read type 
merged_df_minus <- merged_df_minus %>%
  mutate(read = rep("D2minus"))

length_tail <- merged_df_minus %>%
  mutate(length = nchar(V1)) %>%
  filter(length < 25) %>%
  group_by(genotype, rep, length) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(genotype, rep) %>%
  mutate(percentage = count/sum(count) * 100)

length_tail$genotype = factor(length_tail$genotype, levels = c("wt", "urt1.1", "GKheso1", "heso1urt1"))

prettyplot <- length_tail %>%
  ggplot(., aes(x = length, y = percentage, fill = genotype)) + 
  stat_summary(geom = "bar", fun = "mean", position = "dodge") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position = "dodge", size = 0.2) +
  #facet_wrap(~genotype, nrow = 4) + 
  scale_fill_brewer(palette = "RdPu") + 
  cowplot::theme_cowplot()

together <- rbind(merged_df, merged_df_plus, merged_df_minus)
#write.xlsx(x = together, file = "/binf-isilon/PBgrp/xpj980/R_scripts_collected/tas1C_tails.xlsx")

length_tail <- together %>%
  mutate(length = nchar(V1)) %>%
  filter(length < 25) %>%
  group_by(genotype, rep, length, read) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(genotype, rep, read) %>%
  mutate(percentage = count/sum(count) * 100)

length_tail$genotype = factor(length_tail$genotype, levels = c("wt", "urt1.1", "GKheso1", "heso1urt1"))
length_tail$read = factor(length_tail$read, levels = c("D1", "D2plus", "D2minus"))

prettyplot <- length_tail %>%
  ggplot(., aes(x = length, y = percentage, fill = genotype)) + 
  stat_summary(geom = "bar", fun = "mean", position = "dodge") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position = "dodge", size = 0.2) +
  facet_wrap(~read) +
  scale_x_continuous(breaks = c(19:24)) +
  scale_fill_brewer(palette = "Set2") + 
  #scale_fill_manual(values=cbPalette) + 
  cowplot::theme_cowplot()
#ggsave(prettyplot, filename = "/isdata/PBgrp/xpj980/figures/triplemutants/length_TAS1C3D1D2_all.svg", device = "svg", units = "cm", width = 15, height = 10, dpi = 300)

# RP10M 

left_joined <- merge(together, total_alignments, by = c("genotype", "rep"))

length_tail <- left_joined %>%
  mutate(length = nchar(V1)) %>%
  filter(length < 25) %>%
  group_by(genotype, rep, length, read, totalread) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(genotype, rep, read) %>%
  mutate(RP10M = ifelse(count == 0, 0, 10000000 * count / totalread)) %>%
  ungroup() 

length_tail$genotype = factor(length_tail$genotype, levels = c("wt", "urt1.1", "GKheso1", "heso1urt1"))
length_tail$read = factor(length_tail$read, levels = c("D1", "D2plus", "D2minus"))

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

prettyplot <- length_tail %>%
  ggplot(., aes(x = length, y = RP10M, fill = genotype)) + 
  stat_summary(geom = "bar", fun = "mean", position = "dodge") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position = "dodge", size = 0.2) +
  facet_wrap(~read, scales = "free") + 
  scale_x_continuous(breaks = c(19:24)) +
  #scale_fill_brewer(palette = "Set2") + 
  scale_fill_manual(values=cbPalette) + 
  cowplot::theme_cowplot()
#ggsave(prettyplot, filename = "/isdata/PBgrp/xpj980/figures/triplemutants/length_TAS1C3D1D2_all_RP10M.svg", device = "svg", units = "cm", width = 15, height = 10, dpi = 300)

### simple versions only with WT and heso1 ###

length_tail <- together %>%
  filter(genotype %in% c("wt", "GKheso1")) %>%
  mutate(length = nchar(V1)) %>%
  filter(length < 25) %>%
  group_by(genotype, rep, length, read) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(genotype, rep, read) %>%
  mutate(percentage = count/sum(count) * 100)
#write.xlsx(x = length_tail, file = "/binf-isilon/PBgrp/xpj980/R_scripts_collected/Figure7D.xlsx")

length_tail$genotype = factor(length_tail$genotype, levels = c("wt", "GKheso1"))
length_tail$read = factor(length_tail$read, levels = c("D1", "D2plus", "D2minus"))
#### Figure 7D #####
prettyplot <- length_tail %>%
  ggplot(., aes(x = length, y = percentage, fill = genotype)) + 
  stat_summary(geom = "bar", fun = "mean", position = "dodge") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position = "dodge", size = 0.2) +
  facet_wrap(~read) +
  #facet_grid(read~length)
  scale_x_continuous(breaks = c(19:24)) +
  #scale_fill_brewer(palette = "Set2") + 
  scale_fill_manual(values=cbPalette) + 
  cowplot::theme_cowplot()
#ggsave(prettyplot, filename = "/isdata/PBgrp/xpj980/figures/triplemutants/length_TAS1C3D1D2_2geno.svg", device = "svg", units = "cm", width = 15, height = 10, dpi = 300)

# and the numbers to put in the figure: 

length_tail <- together %>%
  filter(genotype %in% c("wt", "GKheso1")) %>%
  mutate(length = nchar(V1)) %>%
  filter(length < 25) %>%
  group_by(genotype, rep, length, read) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(genotype, rep, read) %>%
  summarise(sum = sum(count))

# RP10M 

left_joined <- merge(together, total_alignments, by = c("genotype", "rep"))

length_tail <- left_joined %>%
  filter(genotype %in% c("wt", "GKheso1")) %>%
  mutate(length = nchar(V1)) %>%
  filter(length < 25) %>%
  group_by(genotype, rep, length, read, totalread) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(genotype, rep, read) %>%
  mutate(RP10M = ifelse(count == 0, 0, 10000000 * count / totalread)) %>%
  ungroup() 

length_tail$genotype = factor(length_tail$genotype, levels = c("wt", "GKheso1"))
length_tail$read = factor(length_tail$read, levels = c("D1", "D2plus", "D2minus"))

# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

prettyplot <- length_tail %>%
  ggplot(., aes(x = length, y = RP10M, fill = genotype)) + 
  stat_summary(geom = "bar", fun = "mean", position = "dodge") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position = "dodge", size = 0.2) +
  facet_wrap(~read, scales = "free") + 
  scale_x_continuous(breaks = c(19:24)) +
  #scale_fill_brewer(palette = "Set2") + 
  scale_fill_manual(values=cbPalette) + 
  cowplot::theme_cowplot()
#ggsave(prettyplot, filename = "/isdata/PBgrp/xpj980/figures/triplemutants/length_TAS1C3D1D2_2geno_RP10M.svg", device = "svg", units = "cm", width = 15, height = 10, dpi = 300)

#### composition ####

get_tail <- function(x) {
  tail <- str_remove(x, "AATGGTCTATTCGCTTGTA")
  return(list(tail = tail))
}

tail <- get_tail(merged_df$V1)

DF2 <- merged_df %>% 
  mutate( 
    tail = get_tail(V1)$tail, 
    tail_len = nchar(tail))

### insert UU-tail bar ###

tail_cat_df <- DF2 %>% 
  mutate(tail_cat = case_when(tail == "TT" ~ "UU-tail", 
                              tail_len == 0 ~ "no_tail", 
                              tail_len == 1 ~ "1-nt", 
                              tail_len == 2 ~ "2-nt", 
                              tail_len == 3 ~ "3-nt", 
                              tail_len > 3 ~ "longer")) %>%
  group_by(genotype, rep, tail_cat) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(genotype, rep) %>%
  mutate(percentage = count/sum(count) * 100)

tail_cat_df$genotype = factor(tail_cat_df$genotype, levels = c("wt", "urt1.1", "GKheso1", "heso1urt1"))
tail_cat_df$tail_cat = factor(tail_cat_df$tail_cat, levels = c("no_tail", "1-nt", "2-nt", "UU-tail", "3-nt", "longer"))

test <- ggplot(tail_cat_df, aes(x = tail_cat, y = percentage, fill = tail_cat)) + 
  stat_summary(geom = "bar", fun = "mean", position = "dodge") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position = "dodge", size = 0.2) +
  facet_wrap(~genotype, nrow = 1) + 
  scale_fill_brewer(palette = "Blues") + 
  cowplot::theme_cowplot()
#ggsave(test, filename = "/isdata/PBgrp/xpj980/figures/triplemutants/TAS1C3D1_percentage.svg", device = "svg", units = "cm", width = 15, height = 10, dpi = 300)

### RP10M ###

featureCounts_sum <- read.table("/binf-isilon/PBgrp/xpj980/datascratch/triplemutants/04.STAR_mappings/old_style_mapping/featureCountsCarlottagtfmodifiedcountall.summary", header = T, row.names = 1)  %>%
  as.data.frame()
new_names <- names(featureCounts_sum) %>% str_remove(".fastq.gz.trimmed.fq.gzAligned.out.sam")

featureCounts_renamed <- featureCounts_sum
names(featureCounts_renamed) <- new_names

total_alignments <- featureCounts_renamed[1,] %>%
  gather(key = "genotype", value = "totalread") %>% 
  separate(genotype, c("genotype", "replicate"), sep="_") %>%
  filter(genotype %in% c("wt", "urt1.1", "GKheso1", "heso1urt1")) %>%
  rename(replicate = "rep")

left_joined <- merge(tail_cat_df, total_alignments, by = c("genotype", "rep")) %>%
  mutate(RP10M = 10000000*.$count/.$totalread)

left_joined$genotype = factor(left_joined$genotype, levels = c("wt", "urt1.1", "GKheso1", "heso1urt1"))
left_joined$tail_cat = factor(left_joined$tail_cat, levels = c("no_tail", "1-nt", "2-nt", "UU-tail", "3-nt", "longer"))

test <- ggplot(left_joined, aes(x = tail_cat, y = RP10M, fill = tail_cat)) + 
  stat_summary(geom = "bar", fun = "mean", position = "dodge") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position = "dodge", size = 0.2) +
  facet_wrap(~genotype, nrow = 1) + 
  scale_fill_brewer(palette = "Blues") + 
  cowplot::theme_cowplot()

#ggsave(test, filename = "/isdata/PBgrp/xpj980/figures/triplemutants/TAS1C3D1_RP10M.svg", device = "svg", units = "cm", width = 15, height = 10, dpi = 300)

### second way where UU is counted twice to show % of di-nucleotides AND % of them that are UU  

tail_cat_df <- DF2 %>% 
  mutate(tail_cat = case_when( 
    tail_len == 0 ~ "no_tail", 
    tail_len == 1 ~ "1-nt", 
    tail_len == 2 ~ "2-nt", 
    tail_len == 3 ~ "3-nt", 
    tail_len > 3 ~ "longer")) %>%
  group_by(genotype, rep, tail_cat) %>%
  mutate(total = n()) %>%
  mutate(tail_cat2 = case_when(tail == "TT" ~ "UU-tail", 
                               T ~tail_cat)) %>%
  group_by(genotype, rep, tail_cat2) %>%
  mutate(total2 = n()) %>%
  mutate(total3 = case_when(tail_cat2 == "UU-tail" ~ total2, 
                            T ~ total)) %>%
  group_by(genotype, rep, tail_cat2, total3) %>%
  summarise() %>%
  ungroup() %>%
  group_by(genotype, rep) %>%
  mutate(percentage = total3/sum(total3) * 100) %>%
  ungroup()
#write.xlsx(x = tail_cat_df, file = "/binf-isilon/PBgrp/xpj980/R_scripts_collected/Figure7C.xlsx")

tail_cat_df$genotype = factor(tail_cat_df$genotype, levels = c("wt", "urt1.1", "GKheso1", "heso1urt1"))
tail_cat_df$tail_cat2 = factor(tail_cat_df$tail_cat2, levels = c("no_tail", "1-nt", "2-nt", "UU-tail", "3-nt", "longer"))

left_joined <- merge(tail_cat_df, total_alignments, by = c("genotype", "rep")) %>%
  mutate(RP10M = 10000000*.$total3/.$totalread)

left_joined$genotype = factor(left_joined$genotype, levels = c("wt", "urt1.1", "GKheso1", "heso1urt1"))
left_joined$tail_cat2 = factor(left_joined$tail_cat2, levels = c("no_tail", "1-nt", "2-nt", "UU-tail", "3-nt", "longer"))

#### Figure 7C #####
test <- ggplot(left_joined, aes(x = tail_cat2, y = RP10M, fill = tail_cat2)) + 
  stat_summary(geom = "bar", fun = "mean", position = "dodge") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position = "dodge", size = 0.2) +
  facet_wrap(~genotype, nrow = 1) + 
  scale_fill_brewer(palette = "Blues") + 
  cowplot::theme_cowplot()
#ggsave(test, filename = "/isdata/PBgrp/xpj980/figures/triplemutants/TAS1C3D1_RP10M_New2.svg", device = "svg", units = "cm", width = 15, height = 10, dpi = 300)

## and in percentage 

test <- ggplot(tail_cat_df, aes(x = tail_cat2, y = percentage, fill = tail_cat2)) + 
  stat_summary(geom = "bar", fun = "mean", position = "dodge") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position = "dodge", size = 0.2) +
  facet_wrap(~genotype, nrow = 1) + 
  scale_fill_brewer(palette = "Blues") + 
  cowplot::theme_cowplot()
#ggsave(test, filename = "/isdata/PBgrp/xpj980/figures/triplemutants/TAS1C3D1_percentage_newbutbad.svg", device = "svg", units = "cm", width = 15, height = 10, dpi = 300)

### 1nt and 2nt composition ### 

one_nt <- DF2 %>% 
  filter(tail_len == 1) %>%
  group_by(genotype, rep, tail) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(genotype, rep) %>%
  mutate(percentage = count/sum(count) * 100)

one_nt$genotype = factor(one_nt$genotype, levels = c("wt", "urt1.1", "GKheso1", "heso1urt1"))

test <- ggplot(one_nt, aes(x = tail, y = percentage, fill = tail)) + 
  stat_summary(geom = "bar", fun = "mean", position = "dodge") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position = "dodge", size = 0.2) +
  facet_wrap(~genotype, nrow = 1) + 
  scale_fill_manual(values=rev(cbPalette)) + 
  #scale_fill_brewer(palette = "RdPu") + 
  cowplot::theme_cowplot()
#ggsave(test, filename = "/isdata/PBgrp/xpj980/figures/triplemutants/TAS1C3D1_onent.svg", device = "svg", units = "cm", width = 15, height = 10, dpi = 300)

two_nt <- DF2 %>% 
  filter(tail_len == 2) %>%
  group_by(genotype, rep, tail) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(genotype, rep) %>%
  mutate(percentage = count/sum(count) * 100)

two_nt$genotype = factor(two_nt$genotype, levels = c("wt", "urt1.1", "GKheso1", "heso1urt1"))

test <- ggplot(two_nt, aes(x = tail, y = percentage, fill = tail)) + 
  stat_summary(geom = "bar", fun = "mean", position = "dodge") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position = "dodge", size = 0.2) +
  facet_wrap(~genotype, nrow = 1) + 
  #scale_fill_brewer(palette = "RdPu") + 
  cowplot::theme_cowplot()
#ggsave(test, filename = "/isdata/PBgrp/xpj980/figures/triplemutants/TAS1C3D1_twont.svg", device = "svg", units = "cm", width = 15, height = 10, dpi = 300)

### in RPM and with different y-axis 

two_nt <- DF2 %>% 
  filter(tail_len == 2) %>%
  group_by(genotype, rep, tail) %>%
  summarise(count = n())

left_joined <- merge(two_nt, total_alignments, by = c("genotype", "rep")) %>%
  mutate(RP10M = 10000000*.$count/.$totalread)

left_joined$genotype = factor(left_joined$genotype, levels = c("wt", "urt1.1", "GKheso1", "heso1urt1"))

test <- ggplot(left_joined, aes(x = tail, y = RP10M, fill = rep)) + 
  stat_summary(geom = "bar", fun = "mean", position = "dodge") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position = "dodge", size = 0.2) +
  facet_wrap(~genotype, nrow = 1, scales = "free") + 
  #scale_fill_brewer(palette = "RdPu") + 
  cowplot::theme_cowplot()
#ggsave(test, filename = "/isdata/PBgrp/xpj980/figures/triplemutants/TAS1C3D1_twontRP10M_free.svg", device = "svg", units = "cm", width = 15, height = 10, dpi = 300)

### two nt sum per rep

two_nt <- DF2 %>% 
  filter(tail_len == 2) %>%
  group_by(genotype, tail) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(genotype) %>%
  mutate(percentage = count/sum(count) * 100)
#write.xlsx(x = two_nt, file = "/binf-isilon/PBgrp/xpj980/R_scripts_collected/EV6B.xlsx")

two_nt$genotype = factor(two_nt$genotype, levels = c("wt", "urt1.1", "GKheso1", "heso1urt1"))

#### EV6B ####
test <- ggplot(two_nt, aes(x = tail, y = percentage, fill = tail)) + 
  geom_bar(stat = "identity") +
  #stat_summary(geom = "bar", fun = "mean", position = "dodge") + 
  #stat_summary(fun.data = mean_se, geom = "errorbar", position = "dodge", size = 0.2) +
  facet_wrap(~genotype, nrow = 1) + 
  #scale_fill_brewer(palette = "RdPu") + 
  cowplot::theme_cowplot()
#ggsave(test, filename = "/isdata/PBgrp/xpj980/figures/triplemutants/TAS1C3D1_twont_repSum.svg", device = "svg", units = "cm", width = 15, height = 10, dpi = 300)

### composition 

test <- DF2 %>% 
  mutate(frac.t = str_count(DF2$tail, "T")/tail_len, 
         frac.c = str_count(DF2$tail, "C")/tail_len,
         frac.g = str_count(DF2$tail, "G")/tail_len,
         frac.a = str_count(DF2$tail, "A")/tail_len) %>% 
  filter(., tail == "" | frac.t > 0.5 | frac.c > 0.5 | frac.g > 0.5 | frac.a > 0.5)

test2 <- test %>% 
  group_by(genotype, tail) %>%
  filter(tail_len < 20) %>% 
  summarise(count = n()) %>% 
  group_by(genotype) %>% 
  mutate(ncount = count/sum(count)) %>% 
  filter(tail != "") %>% 
  mutate(frac.u = str_count(tail, "T")/nchar(tail), 
         frac.c = str_count(tail, "C")/nchar(tail),
         frac.g = str_count(tail, "G")/nchar(tail),
         frac.a = str_count(tail, "A")/nchar(tail)) %>% 
  gather(base, basecount, contains("frac"))

test2$genotype = factor(test2$genotype, levels = c("wt", "urt1.1", "GKheso1", "heso1urt1"))    

plot <- ggplot(test2, aes(x = genotype, y = (basecount*ncount)*100, fill = base %>% str_remove("frac.") %>% toupper())) +
  geom_bar(stat="identity", width = 0.7) +
  labs(x = "genotype", y = "% of reads with an ear", fill = "base") +
  cowplot::theme_cowplot() +
  cowplot::background_grid() +
  #facet_wrap(~ rep) + 
  #geom_text(aes(label = paste0(Round_off,"%")), 
  #position = position_stack(vjust = 0.5), size = 4) +
  scale_fill_manual(values = wes_palette("GrandBudapest1", n = 4)) +
  ggtitle("Proportion of reads with an ear and composition of the ear") 
#ggsave(plot, filename = "/isdata/PBgrp/xpj980/figures/triplemutants/ear_TAS1C3D1.svg", device = "svg", units = "cm", width = 15, height = 15, dpi = 300)


### make table with all TAS transcripts ####
#### TAS1C ####
wt_R1_TAS1C <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1C/wt_R1_seq.txt", header = F, sep = "\t") %>%
  mutate(genotype = rep("wt")) %>%
  mutate(rep = rep("R1")) %>%
  mutate(TAS = rep("TAS1C"))

wt_R2_TAS1C <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1C/wt_R2_seq.txt", header = F, sep = "\t") %>%
  mutate(genotype = rep("wt")) %>%
  mutate(rep = rep("R2")) %>%
  mutate(TAS = rep("TAS1C"))

wt_R3_TAS1C <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1C/wt_R3_seq.txt", header = F, sep = "\t") %>%
  mutate(genotype = rep("wt")) %>%
  mutate(rep = rep("R3")) %>%
  mutate(TAS = rep("TAS1C"))

GKheso1_R1_TAS1C <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1C/GKheso1_R1_seq.txt", header = F, sep = "\t") %>%
  mutate(genotype = rep("heso1")) %>%
  mutate(rep = rep("R1")) %>%
  mutate(TAS = rep("TAS1C"))

GKheso1_R2_TAS1C <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1C/GKheso1_R2_seq.txt", header = F, sep = "\t") %>%
  mutate(genotype = rep("heso1")) %>%
  mutate(rep = rep("R2")) %>%
  mutate(TAS = rep("TAS1C"))

GKheso1_R3_TAS1C <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1C/GKheso1_R3_seq.txt", header = F, sep = "\t") %>%
  mutate(genotype = rep("heso1")) %>%
  mutate(rep = rep("R3")) %>%
  mutate(TAS = rep("TAS1C"))

heso1urt1_R1_TAS1C <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1C/heso1urt1_R1_seq.txt", header = F, sep = "\t") %>%
  mutate(genotype = rep("heso1urt1")) %>%
  mutate(rep = rep("R1")) %>%
  mutate(TAS = rep("TAS1C"))

heso1urt1_R2_TAS1C <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1C/heso1urt1_R2_seq.txt", header = F, sep = "\t") %>%
  mutate(genotype = rep("heso1urt1")) %>%
  mutate(rep = rep("R2")) %>%
  mutate(TAS = rep("TAS1C"))

heso1urt1_R3_TAS1C <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1C/heso1urt1_R3_seq.txt", header = F, sep = "\t") %>%
  mutate(genotype = rep("heso1urt1")) %>%
  mutate(rep = rep("R3")) %>%
  mutate(TAS = rep("TAS1C"))

urt1.1_R1_TAS1C <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1C/urt1_R1_seq.txt", header = F, sep = "\t") %>%
  mutate(genotype = rep("urt1")) %>%
  mutate(rep = rep("R1")) %>%
  mutate(TAS = rep("TAS1C"))

urt1.1_R2_TAS1C <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1C/urt1_R2_seq.txt", header = F, sep = "\t") %>%
  mutate(genotype = rep("urt1")) %>%
  mutate(rep = rep("R2")) %>%
  mutate(TAS = rep("TAS1C"))

urt1.1_R3_TAS1C <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1C/urt1_R3_seq.txt", header = F, sep = "\t") %>%
  mutate(genotype = rep("urt1")) %>%
  mutate(rep = rep("R3")) %>%
  mutate(TAS = rep("TAS1C"))

TAS1C_dfs <- ls(pattern = "_TAS1C")

# Initialize an empty data frame to store the merged data
TAS1C_merged_df <- data.frame()

# Loop through each data frame and merge it using rbind
for (df_name in TAS1C_dfs) {
  df <- get(df_name)  # Get the data frame by name
  TAS1C_merged_df <- rbind(TAS1C_merged_df, df)  # Merge with rbind
}
write.xlsx(TAS1C_merged_df, file = "/isdata/PBgrp/xpj980/datascratch/triplemutants/TAS1C.xlsx")

### TAS1B ####

wt_R1_TAS1B <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1B/wt_R1_TASB_seq.txt", header = F, sep = "\t") %>%
  mutate(genotype = rep("wt")) %>%
  mutate(rep = rep("R1")) %>%
  mutate(TAS = rep("TAS1B"))

wt_R2_TAS1B <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1B/wt_R2_TASB_seq.txt", header = F, sep = "\t") %>%
  mutate(genotype = rep("wt")) %>%
  mutate(rep = rep("R2")) %>%
  mutate(TAS = rep("TAS1B"))

wt_R3_TAS1B <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1B/wt_R3_TASB_seq.txt", header = F, sep = "\t") %>%
  mutate(genotype = rep("wt")) %>%
  mutate(rep = rep("R3")) %>%
  mutate(TAS = rep("TAS1B"))

GKheso1_R1_TAS1B <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1B/GKHeso_R1_TASB_seq.txt", header = F, sep = "\t") %>%
  mutate(genotype = rep("heso1")) %>%
  mutate(rep = rep("R1")) %>%
  mutate(TAS = rep("TAS1B"))

GKheso1_R2_TAS1B <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1B/GKHeso_R2_TASB_seq.txt", header = F, sep = "\t") %>%
  mutate(genotype = rep("heso1")) %>%
  mutate(rep = rep("R2")) %>%
  mutate(TAS = rep("TAS1B"))

GKheso1_R3_TAS1B <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1B/GKHeso_R3_TASB_seq.txt", header = F, sep = "\t") %>%
  mutate(genotype = rep("heso1")) %>%
  mutate(rep = rep("R3")) %>%
  mutate(TAS = rep("TAS1B"))

heso1urt1_R1_TAS1B <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1B/heso1urt1_R1_TASB_seq.txt", header = F, sep = "\t") %>%
  mutate(genotype = rep("heso1urt1")) %>%
  mutate(rep = rep("R1")) %>%
  mutate(TAS = rep("TAS1B"))

heso1urt1_R2_TAS1B <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1B/heso1urt1_R2_TASB_seq.txt", header = F, sep = "\t") %>%
  mutate(genotype = rep("heso1urt1")) %>%
  mutate(rep = rep("R2")) %>%
  mutate(TAS = rep("TAS1B"))

heso1urt1_R3_TAS1B <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1B/heso1urt1_R3_TASB_seq.txt", header = F, sep = "\t") %>%
  mutate(genotype = rep("heso1urt1")) %>%
  mutate(rep = rep("R3")) %>%
  mutate(TAS = rep("TAS1B"))

urt1.1_R1_TAS1B <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1B/urt1_R1_TASB_seq.txt", header = F, sep = "\t") %>%
  mutate(genotype = rep("urt1")) %>%
  mutate(rep = rep("R1")) %>%
  mutate(TAS = rep("TAS1B"))

urt1.1_R2_TAS1B <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1B/urt1_R2_TASB_seq.txt", header = F, sep = "\t") %>%
  mutate(genotype = rep("urt1")) %>%
  mutate(rep = rep("R2")) %>%
  mutate(TAS = rep("TAS1B"))

urt1.1_R3_TAS1B <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1B/urt1_R3_TASB_seq.txt", header = F, sep = "\t") %>%
  mutate(genotype = rep("urt1")) %>%
  mutate(rep = rep("R3")) %>%
  mutate(TAS = rep("TAS1B"))

TAS1B_dfs <- ls(pattern = "_TAS1B")

# Initialize an empty data frame to store the merged data
TAS1B_merged_df <- data.frame()

# Loop through each data frame and merge it using rbind
for (df_name in TAS1B_dfs) {
  df <- get(df_name)  # Get the data frame by name
  TAS1B_merged_df <- rbind(TAS1B_merged_df, df)  # Merge with rbind
}
write.xlsx(TAS1B_merged_df, file = "/isdata/PBgrp/xpj980/datascratch/triplemutants/TAS1B.xlsx")

#### TAS1A ####

wt_R1_TAS1A <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1A/wt_R1_TASA_seq.txt", header = F, sep = "\t") %>%
  mutate(genotype = rep("wt")) %>%
  mutate(rep = rep("R1")) %>%
  mutate(TAS = rep("TAS1A"))

wt_R2_TAS1A <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1A/wt_R2_TASA_seq.txt", header = F, sep = "\t") %>%
  mutate(genotype = rep("wt")) %>%
  mutate(rep = rep("R2")) %>%
  mutate(TAS = rep("TAS1A"))

wt_R3_TAS1A <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1A/wt_R3_TASA_seq.txt", header = F, sep = "\t") %>%
  mutate(genotype = rep("wt")) %>%
  mutate(rep = rep("R3")) %>%
  mutate(TAS = rep("TAS1A"))

GKheso1_R1_TAS1A <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1A/GKheso_R1_TASA_seq.txt", header = F, sep = "\t") %>%
  mutate(genotype = rep("heso1")) %>%
  mutate(rep = rep("R1")) %>%
  mutate(TAS = rep("TAS1A"))

GKheso1_R2_TAS1A <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1A/GKheso_R2_TASA_seq.txt", header = F, sep = "\t") %>%
  mutate(genotype = rep("heso1")) %>%
  mutate(rep = rep("R2")) %>%
  mutate(TAS = rep("TAS1A"))

GKheso1_R3_TAS1A <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1A/GKheso_R3_TASA_seq.txt", header = F, sep = "\t") %>%
  mutate(genotype = rep("heso1")) %>%
  mutate(rep = rep("R3")) %>%
  mutate(TAS = rep("TAS1A"))

heso1urt1_R1_TAS1A <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1A/heso1urt1_R1_TASA_seq.txt", header = F, sep = "\t") %>%
  mutate(genotype = rep("heso1urt1")) %>%
  mutate(rep = rep("R1")) %>%
  mutate(TAS = rep("TAS1A"))

heso1urt1_R2_TAS1A <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1A/heso1urt1_R2_TASA_seq.txt", header = F, sep = "\t") %>%
  mutate(genotype = rep("heso1urt1")) %>%
  mutate(rep = rep("R2")) %>%
  mutate(TAS = rep("TAS1A"))

heso1urt1_R3_TAS1A <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1A/heso1urt1_R3_TASA_seq.txt", header = F, sep = "\t") %>%
  mutate(genotype = rep("heso1urt1")) %>%
  mutate(rep = rep("R3")) %>%
  mutate(TAS = rep("TAS1A"))

urt1.1_R1_TAS1A <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1A/urt1_R1_TASA_seq.txt", header = F, sep = "\t") %>%
  mutate(genotype = rep("urt1")) %>%
  mutate(rep = rep("R1")) %>%
  mutate(TAS = rep("TAS1A"))

urt1.1_R2_TAS1A <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1A/urt1_R2_TASA_seq.txt", header = F, sep = "\t") %>%
  mutate(genotype = rep("urt1")) %>%
  mutate(rep = rep("R2")) %>%
  mutate(TAS = rep("TAS1A"))

urt1.1_R3_TAS1A <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS1A/urt1_R3_TASA_seq.txt", header = F, sep = "\t") %>%
  mutate(genotype = rep("urt1")) %>%
  mutate(rep = rep("R3")) %>%
  mutate(TAS = rep("TAS1A"))

TAS1A_dfs <- ls(pattern = "_TAS1A")

# Initialize an empty data frame to store the merged data
TAS1A_merged_df <- data.frame()

# Loop through each data frame and merge it using rbind
for (df_name in TAS1A_dfs) {
  df <- get(df_name)  # Get the data frame by name
  TAS1A_merged_df <- rbind(TAS1A_merged_df, df)  # Merge with rbind
}
write.xlsx(TAS1A_merged_df, file = "/isdata/PBgrp/xpj980/datascratch/triplemutants/TAS1A.xlsx")

#### TAS2 ####

wt_R1_TAS2 <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS2/wt_R1_TAS2_seq.txt", header = F, sep = "\t") %>%
  mutate(genotype = rep("wt")) %>%
  mutate(rep = rep("R1")) %>%
  mutate(TAS = rep("TAS2"))

wt_R2_TAS2 <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS2/wt_R2_TAS2_seq.txt", header = F, sep = "\t") %>%
  mutate(genotype = rep("wt")) %>%
  mutate(rep = rep("R2")) %>%
  mutate(TAS = rep("TAS2"))

wt_R3_TAS2 <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS2/wt_R3_TAS2_seq.txt", header = F, sep = "\t") %>%
  mutate(genotype = rep("wt")) %>%
  mutate(rep = rep("R3")) %>%
  mutate(TAS = rep("TAS2"))

GKheso1_R1_TAS2 <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS2/GKheso_R1_TAS2_seq.txt", header = F, sep = "\t") %>%
  mutate(genotype = rep("heso1")) %>%
  mutate(rep = rep("R1")) %>%
  mutate(TAS = rep("TAS2"))

GKheso1_R2_TAS2 <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS2/GKheso_R2_TAS2_seq.txt", header = F, sep = "\t") %>%
  mutate(genotype = rep("heso1")) %>%
  mutate(rep = rep("R2")) %>%
  mutate(TAS = rep("TAS2"))

GKheso1_R3_TAS2 <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS2/GKheso_R3_TAS2_seq.txt", header = F, sep = "\t") %>%
  mutate(genotype = rep("heso1")) %>%
  mutate(rep = rep("R3")) %>%
  mutate(TAS = rep("TAS2"))

heso1urt1_R1_TAS2 <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS2/heso1urt1_R1_TAS2_seq.txt", header = F, sep = "\t") %>%
  mutate(genotype = rep("heso1urt1")) %>%
  mutate(rep = rep("R1")) %>%
  mutate(TAS = rep("TAS2"))

heso1urt1_R2_TAS2 <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS2/heso1urt1_R2_TAS2_seq.txt", header = F, sep = "\t") %>%
  mutate(genotype = rep("heso1urt1")) %>%
  mutate(rep = rep("R2")) %>%
  mutate(TAS = rep("TAS2"))

heso1urt1_R3_TAS2 <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS2/heso1urt1_R3_TAS2_seq.txt", header = F, sep = "\t") %>%
  mutate(genotype = rep("heso1urt1")) %>%
  mutate(rep = rep("R3")) %>%
  mutate(TAS = rep("TAS2"))

urt1.1_R1_TAS2 <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS2/urt1_R1_TAS2_seq.txt", header = F, sep = "\t") %>%
  mutate(genotype = rep("urt1")) %>%
  mutate(rep = rep("R1")) %>%
  mutate(TAS = rep("TAS2"))

urt1.1_R2_TAS2 <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS2/urt1_R2_TAS2_seq.txt", header = F, sep = "\t") %>%
  mutate(genotype = rep("urt1")) %>%
  mutate(rep = rep("R2")) %>%
  mutate(TAS = rep("TAS2"))

urt1.1_R3_TAS2 <- read.table("/isdata/PBgrp/xpj980/datascratch/triplemutants/03.trimmed/txt_files/TAS2/urt1_R3_TAS2_seq.txt", header = F, sep = "\t") %>%
  mutate(genotype = rep("urt1")) %>%
  mutate(rep = rep("R3")) %>%
  mutate(TAS = rep("TAS2"))

TAS2_dfs <- ls(pattern = "_TAS2")

# Initialize an empty data frame to store the merged data
TAS2_merged_df <- data.frame()

# Loop through each data frame and merge it using rbind
for (df_name in TAS2_dfs) {
  df <- get(df_name)  # Get the data frame by name
  TAS2_merged_df <- rbind(TAS2_merged_df, df)  # Merge with rbind
}
#write.xlsx(TAS2_merged_df, file = "/isdata/PBgrp/xpj980/datascratch/triplemutants/TAS2.xlsx")

TAS_merged <- rbind(TAS1A_merged_df, TAS1B_merged_df, TAS1C_merged_df, TAS2_merged_df)

TAS_summarized <- TAS_merged %>%
  group_by(genotype, rep, TAS) %>%
  dplyr::count() %>%
  ungroup() %>%
  add_row(genotype = "heso1", rep = "R1", TAS = "TAS1A", n = 0) %>%
  add_row(genotype = "heso1", rep = "R1", TAS = "TAS1B", n = 0) %>%
  arrange(., -n)
write.xlsx(TAS2_merged_df, file = "/isdata/PBgrp/xpj980/datascratch/triplemutants/TAS_summarized.xlsx")

write.xlsx(TAS_summarized, file = "/isdata/PBgrp/xpj980/datascratch/triplemutants/TAS_summarized_counts.xlsx")

# also put AFB2 rows

TAS_summarized <- TAS_merged %>%
  group_by(genotype, rep, TAS) %>%
  dplyr::count() %>%
  ungroup() %>%
  add_row(genotype = "heso1", rep = "R1", TAS = "TAS1A", n = 0) %>%
  add_row(genotype = "heso1", rep = "R1", TAS = "TAS1B", n = 0) %>%
  add_row(genotype = "wt", rep = "R1", TAS = "AFB2", n = 0) %>%
  add_row(genotype = "wt", rep = "R2", TAS = "AFB2", n = 0) %>%
  add_row(genotype = "wt", rep = "R3", TAS = "AFB2", n = 0) %>%
  add_row(genotype = "urt1", rep = "R1", TAS = "AFB2", n = 0) %>%
  add_row(genotype = "urt1", rep = "R2", TAS = "AFB2", n = 0) %>%
  add_row(genotype = "urt1", rep = "R3", TAS = "AFB2", n = 2) %>%
  add_row(genotype = "heso1", rep = "R1", TAS = "AFB2", n = 0) %>%
  add_row(genotype = "heso1", rep = "R2", TAS = "AFB2", n = 0) %>%
  add_row(genotype = "heso1", rep = "R3", TAS = "AFB2", n = 1) %>%
  add_row(genotype = "heso1urt1", rep = "R1", TAS = "AFB2", n = 2) %>%
  add_row(genotype = "heso1urt1", rep = "R2", TAS = "AFB2", n = 0) %>%
  add_row(genotype = "heso1urt1", rep = "R3", TAS = "AFB2", n = 0) %>%
  arrange(., -n)

TAS_summarized$genotype = factor(TAS_summarized$genotype, levels = c("wt", "urt1", "heso1", "heso1urt1"))    

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#### EV6A ####
as_plot <- TAS_summarized %>%
  ggplot(., aes(x = genotype, y = n, fill = TAS)) + 
  stat_summary(geom = "bar", fun = "mean", position = "dodge") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position = "dodge", size = 0.2) +
  #facet_wrap(~TAS) +
  #scale_x_continuous(breaks = c(19:24)) +
  scale_fill_manual(values = cbPalette) + 
  cowplot::theme_cowplot()
#ggsave(as_plot, filename = "/isdata/PBgrp/xpj980/figures/triplemutants/TASreadsTableasPlot.svg", device = "svg", units = "cm", width = 15, height = 15, dpi = 300)

