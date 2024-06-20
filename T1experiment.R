setwd("/binf-isilon/PBgrp/xpj980/datascratch/triplemutants/")

library(magrittr)
#install.packages("tidyverse")
library(tidyverse)
library(reshape2)
library(stringr)
library(Biostrings)
library(rtracklayer)
library(GenomicRanges)
library(devtools)
library(RColorBrewer)
library(forcats)
library(readxl)
#install.packages("seqinr")
library(seqinr)
library(stringi)
library(openxlsx)
#install.packages("wesanderson")
library(wesanderson)
#install.packages("cowplot")
library(cowplot)

##### T1 experiment ####
#### performed twice (count 1 or 2)

number <- c(11,175,8,253,67,30,57,48,47,68,61,102)
phenotype <- c("sick", "healthy", "sick", "healthy", "sick", "healthy", "sick", "healthy", "sick", "healthy", "sick", "healthy")
genotype <- c("empty", "empty","empty", "empty","HESO1WT", "HESO1WT","HESO1WT", "HESO1WT", "HESO1DADA", "HESO1DADA","HESO1DADA", "HESO1DADA")
count <- c(1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 2, 2)

test <- data_frame(number, phenotype, genotype, count)

fraction <- test %>% 
  group_by(genotype, count) %>% 
  mutate(total = sum(number)) %>% 
  mutate(percentage = (number/total)*100)

fraction$genotype <- factor(fraction$genotype, levels = c("empty", "HESO1WT", "HESO1DADA"))
fraction$phenotype <- factor(fraction$phenotype, levels = c("sick", "healthy"))


(plot <- ggplot(fraction, aes(x = genotype, y = percentage, fill = phenotype)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = "genotype", y = "% of phenotype", fill = "phenotype") +
    facet_wrap(~count) + 
    cowplot::theme_cowplot() +
    scale_fill_manual(values = c("#712B75", "#019267")) +
    ggtitle("T1 experiment - reintroducing HESO1 in rrp4urt1heso1")) 
#ggsave(plot, file = "/binf-isilon/PBgrp/xpj980/figures/triplemutants/T1experiment.svg", units = "cm", width = 32, height = 9, dpi = 300)

rep2 <- fraction %>% 
  ungroup() %>%
  filter(., count == 2) %>% 
  select(number, phenotype, genotype) %>%
  mutate(count = rep(1))

chisq.test(matrix(c(11,67,47,175,30,68), 3,2), simulate.p.value = T)


#### new T1 experiment (more phenotypes) #####

## anthocyanin

number <- c(24,106,4,101,2,137)
phenotype <- c("sick", "healthy","sick", "healthy","sick", "healthy")
genotype <- c("HESO1WT", "HESO1WT", "HESO1DADA", "HESO1DADA","empty", "empty")

test <- data_frame(number, phenotype, genotype)

fraction <- test %>% 
  group_by(genotype) %>% 
  mutate(total = sum(number)) %>% 
  mutate(percentage = (number/total)*100)

fraction$genotype <- factor(fraction$genotype, levels = c("empty", "HESO1WT", "HESO1DADA"))
fraction$phenotype <- factor(fraction$phenotype, levels = c("sick", "healthy"))


(plot <- ggplot(fraction, aes(x = genotype, y = percentage, fill = phenotype)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = "genotype", y = "% of phenotype", fill = "phenotype") +
    cowplot::theme_cowplot() +
    scale_fill_manual(values = c("#712B75", "#019267")) +
    ggtitle("T1 experiment - reintroducing HESO1 in rrp4urt1heso1")) 
ggsave(plot, file = "/binf-isilon/PBgrp/xpj980/figures/triplemutants/T1experiment_anthocyanin.svg", units = "cm", width = 32, height = 9, dpi = 300)

## bleach

number <- c(67,63,67,38,55,84)
phenotype <- c("sick", "healthy","sick", "healthy","sick", "healthy")
genotype <- c("HESO1WT", "HESO1WT", "HESO1DADA", "HESO1DADA","empty", "empty")

test <- data_frame(number, phenotype, genotype)

fraction <- test %>% 
  group_by(genotype) %>% 
  mutate(total = sum(number)) %>% 
  mutate(percentage = (number/total)*100)

fraction$genotype <- factor(fraction$genotype, levels = c("empty", "HESO1WT", "HESO1DADA"))
fraction$phenotype <- factor(fraction$phenotype, levels = c("sick", "healthy"))


(plot <- ggplot(fraction, aes(x = genotype, y = percentage, fill = phenotype)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = "genotype", y = "% of phenotype", fill = "phenotype") +
    cowplot::theme_cowplot() +
    scale_fill_manual(values = c("#E6B325", "#019267")) +
    ggtitle("T1 experiment - reintroducing HESO1 in rrp4urt1heso1")) 
ggsave(plot, file = "/binf-isilon/PBgrp/xpj980/figures/triplemutants/T1experiment_bleached.svg", units = "cm", width = 32, height = 9, dpi = 300)

## shoes 

number <- c(32,98,15,90,40,99)
phenotype <- c("sick", "healthy","sick", "healthy","sick", "healthy")
genotype <- c("HESO1WT", "HESO1WT", "HESO1DADA", "HESO1DADA","empty", "empty")

test <- data_frame(number, phenotype, genotype)

fraction <- test %>% 
  group_by(genotype) %>% 
  mutate(total = sum(number)) %>% 
  mutate(percentage = (number/total)*100)

fraction$genotype <- factor(fraction$genotype, levels = c("empty", "HESO1WT", "HESO1DADA"))
fraction$phenotype <- factor(fraction$phenotype, levels = c("sick", "healthy"))


(plot <- ggplot(fraction, aes(x = genotype, y = percentage, fill = phenotype)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = "genotype", y = "% of phenotype", fill = "phenotype") +
    cowplot::theme_cowplot() +
    scale_fill_manual(values = c("#355764", "#019267")) +
    ggtitle("T1 experiment - reintroducing HESO1 in rrp4urt1heso1")) 
ggsave(plot, file = "/binf-isilon/PBgrp/xpj980/figures/triplemutants/T1experiment_shoes.svg", units = "cm", width = 32, height = 9, dpi = 300)

## dwarf 

number <- c(103,27,87,18,67,72)
phenotype <- c("sick", "healthy","sick", "healthy","sick", "healthy")
genotype <- c("HESO1WT", "HESO1WT", "HESO1DADA", "HESO1DADA","empty", "empty")

test <- data_frame(number, phenotype, genotype)

fraction <- test %>% 
  group_by(genotype) %>% 
  mutate(total = sum(number)) %>% 
  mutate(percentage = (number/total)*100)

fraction$genotype <- factor(fraction$genotype, levels = c("empty", "HESO1WT", "HESO1DADA"))
fraction$phenotype <- factor(fraction$phenotype, levels = c("sick", "healthy"))


(plot <- ggplot(fraction, aes(x = genotype, y = percentage, fill = phenotype)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = "genotype", y = "% of phenotype", fill = "phenotype") +
    cowplot::theme_cowplot() +
    scale_fill_manual(values = c("#0F3D3E", "#019267")) +
    ggtitle("T1 experiment - reintroducing HESO1 in rrp4urt1heso1")) 
ggsave(plot, file = "/binf-isilon/PBgrp/xpj980/figures/triplemutants/T1experiment_dwarfs.svg", units = "cm", width = 32, height = 9, dpi = 300)

#### new T1 experiment (more phenotypes) DAG 38 - same plants but counted later 

## anthocyanin

number <- c(58,57,24,62,7,129)
phenotype <- c("sick", "healthy","sick", "healthy","sick", "healthy")
genotype <- c("HESO1WT", "HESO1WT", "HESO1DADA", "HESO1DADA","empty", "empty")

test <- data_frame(number, phenotype, genotype)

fraction <- test %>% 
  group_by(genotype) %>% 
  mutate(total = sum(number)) %>% 
  mutate(percentage = (number/total)*100)

fraction$genotype <- factor(fraction$genotype, levels = c("empty", "HESO1WT", "HESO1DADA"))
fraction$phenotype <- factor(fraction$phenotype, levels = c("sick", "healthy"))


(plot <- ggplot(fraction, aes(x = genotype, y = percentage, fill = phenotype)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = "genotype", y = "% of phenotype", fill = "phenotype") +
    cowplot::theme_cowplot() +
    scale_fill_manual(values = c("#712B75", "#019267")) +
    ggtitle("T1 experiment - reintroducing HESO1 in rrp4urt1heso1")) 
ggsave(plot, file = "/binf-isilon/PBgrp/xpj980/figures/triplemutants/T1experiment_anthocyanin_DAG38.svg", units = "cm", width = 32, height = 9, dpi = 300)


data <- matrix(c(57, 58, 129, 7), nrow = 2, byrow = TRUE)

# Name the rows and columns for clarity
rownames(data) <- c("Heso1", "WT")
colnames(data) <- c("Healthy", "Sick")

# Print the contingency table
print(data)

# Perform the chi-square test
chi_square_test <- chisq.test(data)


## bleach

number <- c(51,64,42,44,36,100)
phenotype <- c("sick", "healthy","sick", "healthy","sick", "healthy")
genotype <- c("HESO1WT", "HESO1WT", "HESO1DADA", "HESO1DADA","empty", "empty")

test <- data_frame(number, phenotype, genotype)

fraction <- test %>% 
  group_by(genotype) %>% 
  mutate(total = sum(number)) %>% 
  mutate(percentage = (number/total)*100)

fraction$genotype <- factor(fraction$genotype, levels = c("empty", "HESO1WT", "HESO1DADA"))
fraction$phenotype <- factor(fraction$phenotype, levels = c("sick", "healthy"))


(plot <- ggplot(fraction, aes(x = genotype, y = percentage, fill = phenotype)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = "genotype", y = "% of phenotype", fill = "phenotype") +
    cowplot::theme_cowplot() +
    scale_fill_manual(values = c("#E6B325", "#019267")) +
    ggtitle("T1 experiment - reintroducing HESO1 in rrp4urt1heso1")) 
ggsave(plot, file = "/binf-isilon/PBgrp/xpj980/figures/triplemutants/T1experiment_bleached_DAG38.svg", units = "cm", width = 32, height = 9, dpi = 300)

## shoes 

number <- c(24,91,17,69,31,105)
phenotype <- c("sick", "healthy","sick", "healthy","sick", "healthy")
genotype <- c("HESO1WT", "HESO1WT", "HESO1DADA", "HESO1DADA","empty", "empty")

test <- data_frame(number, phenotype, genotype)

fraction <- test %>% 
  group_by(genotype) %>% 
  mutate(total = sum(number)) %>% 
  mutate(percentage = (number/total)*100)

fraction$genotype <- factor(fraction$genotype, levels = c("empty", "HESO1WT", "HESO1DADA"))
fraction$phenotype <- factor(fraction$phenotype, levels = c("sick", "healthy"))


(plot <- ggplot(fraction, aes(x = genotype, y = percentage, fill = phenotype)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = "genotype", y = "% of phenotype", fill = "phenotype") +
    cowplot::theme_cowplot() +
    scale_fill_manual(values = c("#355764", "#019267")) +
    ggtitle("T1 experiment - reintroducing HESO1 in rrp4urt1heso1")) 
ggsave(plot, file = "/binf-isilon/PBgrp/xpj980/figures/triplemutants/T1experiment_shoes_DAG38.svg", units = "cm", width = 32, height = 9, dpi = 300)

## dwarf 

number <- c(79,36,57,29,36,100)
phenotype <- c("sick", "healthy","sick", "healthy","sick", "healthy")
genotype <- c("HESO1WT", "HESO1WT", "HESO1DADA", "HESO1DADA","empty", "empty")

test <- data_frame(number, phenotype, genotype)

fraction <- test %>% 
  group_by(genotype) %>% 
  mutate(total = sum(number)) %>% 
  mutate(percentage = (number/total)*100)

fraction$genotype <- factor(fraction$genotype, levels = c("empty", "HESO1WT", "HESO1DADA"))
fraction$phenotype <- factor(fraction$phenotype, levels = c("sick", "healthy"))


(plot <- ggplot(fraction, aes(x = genotype, y = percentage, fill = phenotype)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = "genotype", y = "% of phenotype", fill = "phenotype") +
    cowplot::theme_cowplot() +
    scale_fill_manual(values = c("#0F3D3E", "#019267")) +
    ggtitle("T1 experiment - reintroducing HESO1 in rrp4urt1heso1")) 
ggsave(plot, file = "/binf-isilon/PBgrp/xpj980/figures/triplemutants/T1experiment_dwarfs_DAG38.svg", units = "cm", width = 32, height = 9, dpi = 300)

data <- matrix(c(36, 79, 100, 36), nrow = 2, byrow = TRUE)

#### STATISTICS FOR PAPER ####
# Name the rows and columns for clarity
rownames(data) <- c("Heso1", "WT")
colnames(data) <- c("Healthy", "Sick")

# Print the contingency table
print(data)

# Perform the chi-square test
chi_square_test <- chisq.test(data)

