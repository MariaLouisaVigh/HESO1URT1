setwd("/binf-isilon/PBgrp/xpj980/datascratch/3RACE/RACEseq3/")

args = commandArgs(trailingOnly=TRUE)

#### load libraries ####
library(dplyr)
library(stringr)
library(purrr)
library(readr)
library(tidyverse)
library(wesanderson)
#install.packages("spgs")
library(spgs)
#install.packages("ggridges")
library(ggridges)
library(openxlsx)

is_similar <- function(r1,
                       r2,
                       max_diff=2) {
  
  r1s<-str_split(r1,"", simplify = T)
  r2s<-str_split(r2,"", simplify = T)
  sum(r1s != r2s)<=max_diff
}


similar_to_any <- function(bar1, barcodes, max_diff=2) {
  sum <- map_lgl(barcodes, function(bar2) {
    is_similar(bar1,bar2, max_diff = max_diff)
  }) %>% sum()
  sum>1
}

similar_to_any2 <- function(barcodes, max_diff=2) {
  
  map_lgl(barcodes, similar_to_any, barcodes=barcodes, max_diff=max_diff)
  
}


how_many_similar <- function(barcodes) {
  sum(similar_to_any2(barcodes))
}

mean_n_similar <- function(n,
                           reps=10,
                           size=15) {
  replicate(reps, {
    barcodes <- replicate(n, sample(c("A","T","C","G"),size = size, replace = T) %>% str_c(collapse = ""), simplify = T)
    how_many_similar(barcodes)
  }) %>% mean()
}

# test how many similar for n = 10 to 100 
# map_dbl(c(1:10)*10, mean_n_similar, size=15, reps=10)

####RACESeq set1 ####
#### AFB2 5CF###### 
# 
# file_names <- list.files("/binf-isilon/PBgrp/xpj980/datascratch/3RACE/02.merged_lanes/textfiles_deduplicator/AFB2_NIBBLE", full.names = T, pattern = "ends")
# 
# DF <- tibble()
# for (f in file_names) {
#   print(f)
#   dat <- read_lines(f) %>% 
#     as_tibble() %>%
#     mutate(name = f %>% str_replace(".*//(.*)_read2.*", "\\1"),
#            length=nchar(value))
#   DF <- rbind(DF, dat)
# }
# 
# df <- DF %>%
#   mutate(barcode=str_sub(value,1,15),
#          delim=str_sub(value,16,20),
#          read=str_sub(value,21,-1L)) %>% 
#   filter(length > 19)
# 
# deduplicated <- df %>%
#   group_by(read, name) %>%
#   mutate(is_dup=similar_to_any2(barcode, max_diff = 2)) %>%
#   arrange(read) %>%
#   group_by(read, name, is_dup) %>%
#   mutate(n=row_number()) %>% 
#   filter(n==1 | (!is_dup))
#save(deduplicated, file ="/binf-isilon/PBgrp/xpj980/datascratch/3RACE/deduplicatedDFfirstlib.Rdat")
load(file ="/binf-isilon/PBgrp/xpj980/datascratch/3RACE/deduplicatedDFfirstlib.Rdat")

deduplicated$name2 <- str_remove(deduplicated$name, "/binf-isilon/PBgrp/xpj980/datascratch/3RACE/02.merged_lanes/textfiles_deduplicator/AFB2_NIBBLE/") 

deduplicated$name <- str_remove(deduplicated$name2, "_read2.ends.txt") 

get_tail <- function(x) {
  tail <- str_remove(x, "AACAATGCGA|AACAATGCG|AACAATGC|AACAATG|AACAAT|AACAA|AACA|AAC|AA|A|")
  trim <- -(10-nchar(x)+nchar(tail))
  return(list(tail = tail, trim = trim))
}

DF2 <- deduplicated %>% 
  ungroup() %>%
  mutate(revcomp = reverseComplement(deduplicated$read, case = "upper"), 
         tail = get_tail(revcomp)$tail, 
         trim = get_tail(revcomp)$trim, 
         tail_len = nchar(tail)) %>% 
  separate(name, c("sample", "rep"), sep="_")

rem_seq_errors <- DF2 %>% 
  mutate(to_discard = case_when(-(trim) == tail_len ~"discard", 
                                T ~ "keep"), 
         to_discard = case_when(trim == 0 & tail_len == 0 ~ "keep", 
                                T ~ to_discard), 
         tail_len = as.double(tail_len)) %>% 
  filter(to_discard != "discard") 

### to gather with newer replicates ###

rem_seq_errors$lib <- paste0("lib1")  

rem_seq_errorslib1 <- rem_seq_errors

add_pos <- rem_seq_errors %>%
  mutate(position = case_when(trim < 0 & tail_len == 0 ~ trim, 
                              T ~ tail_len))
to_plot <- add_pos %>% 
  group_by(sample, rep, position) %>% 
  filter(tail_len <20) %>%
  summarize(count=n()) %>% 
  ungroup() %>% 
  group_by(sample, rep) %>% 
  mutate(ncount=count/sum(count)) 

to_plot2 <- rem_seq_errors %>% 
  group_by(sample, rep, trim, tail_len) %>% 
  filter(tail_len <10) %>%
  summarize(count=n()) %>% 
  ungroup() %>% 
  group_by(sample, rep) %>% 
  mutate(ncount=count/sum(count)) 

### BLAKE MEYERS PLOT ###

to_plot2 <- rem_seq_errors %>% 
  group_by(sample, trim, tail_len) %>% 
  filter(tail_len <10) %>%
  summarize(count=n()) %>% 
  ungroup() %>% 
  group_by(sample) %>% 
  mutate(ncount=count/sum(count)) 

to_plot2$sample = factor(to_plot2$sample, levels = c("WT", "rdr6", "urt1", "heso1", "heso1urt1"))

blake_plot_avg <- ggplot(to_plot2, aes(x = trim, y = tail_len)) + 
  geom_point(aes(color = sample, size = ncount)) +
  scale_colour_manual(values = rev(wes_palette("Darjeeling2", n = 5))) +
  #scale_colour_manual(values=c("rep1"="chocolate3","rep2"="cyan4", "rep3" = "firebrick")) + 
  scale_size_area(limits = c(0,1), breaks = c(0,0.001,0.01,0.1,1)) +
  facet_wrap(~sample, nrow = 1) +  
  scale_x_continuous(breaks = c(-10:0)) + 
  scale_y_continuous(breaks = c(0:10)) + 
  cowplot::theme_cowplot() +
  cowplot::background_grid() +
  ggtitle("Nibbling and tailing AFB2 5'CF") +
  xlab("nibbling") +
  ylab("tailing") 
#ggsave(blake_plot_avg, file = "/binf-isilon/PBgrp/xpj980/datascratch/3RACE/05.figures/BlakeMeyersavg_DEDUP.svg", units = "cm", width = 32, height = 9, dpi = 300)

# keep reps seperated 

to_plot2 <- rem_seq_errors %>% 
  group_by(sample, rep, trim, tail_len) %>% 
  filter(tail_len <10) %>%
  summarize(count=n()) %>% 
  ungroup() %>% 
  group_by(sample, rep) %>% 
  mutate(ncount=count/sum(count)) 

to_plot2$sample = factor(to_plot2$sample, levels = c("WT", "rdr6", "urt1", "heso1", "heso1urt1"))

blake_plot_avg <- ggplot(to_plot2, aes(x = trim, y = tail_len)) + 
  geom_point(aes(color = sample, size = ncount)) +
  scale_colour_manual(values = rev(wes_palette("Darjeeling2", n = 5))) +
  #scale_colour_manual(values=c("rep1"="chocolate3","rep2"="cyan4", "rep3" = "firebrick")) + 
  scale_size_area(limits = c(0,1), breaks = c(0,0.001,0.01,0.1,1)) +
  #facet_wrap(~sample, nrow = 1) +  
  facet_grid(sample ~ rep) +
  scale_x_continuous(breaks = c(-10:0)) + 
  scale_y_continuous(breaks = c(0:10)) + 
  cowplot::theme_cowplot() +
  cowplot::background_grid() +
  ggtitle("Nibbling and tailing AFB2 5'CF") +
  xlab("nibbling") +
  ylab("tailing") 
#ggsave(blake_plot_avg, file = "/binf-isilon/PBgrp/xpj980/datascratch/3RACE/05.figures/BlakeMeyersavg_DEDUP_reps.svg", units = "cm", width = 32, height = 9, dpi = 300)

how_many <- rem_seq_errors %>% 
  group_by(sample) %>% 
  summarise(n = n())

how_many_reps <- rem_seq_errors %>% 
  group_by(sample, rep) %>% 
  summarise(n = n())

## U-tail only 

Utail <- c("T", "TT", "TTT", "TTTT", "TTTTT", "TTTTTT", "TTTTTTT", "TTTTTTTT", "TTTTTTTTT", "TTTTTTTTTT")

to_plot3 <- rem_seq_errors %>% 
  filter(tail_len < 10) %>% 
  filter(tail %in% c("", Utail)) %>%
  group_by(sample, trim, tail_len) %>%
  summarize(count=n()) %>% 
  ungroup() %>% 
  group_by(sample) %>% 
  mutate(ncount=count/sum(count)) 

to_plot3$sample = factor(to_plot3$sample, levels = c("WT", "rdr6", "urt1", "heso1", "heso1urt1"))

blake_plot_onlyU <- ggplot(to_plot3, aes(x = trim, y = tail_len)) + 
  geom_point(aes(color = sample, size = ncount)) +
  scale_colour_manual(values = rev(wes_palette("Darjeeling2", n = 5))) +
  #scale_colour_manual(values=c("rep1"="chocolate3","rep2"="cyan4", "rep3" = "firebrick")) + 
  scale_size_area(limits = c(0,1), breaks = c(0,0.001,0.01,0.1,1)) +
  facet_wrap(~sample, nrow = 1) +  
  scale_x_continuous(breaks = c(-10:0)) + 
  scale_y_continuous(breaks = c(0:10)) + 
  cowplot::theme_cowplot() +
  cowplot::background_grid() +
  ggtitle("Nibbling and U-tailing AFB2 5'CF") +
  xlab("nibbling") +
  ylab("tailing") 
#ggsave(blake_plot_onlyU, file = "/binf-isilon/PBgrp/xpj980/datascratch/3RACE/05.figures/BlakeMeyersavgUtail_DEDUP.svg", units = "cm", width = 32, height = 9, dpi = 300)

how_many <- rem_seq_errors %>% 
  filter(tail_len < 10) %>% 
  filter(tail %in% c("", Utail)) %>%
  group_by(sample) %>% 
  summarise(n = n()) 

how_many_reps <- rem_seq_errors %>% 
  filter(tail_len < 10) %>% 
  filter(tail %in% c("", Utail)) %>%
  group_by(sample, rep) %>% 
  summarise(n = n()) 

# keep replicates seperate and only look at U tails 
to_plot3 <- rem_seq_errors %>% 
  filter(tail_len < 10) %>% 
  filter(tail %in% c("", Utail)) %>%
  group_by(sample, rep, trim, tail_len) %>%
  summarize(count=n()) %>% 
  ungroup() %>% 
  group_by(sample, rep) %>% 
  mutate(ncount=count/sum(count)) 

to_plot3$sample = factor(to_plot3$sample, levels = c("WT", "rdr6", "urt1", "heso1", "heso1urt1"))

blake_plot_onlyU <- ggplot(to_plot3, aes(x = trim, y = tail_len)) + 
  geom_point(aes(color = sample, size = ncount)) +
  scale_colour_manual(values = rev(wes_palette("Darjeeling2", n = 5))) +
  #scale_colour_manual(values=c("rep1"="chocolate3","rep2"="cyan4", "rep3" = "firebrick")) + 
  scale_size_area(limits = c(0,1), breaks = c(0,0.001,0.01,0.1,1)) +
  facet_grid(rep~sample) +  
  scale_x_continuous(breaks = c(-10:0)) + 
  scale_y_continuous(breaks = c(0:10)) + 
  cowplot::theme_cowplot() +
  cowplot::background_grid() +
  ggtitle("Nibbling and U-tailing AFB2 5'CF") +
  xlab("nibbling") +
  ylab("tailing") 
#ggsave(blake_plot_onlyU, filename = "/binf-isilon/PBgrp/xpj980/datascratch/3RACE/05.figures/Meyerplot_reps_DEDUP.svg", units = "cm", width = 32, height = 9, dpi = 300)

## the tail should at least have 50% of same base to be a true tail 

test <- DF2 %>% 
  mutate(frac.t = str_count(DF2$tail, "T")/tail_len, 
         frac.c = str_count(DF2$tail, "C")/tail_len,
         frac.g = str_count(DF2$tail, "G")/tail_len,
         frac.a = str_count(DF2$tail, "A")/tail_len) %>% 
  filter(., tail == "" | frac.t > 0.5 | frac.c > 0.5 | frac.g > 0.5 | frac.a > 0.5)

how_many <- test %>% 
  group_by(sample) %>% 
  summarise(n = n())

how_many_reps <- test %>% 
  group_by(sample, rep) %>% 
  summarise(n = n())

test_to_plot <- test %>% 
  group_by(sample, trim, tail_len) %>% 
  filter(tail_len <10) %>%
  summarize(count=n()) %>% 
  ungroup() %>% 
  group_by(sample) %>% 
  mutate(ncount=count/sum(count)) 

test_to_plot$sample = factor(test_to_plot$sample, levels = c("WT", "rdr6", "urt1", "heso1", "heso1urt1"))

blake_plot_avg_test <- ggplot(test_to_plot, aes(x = trim, y = tail_len)) + 
  geom_point(aes(color = sample, size = ncount)) +
  scale_colour_manual(values = rev(wes_palette("Darjeeling2", n = 5))) +
  #scale_colour_manual(values=c("rep1"="chocolate3","rep2"="cyan4", "rep3" = "firebrick")) + 
  scale_size_area(limits = c(0,1), breaks = c(0,0.001,0.01,0.1,1)) +
  facet_wrap(~sample, nrow = 1) +  
  scale_x_continuous(breaks = c(-10:0)) + 
  scale_y_continuous(breaks = c(0:10)) + 
  cowplot::theme_cowplot() +
  cowplot::background_grid() +
  ggtitle("Nibbling and tailing AFB2 5'CF") +
  xlab("nibbling") +
  ylab("tailing") 
#ggsave(blake_plot_avg_test, file = "/binf-isilon/PBgrp/xpj980/datascratch/3RACE/05.figures/BlakeMeyersavgtailreq_DEDUP.svg", units = "cm", width = 32, height = 9, dpi = 300)

### Blake Meyers plot zoomed out ### 
test_to_plot <- test %>% 
  group_by(sample, trim, tail_len) %>% 
  summarize(count=n()) %>% 
  ungroup() %>% 
  group_by(sample) %>% 
  mutate(ncount=count/sum(count)) 

test_to_plot$sample = factor(test_to_plot$sample, levels = c("WT", "rdr6", "urt1", "heso1", "heso1urt1"))

blake_plot_avg_test <- ggplot(test_to_plot, aes(x = trim, y = tail_len)) + 
  geom_point(aes(color = sample, size = ncount)) +
  scale_colour_manual(values = rev(wes_palette("Darjeeling2", n = 5))) +
  #scale_colour_manual(values=c("rep1"="chocolate3","rep2"="cyan4", "rep3" = "firebrick")) + 
  scale_size_area(limits = c(0,1), breaks = c(0,0.001,0.01,0.1,1)) +
  facet_wrap(~sample, nrow = 1) +  
  scale_x_continuous(breaks = c(-10:0)) + 
  cowplot::theme_cowplot() +
  cowplot::background_grid() +
  ggtitle("Nibbling and tailing AFB2 5'CF") +
  xlab("nibbling") +
  ylab("tailing") 
#ggsave(blake_plot_avg_test, file = "/binf-isilon/PBgrp/xpj980/datascratch/3RACE/05.figures/BlakeMeyersavgtailreqnozoom_DEDUP.svg", units = "cm", width = 32, height = 9, dpi = 300)

### TAIL COMPOSITION ###

test2 <- test %>% 
  filter(tail_len < 20) %>%
  mutate(frac.t = str_count(tail, "T"), 
         frac.c = str_count(tail, "C"),
         frac.g = str_count(tail, "G"),
         frac.a = str_count(tail, "A")) %>% 
  group_by(sample) %>% 
  summarise(Uridine=sum(frac.t), Cytosine=sum(frac.c), Guanosine=sum(frac.g), Adenosine=sum(frac.a)) %>% 
  gather(base, count, Uridine, Cytosine, Guanosine, Adenosine) %>%
  group_by(sample) %>% 
  mutate(total = sum(count), 
         perc = count/total, 
         Round_off = round(perc ,digit=2))

test2$sample = factor(test2$sample, levels = c("WT", "rdr6", "urt1", "heso1", "heso1urt1"))  
plot <- ggplot(test2, aes(x = factor(sample), y = Round_off*100, fill = factor(base))) +
  geom_bar(stat="identity", width = 0.7) +
  labs(x = "genotype", y = "percent", fill = "base") +
  cowplot::theme_cowplot() +
  geom_text(aes(label = paste0(Round_off*100,"%")), 
            position = position_stack(vjust = 0.5), size = 4) +
  scale_fill_manual(values = wes_palette("GrandBudapest1", n = 4)) +
  ggtitle("Tail composition of tails 1-20 nt long") +
  xlab("genotype") 
#ggsave(plot, file = "/binf-isilon/PBgrp/xpj980/datascratch/3RACE/05.figures/tailcomposition.svg", units = "in", width = 30, height = 30, dpi = 300)

### with plot showing how many % of total reads even have a tail 
test2 <- test %>% 
  group_by(sample, tail) %>%
  filter(tail_len < 20) %>% 
  summarise(count = n()) %>% 
  group_by(sample) %>% 
  mutate(ncount = count/sum(count)) %>% 
  filter(tail != "") %>% 
  mutate(frac.u = str_count(tail, "T")/nchar(tail), 
         frac.c = str_count(tail, "C")/nchar(tail),
         frac.g = str_count(tail, "G")/nchar(tail),
         frac.a = str_count(tail, "A")/nchar(tail)) %>% 
  gather(base, basecount, contains("frac"))

test2$sample = factor(test2$sample, levels = c("WT", "rdr6", "urt1", "heso1", "heso1urt1"))    
plot <- ggplot(test2, aes(x = factor(sample), y = (basecount*ncount)*100, fill = base %>% str_remove("frac.") %>% toupper())) +
  geom_bar(stat="identity", width = 0.7) +
  labs(x = "genotype", y = "% of reads with a tail", fill = "base") +
  cowplot::theme_cowplot() +
  cowplot::background_grid() +
  #geom_text(aes(label = paste0(Round_off,"%")), 
  #position = position_stack(vjust = 0.5), size = 4) +
  scale_fill_manual(values = wes_palette("GrandBudapest1", n = 4)) +
  ggtitle("Proportion of reads with a tail and composition of the tails (1-20 nt long)") 
#ggsave(plot, file = "/binf-isilon/PBgrp/xpj980/datascratch/3RACE/05.figures/test.svg", units = "in", width = 30, height = 30, dpi = 300)

# again, plot replicates individually 
test2 <- test %>% 
  group_by(rep, sample, tail) %>%
  filter(tail_len < 20) %>% 
  summarise(count = n()) %>% 
  group_by(rep, sample) %>% 
  mutate(ncount = count/sum(count)) %>% 
  filter(tail != "") %>% 
  mutate(frac.u = str_count(tail, "T")/nchar(tail), 
         frac.c = str_count(tail, "C")/nchar(tail),
         frac.g = str_count(tail, "G")/nchar(tail),
         frac.a = str_count(tail, "A")/nchar(tail)) %>% 
  gather(base, basecount, contains("frac"))

test2$sample = factor(test2$sample, levels = c("WT", "rdr6", "urt1", "heso1", "heso1urt1"))    
plot <- ggplot(test2, aes(x = factor(sample), y = (basecount*ncount)*100, fill = base %>% str_remove("frac.") %>% toupper())) +
  geom_bar(stat="identity", width = 0.7) +
  labs(x = "genotype", y = "% of reads with a tail", fill = "base") +
  cowplot::theme_cowplot() +
  cowplot::background_grid() +
  facet_wrap(~rep) +
  #geom_text(aes(label = paste0(Round_off,"%")), 
  #position = position_stack(vjust = 0.5), size = 4) +
  scale_fill_manual(values = wes_palette("GrandBudapest1", n = 4)) +
  ggtitle("Proportion of reads with a tail and composition of the tails (1-20 nt long)") 
#ggsave(plot, file = "/binf-isilon/PBgrp/xpj980/datascratch/3RACE/05.figures/test_reps_DEDUP.svg", units = "in", width = 30, height = 30, dpi = 300)

#### AFB2 polyA ####

file_names <- list.files("/binf-isilon/PBgrp/xpj980/datascratch/3RACE/02.merged_lanes/textfiles_deduplicator/AFB2polyA/", full.names = T, pattern = "ends")

DF <- tibble()
for (f in file_names) {
  print(f)
  dat <- read_lines(f) %>% 
    as_tibble() %>%
    mutate(name = f %>% str_replace(".*//(.*)_read2.*", "\\1"),
           length=nchar(value))
  DF <- rbind(DF, dat)
}

df <- DF %>%
  mutate(barcode=str_sub(value,1,15),
         delim=str_sub(value,16,20),
         read=str_sub(value,21,-1L)) %>%
  filter(length > 19)

deduplicated <- df %>%
  group_by(read, name) %>%
  mutate(is_dup=similar_to_any2(barcode, max_diff = 2)) %>%
  arrange(read) %>%
  group_by(read, name, is_dup) %>%
  mutate(n=row_number()) %>% 
  filter(n==1 | (!is_dup))

filter_tails <- deduplicated %>%
  ungroup() %>%
  separate(name, c("sample", "rep"), sep="_") %>%
  mutate(revcomp = reverseComplement(deduplicated$read, case = "upper"), 
         tail_len = nchar(revcomp))
tails_filtered <- filter_tails %>% 
  mutate(frac.t = str_count(filter_tails$revcomp, "T")/tail_len, 
         frac.c = str_count(filter_tails$revcomp, "C")/tail_len,
         frac.g = str_count(filter_tails$revcomp, "G")/tail_len, 
         frac.a = str_count(filter_tails$revcomp, "A")/tail_len) %>% 
  filter(., revcomp == "" | frac.t > 0.5 | frac.c > 0.5 | frac.g > 0.5 | frac.a > 0.5)

### base composition ###

## a distribution of the types of tails 
uendings <- tails_filtered %>% 
  mutate(tailtype = if_else(frac.a ==1, "Atail", "no"), 
         tailtype = if_else(str_detect(revcomp, "T$"), "Uending", tailtype),
         tailtype = if_else(revcomp == "", "tailless", tailtype),
         tailtype = case_when(tailtype=="no" ~"other", 
                              T ~ tailtype)) 
correctNA <- uendings %>% 
  mutate(tailtype = case_when(is.na(uendings$tailtype) ~ "tailless", 
                              T ~ tailtype)) 
#filter(tail_len < 40)

test <- correctNA %>% 
  group_by(sample, tailtype) %>% 
  summarise(count = n())   

plotting <- test %>% 
  group_by(sample) %>% 
  mutate(perc = count/sum(count)*100, 
         sample = factor(sample, c("WT", "rdr6", "heso1", "urt1", "heso1urt1")), 
         tailtype = factor(tailtype, c("tailless", "other", "Uending", "Atail"))) %>%
  ggplot(., aes(x = tailtype, y = perc, fill = tailtype)) + 
  geom_bar(stat = "identity", position = "dodge") +
  #geom_point() +
  facet_wrap(~sample, nrow = 5) +
  coord_flip() +
  #theme(aspect.ratio = 1) + 
  cowplot::theme_cowplot() +
  cowplot::background_grid() +
  ggtitle("Tailtype on AFB2 3'CF") +
  ylab("percentage of reads") +
  xlab("Type of tail") +
  scale_fill_manual(values = rev(wes_palette("Royal1", n = 4))) +
  theme(aspect.ratio = 1)
#cowplot::save_plot("/binf-isilon/PBgrp/xpj980/datascratch/3RACE/05.figures/tailtypeAFB23CF_DEDUP.svg", plotting, base_asp = 1, nrow = 5)

# and replicates seperated 
test <- correctNA %>% 
  group_by(rep, sample, tailtype) %>% 
  summarise(count = n())   

plotting <- test %>% 
  group_by(rep, sample) %>% 
  mutate(perc = count/sum(count)*100, 
         sample = factor(sample, c("WT", "rdr6", "heso1", "urt1", "heso1urt1")), 
         tailtype = factor(tailtype, c("tailless", "other", "Uending", "Atail"))) %>%
  ggplot(., aes(x = tailtype, y = perc, fill = tailtype)) + 
  geom_bar(stat = "identity", position = "dodge") +
  #geom_point() +
  facet_grid(rep~sample) +
  coord_flip() +
  #theme(aspect.ratio = 1) + 
  cowplot::theme_cowplot() +
  cowplot::background_grid() +
  ggtitle("Tailtype on AFB2 3'CF") +
  ylab("percentage of reads") +
  xlab("Type of tail") +
  scale_fill_manual(values = rev(wes_palette("Royal1", n = 4))) +
  theme(aspect.ratio = 1)
#cowplot::save_plot("/binf-isilon/PBgrp/xpj980/datascratch/3RACE/05.figures/tailtypeAFB23CF_reps_DEDUP.svg", plotting, base_asp = 1, nrow = 5)

# make errorbars instead 

test <- correctNA %>% 
  group_by(sample, rep, tailtype) %>% 
  summarise(count = n())   

plotting <- test %>% 
  group_by(sample, rep) %>% 
  mutate(perc = count/sum(count)*100, 
         sample = factor(sample, c("WT", "rdr6", "heso1", "urt1", "heso1urt1")), 
         tailtype = factor(tailtype, c("tailless", "other", "Uending", "Atail"))) %>%
  ggplot(., aes(x = tailtype, y = perc, fill = tailtype)) + 
  #geom_bar(stat = "identity", position = "dodge") +
  stat_summary(geom = "bar", fun = "mean", position = "dodge") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position = "dodge", size = 0.2) +
  #geom_point() +
  facet_wrap(~sample, nrow = 5) +
  coord_flip() +
  #theme(aspect.ratio = 1) + 
  cowplot::theme_cowplot() +
  cowplot::background_grid() +
  ggtitle("Tailtype on AFB2 3'CF") +
  ylab("percentage of reads") +
  xlab("Type of tail") +
  scale_fill_manual(values = rev(wes_palette("Royal1", n = 4))) +
  theme(aspect.ratio = 1)
#cowplot::save_plot("/binf-isilon/PBgrp/xpj980/datascratch/3RACE/05.figures/tailtypeAFB23CF_DEDUP_errorbars.svg", plotting, base_asp = 1, nrow = 5)

test <- correctNA %>% 
  group_by(sample, rep) %>% 
  summarise(count = n())  

### count Us in Uending tails 

correctNA$sample = factor(correctNA$sample, levels = c("WT", "rdr6", "urt1", "heso1", "heso1urt1"))    

look <- correctNA %>% 
  filter(tailtype == "Uending") %>%
  mutate(onlyU = str_extract(string = revcomp, pattern = regex("T+$")), 
         Uno = str_count(onlyU)) %>% 
  group_by(sample, rep, Uno) %>% 
  summarise(count = n()) %>% 
  group_by(sample, rep) %>% 
  mutate(perc = count/sum(count)*100) %>%
  ggplot(., aes(x = Uno, y = perc, fill = as.factor(Uno))) + 
  #geom_bar(stat = "identity", position = "dodge") + 
  stat_summary(geom = "bar", fun = "mean", position = "dodge") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position = "dodge", size = 0.2) +
  facet_wrap(~sample, nrow = 5) +
  scale_fill_manual(values = wes_palette("IsleofDogs1", n = 6)) +
  cowplot::theme_cowplot() +
  cowplot::background_grid() +
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks = c(1,2,3,4,5,6))
#ggsave(look, file = "/binf-isilon/PBgrp/xpj980/datascratch/3RACE/05.figures/AFB23CFUendingUnumber_DEDUP_errorbars.svg", units = "cm", width = 30, height = 15, dpi = 300)

look <- correctNA %>% 
  filter(tailtype == "Uending") %>%
  mutate(onlyU = str_extract(string = revcomp, pattern = regex("T+$")), 
         Uno = str_count(onlyU)) %>% 
  group_by(sample, rep, Uno) %>% 
  add_tally() %>% 
  ggplot(., aes(x = Uno, y = n, fill = as.factor(Uno))) + 
  geom_bar(stat = "identity", position = "dodge") + 
  facet_grid(rep~sample) +
  scale_fill_manual(values = wes_palette("IsleofDogs1", n = 6)) +
  cowplot::theme_cowplot() +
  cowplot::background_grid() +
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks = c(1,2,3,4,5,6))
#ggsave(look, file = "/binf-isilon/PBgrp/xpj980/datascratch/3RACE/05.figures/AFB23CFUendingUnumber_reps_sep.svg", units = "cm", width = 30, height = 15, dpi = 300)

#### RACESeq set 2 #####
#### AFB2 polyA MIR393 ####

file_names <- list.files("/binf-isilon/PBgrp/xpj980/datascratch/3RACE/RACEseq2/02.merged_lanes/txtfiles_deduplicator/polyA/", full.names = T, pattern = "ends")

DF <- tibble()
for (f in file_names) {
  print(f)
  dat <- read_lines(f) %>% 
    as_tibble() %>%
    mutate(name = f %>% str_replace(".*//(.*)_read2.*", "\\1"),
           length=nchar(value))
  DF <- rbind(DF, dat)
}

df <- DF %>%
  mutate(barcode=str_sub(value,1,15),
         delim=str_sub(value,16,20),
         read=str_sub(value,21,-1L)) %>%
  filter(length > 19)

deduplicated <- df %>%
  group_by(read, name) %>%
  mutate(is_dup=similar_to_any2(barcode, max_diff = 2)) %>%
  arrange(read) %>%
  group_by(read, name, is_dup) %>%
  mutate(n=row_number()) %>% 
  filter(n==1 | (!is_dup))

filter_tails <- deduplicated %>%
  ungroup() %>%
  separate(name, c("sample", "rep"), sep="_") %>%
  mutate(revcomp = reverseComplement(deduplicated$read, case = "upper"), 
         tail_len = nchar(revcomp))
tails_filtered <- filter_tails %>% 
  mutate(frac.t = str_count(filter_tails$revcomp, "T")/tail_len, 
         frac.c = str_count(filter_tails$revcomp, "C")/tail_len,
         frac.g = str_count(filter_tails$revcomp, "G")/tail_len, 
         frac.a = str_count(filter_tails$revcomp, "A")/tail_len) %>% 
  filter(., revcomp == "" | frac.t > 0.5 | frac.c > 0.5 | frac.g > 0.5 | frac.a > 0.5)

test <- filter_tails %>% 
  mutate(frac.t = str_count(filter_tails$revcomp, "T")/tail_len, 
         frac.c = str_count(filter_tails$revcomp, "C")/tail_len,
         frac.g = str_count(filter_tails$revcomp, "G")/tail_len, 
         frac.a = str_count(filter_tails$revcomp, "A")/tail_len) %>% 
  filter(., revcomp == "" | frac.t > 0.5 | frac.c > 0.5 | frac.g > 0.5 | frac.a > 0.5) %>% 
  group_by(sample, rep, length) %>% 
  summarize(count=n()) 

test1 <- test %>%
  group_by(sample, rep) %>%
  mutate(norm = (count/sum(count))*100)

test <- filter_tails %>% 
  mutate(frac.t = str_count(filter_tails$revcomp, "T")/tail_len, 
         frac.c = str_count(filter_tails$revcomp, "C")/tail_len,
         frac.g = str_count(filter_tails$revcomp, "G")/tail_len, 
         frac.a = str_count(filter_tails$revcomp, "A")/tail_len) %>% 
  filter(., revcomp == "" | frac.t > 0.5 | frac.c > 0.5 | frac.g > 0.5 | frac.a > 0.5) %>% 
  group_by(sample, length) %>% 
  summarize(count=n())

# same plot but overlayed
test1 <- test %>%
  group_by(sample) %>%
  mutate(norm = (count/sum(count))*100)
#write.xlsx(x = test1, file = "/binf-isilon/PBgrp/xpj980/R_scripts_collected/histogram_EV3.xlsx")

#### EV3 #####
densityplot <- ggplot(test1, aes(x=length, weight = norm, color = sample, fill = sample, alpha =0.01)) + 
  geom_histogram(binwidth = 1) + 
  # stat_density(kernel = "gaussian", 
  #              geom = "area", 
  #              adjust = 0.5, 
  #              position = "identity") + 
  cowplot::theme_cowplot() +
  facet_wrap(~sample, ncol =1) +
  #xlim(0, 40) +
  ggtitle("distribution of tail length AFB2 3'CF") +
  xlab("length") +
  #ylab("% of reads") +
  scale_color_manual(values = wes_palette("Moonrise3", n = 5)) +
  scale_fill_manual(values = wes_palette("Moonrise3", n = 5)) 
#ggsave(densityplot, file = "/binf-isilon/PBgrp/xpj980/datascratch/3RACE/05.figures/AFB23CFdistriboverlay_HISTO.svg", units = "cm", width = 30, height = 15, dpi = 300)

### base composition ###
#gather frac columns so that I can facet_grid on base 

frac.table <- tails_filtered %>% 
  gather(., base, fraction, 7:10)

frac.table$sample = factor(frac.table$sample, levels = c("WT", "rdr6", "urt1", "heso1", "heso1urt1"))    
frac.table$base = factor(frac.table$base, levels = c("frac.a", "frac.t", "frac.c", "frac.g"))

## Peter asked to see a distibution of the types of tails 
uendings <- tails_filtered %>% 
  mutate(tailtype = if_else(frac.a ==1, "Atail", "no"), 
         tailtype = if_else(str_detect(revcomp, "T$"), "Uending", tailtype),
         tailtype = if_else(revcomp == "", "tailless", tailtype),
         tailtype = case_when(tailtype=="no" ~"other", 
                              T ~ tailtype)) 
correctNA <- uendings %>% 
  mutate(tailtype = case_when(is.na(uendings$tailtype) ~ "tailless", 
                              T ~ tailtype)) 
#filter(tail_len < 40)

test <- correctNA %>% 
  group_by(sample, tailtype) %>% 
  summarise(count = n())   

plotting <- test %>% 
  group_by(sample) %>% 
  mutate(perc = count/sum(count)*100, 
         sample = factor(sample, c("WT", "urt1", "MIR393")), 
         tailtype = factor(tailtype, c("tailless", "other", "Uending", "Atail"))) %>%
  ggplot(., aes(x = tailtype, y = perc, fill = tailtype)) + 
  geom_bar(stat = "identity", position = "dodge") +
  #geom_point() +
  facet_wrap(~sample, nrow = 5) +
  coord_flip() +
  #theme(aspect.ratio = 1) + 
  cowplot::theme_cowplot() +
  cowplot::background_grid() +
  ggtitle("Tailtype on AFB2 3'CF") +
  ylab("percentage of reads") +
  xlab("Type of tail") +
  scale_fill_manual(values = rev(wes_palette("Royal1", n = 4))) +
  theme(aspect.ratio = 1)
#cowplot::save_plot("/binf-isilon/PBgrp/xpj980/datascratch/3RACE/05.figures/tailtypeAFB23CF_DEDUP_MIR393.svg", plotting, base_asp = 1, nrow = 5)

# and replicates seperated 
test <- correctNA %>% 
  group_by(rep, sample, tailtype) %>% 
  summarise(count = n())   

plotting <- test %>% 
  group_by(rep, sample) %>% 
  mutate(perc = count/sum(count)*100, 
         sample = factor(sample, c("WT", "urt1", "MIR393")), 
         tailtype = factor(tailtype, c("tailless", "other", "Uending", "Atail"))) %>%
  ggplot(., aes(x = tailtype, y = perc, fill = tailtype)) + 
  geom_bar(stat = "identity", position = "dodge") +
  #geom_point() +
  facet_grid(rep~sample) +
  coord_flip() +
  #theme(aspect.ratio = 1) + 
  cowplot::theme_cowplot() +
  cowplot::background_grid() +
  ggtitle("Tailtype on AFB2 3'CF") +
  ylab("percentage of reads") +
  xlab("Type of tail") +
  scale_fill_manual(values = rev(wes_palette("Royal1", n = 4))) +
  theme(aspect.ratio = 1)
#cowplot::save_plot("/binf-isilon/PBgrp/xpj980/datascratch/3RACE/05.figures/tailtypeAFB23CF_reps_DEDUP_MIR393.svg", plotting, base_asp = 1, nrow = 5)

# make errorbars instead 
#### Figure 4B #####
test <- correctNA %>% 
  group_by(sample, rep, tailtype) %>% 
  summarise(count = n())   
#write.xlsx(x = test, file = "/binf-isilon/PBgrp/xpj980/R_scripts_collected/Figure4B.xlsx")

plotting <- test %>% 
  group_by(sample, rep) %>% 
  mutate(perc = count/sum(count)*100, 
         sample = factor(sample, c("WT", "urt1", "MIR393")), 
         tailtype = factor(tailtype, c("tailless", "other", "Uending", "Atail"))) %>%
  ggplot(., aes(x = tailtype, y = perc, fill = tailtype)) + 
  #geom_bar(stat = "identity", position = "dodge") +
  stat_summary(geom = "bar", fun = "mean", position = "dodge") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position = "dodge", size = 0.2) +
  #geom_point() +
  facet_wrap(~sample, nrow = 5) +
  coord_flip() +
  #theme(aspect.ratio = 1) + 
  cowplot::theme_cowplot() +
  cowplot::background_grid() +
  ggtitle("Tailtype on AFB2 3'CF") +
  ylab("percentage of reads") +
  xlab("Type of tail") +
  scale_fill_manual(values = rev(wes_palette("Royal1", n = 4))) +
  theme(aspect.ratio = 1)
cowplot::save_plot("/binf-isilon/PBgrp/xpj980/datascratch/3RACE/05.figures/tailtypeAFB23CF_DEDUP_errorbars_MIR393.svg", plotting, base_asp = 1, nrow = 5)

### count Us in Uending tails 

correctNA$sample = factor(correctNA$sample, levels = c("WT", "urt1", "MIR393"))    

look <- correctNA %>% 
  filter(tailtype == "Uending") %>%
  mutate(onlyU = str_extract(string = revcomp, pattern = regex("T+$")), 
         Uno = str_count(onlyU)) %>% 
  group_by(sample, rep, Uno) %>% 
  summarise(count = n()) %>% 
  group_by(sample, rep) %>% 
  mutate(perc = count/sum(count)*100) %>%
  ggplot(., aes(x = Uno, y = perc, fill = as.factor(Uno))) + 
  #geom_bar(stat = "identity", position = "dodge") + 
  stat_summary(geom = "bar", fun = "mean", position = "dodge") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position = "dodge", size = 0.2) +
  facet_wrap(~sample, nrow = 5) +
  scale_fill_manual(values = wes_palette("IsleofDogs1", n = 6)) +
  cowplot::theme_cowplot() +
  cowplot::background_grid() +
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks = c(1,2,3,4,5,6))
#ggsave(look, file = "/binf-isilon/PBgrp/xpj980/datascratch/3RACE/05.figures/AFB23CFUendingUnumber_DEDUP_errorbars_MIR393.svg", units = "cm", width = 30, height = 15, dpi = 300)

look <- correctNA %>% 
  filter(tailtype == "Uending") %>%
  mutate(onlyU = str_extract(string = revcomp, pattern = regex("T+$")), 
         Uno = str_count(onlyU)) %>% 
  group_by(sample, rep, Uno) %>% 
  add_tally() %>% 
  ggplot(., aes(x = Uno, y = n, fill = as.factor(Uno))) + 
  geom_bar(stat = "identity", position = "dodge") + 
  facet_grid(rep~sample) +
  scale_fill_manual(values = wes_palette("IsleofDogs1", n = 6)) +
  cowplot::theme_cowplot() +
  cowplot::background_grid() +
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks = c(1,2,3,4,5,6))
#ggsave(look, file = "/binf-isilon/PBgrp/xpj980/datascratch/3RACE/05.figures/AFB23CFUendingUnumber_reps_sep.svg", units = "cm", width = 30, height = 15, dpi = 300)

#### AFB2 5CF RRP4###### 

file_names <- list.files("/binf-isilon/PBgrp/xpj980/datascratch/3RACE/RACEseq2/02.merged_lanes/txtfiles_deduplicator/AFB2_NIBBLE/", full.names = T, pattern = "ends")

DF <- tibble()
for (f in file_names) {
  print(f)
  dat <- read_lines(f) %>% 
    as_tibble() %>%
    mutate(name = f %>% str_replace(".*//(.*)_read2.*", "\\1"),
           length=nchar(value))
  DF <- rbind(DF, dat)
}

df <- DF %>%
  mutate(barcode=str_sub(value,1,15),
         delim=str_sub(value,16,20),
         read=str_sub(value,21,-1L))

deduplicated <- df %>%
  group_by(read, name) %>%
  mutate(is_dup=similar_to_any2(barcode, max_diff = 2)) %>%
  arrange(read) %>%
  group_by(read, name, is_dup) %>%
  mutate(n=row_number()) %>% 
  filter(n==1 | (!is_dup))

get_tail <- function(x) {
  tail <- str_remove(x, "AACAATGCGA|AACAATGCG|AACAATGC|AACAATG|AACAAT|AACAA|AACA|AAC|AA|A|")
  trim <- -(10-nchar(x)+nchar(tail))
  return(list(tail = tail, trim = trim))
}

DF2 <- deduplicated %>% 
  ungroup() %>%
  mutate(revcomp = reverseComplement(deduplicated$read, case = "upper"), 
         tail = get_tail(revcomp)$tail, 
         trim = get_tail(revcomp)$trim, 
         tail_len = nchar(tail)) %>% 
  separate(name, c("sample", "rep"), sep="_")

rem_seq_errors <- DF2 %>% 
  mutate(to_discard = case_when(-(trim) == tail_len ~"discard", 
                                T ~ "keep"), 
         to_discard = case_when(trim == 0 & tail_len == 0 ~ "keep", 
                                T ~ to_discard), 
         tail_len = as.double(tail_len)) %>% 
  filter(to_discard != "discard") 

add_pos <- rem_seq_errors %>%
  mutate(position = case_when(trim < 0 & tail_len == 0 ~ trim, 
                              T ~ tail_len))

to_plot <- add_pos %>% 
  group_by(sample, rep, position) %>% 
  filter(tail_len <20) %>%
  summarize(count=n()) %>% 
  ungroup() %>% 
  group_by(sample, rep) %>% 
  mutate(ncount=count/sum(count)) 

to_plot2 <- rem_seq_errors %>% 
  group_by(sample, rep, trim, tail_len) %>% 
  filter(tail_len <10) %>%
  summarize(count=n()) %>% 
  ungroup() %>% 
  group_by(sample, rep) %>% 
  mutate(ncount=count/sum(count)) 

### BLAKE MEYERS PLOT ###
blake_plot <- ggplot(to_plot2, aes(x = trim, y = tail_len)) + 
  geom_point(aes(color = rep, size = ncount), shape = 21) +
  scale_colour_manual(values=c("rep1"="chocolate3","rep2"="cyan4", "rep3" = "firebrick")) + 
  scale_size_area(limits = c(0,1), breaks = c(0,0.001,0.01,0.1,1)) +
  facet_wrap(~sample, nrow = 1) + 
  #facet_grid(sample ~ rep) +
  scale_x_continuous(breaks = c(-10:0)) + 
  scale_y_continuous(breaks = c(0:10)) + 
  cowplot::theme_cowplot() +
  cowplot::background_grid() +
  ggtitle("Nibbling and tailing AFB2 5'CF") +
  xlab("nibbling") +
  ylab("tailing") 

#ggsave(blake_plot, file = "/binf-isilon/PBgrp/xpj980/datascratch/3RACE/05.figures/BlakeMeyersrep.svg", units = "cm", width = 32, height = 9, dpi = 300)

to_plot2 <- rem_seq_errors %>% 
  group_by(sample, trim, tail_len) %>% 
  filter(tail_len <10) %>%
  summarize(count=n()) %>% 
  ungroup() %>% 
  group_by(sample) %>% 
  mutate(ncount=count/sum(count)) 

to_plot2$sample = factor(to_plot2$sample, levels = c("WT", "rrp4"))

blake_plot_avg <- ggplot(to_plot2, aes(x = trim, y = tail_len)) + 
  geom_point(aes(color = sample, size = ncount)) +
  scale_colour_manual(values = rev(wes_palette("Darjeeling2", n = 3))) +
  #scale_colour_manual(values=c("rep1"="chocolate3","rep2"="cyan4", "rep3" = "firebrick")) + 
  scale_size_area(limits = c(0,1), breaks = c(0,0.001,0.01,0.1,1)) +
  facet_wrap(~sample, nrow = 1) +  
  scale_x_continuous(breaks = c(-10:0)) + 
  scale_y_continuous(breaks = c(0:10)) + 
  cowplot::theme_cowplot() +
  cowplot::background_grid() +
  ggtitle("Nibbling and tailing AFB2 5'CF") +
  xlab("nibbling") +
  ylab("tailing") 
#ggsave(blake_plot_avg, file = "/binf-isilon/PBgrp/xpj980/datascratch/3RACE/05.figures/BlakeMeyersavg_DEDUP_RRP4.svg", units = "cm", width = 32, height = 9, dpi = 300)

# keep reps seperated 

to_plot2 <- rem_seq_errors %>% 
  group_by(sample, rep, trim, tail_len) %>% 
  filter(tail_len <10) %>%
  summarize(count=n()) %>% 
  ungroup() %>% 
  group_by(sample, rep) %>% 
  mutate(ncount=count/sum(count)) 

to_plot2$sample = factor(to_plot2$sample, levels = c("WT", "rrp4"))

blake_plot_avg <- ggplot(to_plot2, aes(x = trim, y = tail_len)) + 
  geom_point(aes(color = sample, size = ncount)) +
  scale_colour_manual(values = rev(wes_palette("Darjeeling2", n = 3))) +
  #scale_colour_manual(values=c("rep1"="chocolate3","rep2"="cyan4", "rep3" = "firebrick")) + 
  scale_size_area(limits = c(0,1), breaks = c(0,0.001,0.01,0.1,1)) +
  #facet_wrap(~sample, nrow = 1) +  
  facet_grid(sample ~ rep) +
  scale_x_continuous(breaks = c(-10:0)) + 
  scale_y_continuous(breaks = c(0:10)) + 
  cowplot::theme_cowplot() +
  cowplot::background_grid() +
  ggtitle("Nibbling and tailing AFB2 5'CF") +
  xlab("nibbling") +
  ylab("tailing") 
#ggsave(blake_plot_avg, file = "/binf-isilon/PBgrp/xpj980/datascratch/3RACE/05.figures/BlakeMeyersavg_DEDUP_reps_RRP4.svg", units = "cm", width = 32, height = 9, dpi = 300)

how_many <- rem_seq_errors %>% 
  group_by(sample) %>% 
  summarise(n = n())

how_many_reps <- rem_seq_errors %>% 
  group_by(sample, rep) %>% 
  summarise(n = n())

# only U-tail

Utail <- c("T", "TT", "TTT", "TTTT", "TTTTT", "TTTTTT", "TTTTTTT", "TTTTTTTT", "TTTTTTTTT", "TTTTTTTTTT")

to_plot3 <- rem_seq_errors %>% 
  filter(tail_len < 6) %>% 
  filter(trim > -6) %>%
  filter(tail %in% c("", Utail)) %>%
  group_by(sample, trim, tail_len) %>%
  summarize(count=n()) %>% 
  ungroup() %>% 
  group_by(sample) %>% 
  mutate(ncount=count/sum(count)) 

to_plot3$sample = factor(to_plot3$sample, levels = c("WT", "rrp4"))

blake_plot_onlyU <- ggplot(to_plot3, aes(x = trim, y = tail_len)) + 
  geom_point(aes(color = sample, size = ncount)) +
  scale_colour_manual(values = rev(wes_palette("Darjeeling2", n = 3))) +
  #scale_colour_manual(values=c("rep1"="chocolate3","rep2"="cyan4", "rep3" = "firebrick")) + 
  scale_size_area(limits = c(0,1), breaks = c(0.2,0.4,0.6,0.8,1)) +
  facet_wrap(~sample, nrow = 1) +  
  scale_x_continuous(breaks = c(-5:0)) + 
  scale_y_continuous(breaks = c(0:6)) + 
  cowplot::theme_cowplot() +
  cowplot::background_grid() +
  ggtitle("Nibbling and U-tailing AFB2 5'CF") +
  xlab("nibbling") +
  ylab("tailing") 
#ggsave(blake_plot_onlyU, file = "/binf-isilon/PBgrp/xpj980/datascratch/3RACE/05.figures/060124_BlakeMeyersavgUtail_DEDUP_RRP4.svg", units = "cm", width = 32, height = 16, dpi = 300)

how_many <- rem_seq_errors %>% 
  filter(tail_len < 10) %>% 
  filter(tail %in% c("", Utail)) %>%
  group_by(sample) %>% 
  summarise(n = n()) 

how_many_reps <- rem_seq_errors %>% 
  filter(tail_len < 10) %>% 
  filter(tail %in% c("", Utail)) %>%
  group_by(sample, rep, tail) %>% 
  summarise(n = n()) 

# keep replicates seperate and only look at U tails 
to_plot3 <- rem_seq_errors %>% 
  filter(tail_len < 6) %>% 
  filter(trim > -6) %>%
  filter(tail %in% c("", Utail)) %>%
  group_by(sample, rep, trim, tail_len) %>%
  summarize(count=n()) %>% 
  ungroup() %>% 
  group_by(sample, rep) %>% 
  mutate(ncount=count/sum(count)) 
#write.xlsx(x = to_plot3, file = "/binf-isilon/PBgrp/xpj980/R_scripts_collected/Figure3C.xlsx")

to_plot3$sample = factor(to_plot3$sample, levels = c("WT", "rrp4"))

#### Figure 3C #####
blake_plot_onlyU <- ggplot(to_plot3, aes(x = trim, y = tail_len)) + 
  geom_point(aes(color = sample, size = ncount)) +
  scale_colour_manual(values = rev(wes_palette("Darjeeling2", n = 3))) +
  #scale_colour_manual(values=c("rep1"="chocolate3","rep2"="cyan4", "rep3" = "firebrick")) + 
  scale_size_area(limits = c(0,1), breaks = c(0,0.2,0.4,0.6,0.8,1)) +
  facet_grid(rep~sample) +  
  scale_x_continuous(breaks = c(-5:0)) + 
  scale_y_continuous(breaks = c(0:6)) + 
  cowplot::theme_cowplot() +
  cowplot::background_grid() +
  ggtitle("Nibbling and U-tailing AFB2 5'CF") +
  xlab("nibbling") +
  ylab("tailing") 
#ggsave(blake_plot_onlyU, filename = "/binf-isilon/PBgrp/xpj980/datascratch/3RACE/05.figures/060124_Meyerplot_reps_DEDUP_RRP4.svg", units = "cm", width = 32, height = 32, dpi = 300)

## Maybe I can also filter by saying that the tail should at least have 50% of same base? 

test <- DF2 %>% 
  mutate(frac.t = str_count(DF2$tail, "T")/tail_len, 
         frac.c = str_count(DF2$tail, "C")/tail_len,
         frac.g = str_count(DF2$tail, "G")/tail_len,
         frac.a = str_count(DF2$tail, "A")/tail_len) %>% 
  filter(., tail == "" | frac.t > 0.5 | frac.c > 0.5 | frac.g > 0.5 | frac.a > 0.5)

how_many <- test %>% 
  group_by(sample) %>% 
  summarise(n = n())

how_many_reps <- test %>% 
  group_by(sample, rep) %>% 
  summarise(n = n())

test_to_plot <- test %>% 
  group_by(sample, trim, tail_len) %>% 
  filter(tail_len <10) %>%
  summarize(count=n()) %>% 
  ungroup() %>% 
  group_by(sample) %>% 
  mutate(ncount=count/sum(count)) 

test_to_plot$sample = factor(test_to_plot$sample, levels = c("WT", "rrp4"))

blake_plot_avg_test <- ggplot(test_to_plot, aes(x = trim, y = tail_len)) + 
  geom_point(aes(color = sample, size = ncount)) +
  scale_colour_manual(values = rev(wes_palette("Darjeeling2", n = 3))) +
  #scale_colour_manual(values=c("rep1"="chocolate3","rep2"="cyan4", "rep3" = "firebrick")) + 
  scale_size_area(limits = c(0,1), breaks = c(0,0.001,0.01,0.1,1)) +
  facet_wrap(~sample, nrow = 1) +  
  scale_x_continuous(breaks = c(-10:0)) + 
  scale_y_continuous(breaks = c(0:10)) + 
  cowplot::theme_cowplot() +
  cowplot::background_grid() +
  ggtitle("Nibbling and tailing AFB2 5'CF") +
  xlab("nibbling") +
  ylab("tailing") 
#ggsave(blake_plot_avg_test, file = "/binf-isilon/PBgrp/xpj980/datascratch/3RACE/05.figures/BlakeMeyersavgtailreq_DEDUP_RRP4.svg", units = "cm", width = 32, height = 9, dpi = 300)

### Blake Meyers plot zoomed out ### 
test_to_plot <- test %>% 
  group_by(sample, trim, tail_len) %>% 
  summarize(count=n()) %>% 
  ungroup() %>% 
  group_by(sample) %>% 
  mutate(ncount=count/sum(count)) 

test_to_plot$sample = factor(test_to_plot$sample, levels = c("WT", "rrp4"))

blake_plot_avg_test <- ggplot(test_to_plot, aes(x = trim, y = tail_len)) + 
  geom_point(aes(color = sample, size = ncount)) +
  scale_colour_manual(values = rev(wes_palette("Darjeeling2", n = 3))) +
  #scale_colour_manual(values=c("rep1"="chocolate3","rep2"="cyan4", "rep3" = "firebrick")) + 
  scale_size_area(limits = c(0,1), breaks = c(0,0.001,0.01,0.1,1)) +
  facet_wrap(~sample, nrow = 1) +  
  scale_x_continuous(breaks = c(-10:0)) + 
  cowplot::theme_cowplot() +
  cowplot::background_grid() +
  ggtitle("Nibbling and tailing AFB2 5'CF") +
  xlab("nibbling") +
  ylab("tailing") 
#ggsave(blake_plot_avg_test, file = "/binf-isilon/PBgrp/xpj980/datascratch/3RACE/05.figures/BlakeMeyersavgtailreqnozoom_DEDUP.svg", units = "cm", width = 32, height = 9, dpi = 300)

### TAIL COMPOSITION ###

test2 <- test %>% 
  group_by(sample, tail) %>%
  filter(tail_len < 20) %>% 
  summarise(count = n()) %>% 
  group_by(sample) %>% 
  mutate(ncount = count/sum(count)) %>% 
  filter(tail != "") %>% 
  mutate(frac.u = str_count(tail, "T")/nchar(tail), 
         frac.c = str_count(tail, "C")/nchar(tail),
         frac.g = str_count(tail, "G")/nchar(tail),
         frac.a = str_count(tail, "A")/nchar(tail)) %>% 
  gather(base, basecount, contains("frac"))

test2$sample = factor(test2$sample, levels = c("WT", "rrp4"))    
plot <- ggplot(test2, aes(x = factor(sample), y = (basecount*ncount)*100, fill = base %>% str_remove("frac.") %>% toupper())) +
  geom_bar(stat="identity", width = 0.7) +
  labs(x = "genotype", y = "% of reads with a tail", fill = "base") +
  cowplot::theme_cowplot() +
  cowplot::background_grid() +
  #geom_text(aes(label = paste0(Round_off,"%")), 
  #position = position_stack(vjust = 0.5), size = 4) +
  scale_fill_manual(values = wes_palette("GrandBudapest1", n = 4)) +
  ggtitle("Proportion of reads with a tail and composition of the tails (1-20 nt long)") 
#ggsave(plot, file = "/binf-isilon/PBgrp/xpj980/datascratch/3RACE/05.figures/tailcomposition_RRP4.svg", units = "in", width = 30, height = 30, dpi = 300)

# again, plot replicates individually 
test2 <- test %>% 
  group_by(rep, sample, tail) %>%
  filter(tail_len < 20) %>% 
  summarise(count = n()) %>% 
  group_by(rep, sample) %>% 
  mutate(ncount = count/sum(count)) %>% 
  filter(tail != "") %>% 
  mutate(frac.u = str_count(tail, "T")/nchar(tail), 
         frac.c = str_count(tail, "C")/nchar(tail),
         frac.g = str_count(tail, "G")/nchar(tail),
         frac.a = str_count(tail, "A")/nchar(tail)) %>% 
  gather(base, basecount, contains("frac"))

test2$sample = factor(test2$sample, levels = c("WT", "rrp4"))    
plot <- ggplot(test2, aes(x = factor(sample), y = (basecount*ncount)*100, fill = base %>% str_remove("frac.") %>% toupper())) +
  geom_bar(stat="identity", width = 0.7) +
  labs(x = "genotype", y = "% of reads with a tail", fill = "base") +
  cowplot::theme_cowplot() +
  cowplot::background_grid() +
  facet_wrap(~rep) +
  #geom_text(aes(label = paste0(Round_off,"%")), 
  #position = position_stack(vjust = 0.5), size = 4) +
  scale_fill_manual(values = wes_palette("GrandBudapest1", n = 4)) +
  ggtitle("Proportion of reads with a tail and composition of the tails (1-20 nt long)") 
#ggsave(plot, file = "/binf-isilon/PBgrp/xpj980/datascratch/3RACE/05.figures/test_reps_DEDUP_rrp4.svg", units = "in", width = 30, height = 30, dpi = 300)

### how many AFB2 reads per sample 

read_no <- rem_seq_errors %>%
  filter(tail_len < 10) %>%
  group_by(sample, rep) %>% 
  summarize(count=n()) %>% 
  mutate(target = "AFB2")

plot <-  read_no %>% 
  mutate(class = fct_reorder(sample, count, .fun='median')) %>%
  ggplot(aes(x = reorder(class, count), y=count, fill=rep)) + 
  geom_bar(stat="identity", position = "dodge", width=0.2) +
  scale_fill_manual(values=c("rep1"="chocolate3","rep2"="cyan4", "rep3" = "firebrick")) +
  scale_y_continuous(breaks = c(100,200,300,400,500,600,700,800)) + 
  cowplot::theme_cowplot() +
  cowplot::background_grid() +
  ggtitle("amount of reads of AFB2 5'CF") +
  xlab("genotype") +
  ylab("number of reads") 
#ggsave(plot, file = "/binf-isilon/PBgrp/xpj980/datascratch/3RACE/05.figures/AFB25CFreadsno_DEDUP.svg", units = "in", width = 30, height = 30, dpi = 300)

#### RACESeq set3 #####
#### AFB2 5CF new heso reps###### 
# too many reads to run the functions fast. I have therefore saved the deduplicated DF

#file_names <- list.files("/binf-isilon/PBgrp/xpj980/datascratch/3RACE/RACEseq3/AFB2NIBBLE2/", full.names = T, pattern = "ends")

#DF <- tibble()
#for (f in file_names) {
# print(f)
#dat <- read_lines(f) %>% 
# as_tibble() %>%
#mutate(name = f %>% str_replace(".*//(.*)_read2.*", "\\1"),
#      length=nchar(value))
#  DF <- rbind(DF, dat)
#}

#df <- DF %>%
# mutate(barcode=str_sub(value,1,15),
#       delim=str_sub(value,16,20),
#      read=str_sub(value,21,-1L)) %>% 
#filter(length > 19) %>% 
#filter(length < 45)

# there are so many reads that I will split them up in 12 samples and then merge into one deduplicated table afterwards 

#one <- df %>%
# filter(name == "GKHeso1_R1") %>%
#group_by(read, name) %>%
#mutate(is_dup=similar_to_any2(barcode, max_diff = 2)) %>%
#arrange(read) %>%
#group_by(read, name, is_dup) %>%
#mutate(n=row_number()) %>% 
#filter(n==1 | (!is_dup))

#two <-  df %>%
# filter(name == "GKHeso1_R2") %>%
#group_by(read, name) %>%
#mutate(is_dup=similar_to_any2(barcode, max_diff = 2)) %>%
#arrange(read) %>%
#group_by(read, name, is_dup) %>%
#mutate(n=row_number()) %>% 
#filter(n==1 | (!is_dup))

#three <- df %>%
# filter(name == "GKHeso1_R3") %>%
#group_by(read, name) %>%
#mutate(is_dup=similar_to_any2(barcode, max_diff = 2)) %>%
#arrange(read) %>%
#group_by(read, name, is_dup) %>%
#mutate(n=row_number()) %>% 
#filter(n==1 | (!is_dup))

#four <-  df %>%
# filter(name == "GKHeso1_R4") %>%
#group_by(read, name) %>%
#mutate(is_dup=similar_to_any2(barcode, max_diff = 2)) %>%
#arrange(read) %>%
#group_by(read, name, is_dup) %>%
#mutate(n=row_number()) %>% 
#filter(n==1 | (!is_dup))

#five <-  df %>%
# filter(name == "GKHeso1urt1_R1") %>%
#group_by(read, name) %>%
#mutate(is_dup=similar_to_any2(barcode, max_diff = 2)) %>%
#arrange(read) %>%
#group_by(read, name, is_dup) %>%
#mutate(n=row_number()) %>% 
#filter(n==1 | (!is_dup))

#six <- df %>%
# filter(name == "GKHeso1urt1_R2") %>%
#group_by(read, name) %>%
#mutate(is_dup=similar_to_any2(barcode, max_diff = 2)) %>%
#arrange(read) %>%
#group_by(read, name, is_dup) %>%
#mutate(n=row_number()) %>% 
#filter(n==1 | (!is_dup))

#seven <- df %>%
# filter(name == "GKHeso1urt1_R4") %>%
#group_by(read, name) %>%
#mutate(is_dup=similar_to_any2(barcode, max_diff = 2)) %>%
#arrange(read) %>%
#group_by(read, name, is_dup) %>%
#mutate(n=row_number()) %>% 
#filter(n==1 | (!is_dup))

#eight <- df %>%
# filter(name == "urt1_R2") %>%
#group_by(read, name) %>%
#mutate(is_dup=similar_to_any2(barcode, max_diff = 2)) %>%
#arrange(read) %>%
#group_by(read, name, is_dup) %>%
#mutate(n=row_number()) %>% 
#filter(n==1 | (!is_dup))

#nine <- df %>%
# filter(name == "urt1_R3") %>%
#group_by(read, name) %>%
#mutate(is_dup=similar_to_any2(barcode, max_diff = 2)) %>%
#arrange(read) %>%
#group_by(read, name, is_dup) %>%
#mutate(n=row_number()) %>% 
#filter(n==1 | (!is_dup))

#ten <- df %>%
# filter(name == "WT_R2") %>%
#group_by(read, name) %>%
#mutate(is_dup=similar_to_any2(barcode, max_diff = 2)) %>%
#arrange(read) %>%
#group_by(read, name, is_dup) %>%
#mutate(n=row_number()) %>% 
#filter(n==1 | (!is_dup))

#eleven <- df %>%
# filter(name == "WT_R3") %>%
#group_by(read, name) %>%
#mutate(is_dup=similar_to_any2(barcode, max_diff = 2)) %>%
#arrange(read) %>%
#group_by(read, name, is_dup) %>%
#mutate(n=row_number()) %>% 
#filter(n==1 | (!is_dup))

#twelve <- df %>%
# filter(name == "WT_R4") %>%
#  group_by(read, name) %>%
# mutate(is_dup=similar_to_any2(barcode, max_diff = 2)) %>%
#arrange(read) %>%
#group_by(read, name, is_dup) %>%
#mutate(n=row_number()) %>% 
#filter(n==1 | (!is_dup))

#deduplicated <- bind_rows(one, two, three, four, five, six, seven, eight, nine, ten, eleven, twelve)
#save(deduplicated, file ="/binf-isilon/PBgrp/xpj980/datascratch/3RACE/RACEseq3/deduplicatedDF.Rdat")
load(file ="/binf-isilon/PBgrp/xpj980/datascratch/3RACE/RACEseq3/deduplicatedDF.Rdat")

get_tail <- function(x) {
  tail <- str_remove(x, "AACAATGCGA|AACAATGCG|AACAATGC|AACAATG|AACAAT|AACAA|AACA|AAC|AA|A|")
  trim <- -(10-nchar(x)+nchar(tail))
  return(list(tail = tail, trim = trim))
}

DF2 <- deduplicated %>% 
  ungroup() %>%
  mutate(revcomp = reverseComplement(deduplicated$read, case = "upper"), 
         tail = get_tail(revcomp)$tail, 
         trim = get_tail(revcomp)$trim, 
         tail_len = nchar(tail)) %>% 
  separate(name, c("sample", "rep"), sep="_")

rem_seq_errors <- DF2 %>% 
  mutate(to_discard = case_when(-(trim) == tail_len ~"discard", 
                                T ~ "keep"), 
         to_discard = case_when(trim == 0 & tail_len == 0 ~ "keep", 
                                T ~ to_discard), 
         tail_len = as.double(tail_len)) %>% 
  filter(to_discard != "discard") 

add_pos <- rem_seq_errors %>%
  mutate(position = case_when(trim < 0 & tail_len == 0 ~ trim, 
                              T ~ tail_len))

to_plot <- add_pos %>% 
  group_by(sample, rep, position) %>% 
  filter(tail_len <20) %>%
  summarize(count=n()) %>% 
  ungroup() %>% 
  group_by(sample, rep) %>% 
  mutate(ncount=count/sum(count)) 

to_plot2 <- rem_seq_errors %>% 
  group_by(sample, rep, trim, tail_len) %>% 
  filter(tail_len <10) %>%
  summarize(count=n()) %>% 
  ungroup() %>% 
  group_by(sample, rep) %>% 
  mutate(ncount=count/sum(count)) 

### to merge with first library 
rem_seq_errors$lib <- paste0("lib2")

rem_seq_errorslib1select <- rem_seq_errorslib1 %>%
  select(value, sample, rep, length, barcode, delim, read, is_dup, n, revcomp, tail, trim, tail_len, to_discard, lib)

rem_seq_errorslib2 <- rem_seq_errors

gathered <- rbind(rem_seq_errorslib1select, rem_seq_errorslib2)

### BLAKE MEYERS PLOT ###
blake_plot <- ggplot(to_plot2, aes(x = trim, y = tail_len)) + 
  geom_point(aes(color = rep, size = ncount), shape = 21) +
  scale_colour_manual(values=c("rep1"="chocolate3","rep2"="cyan4", "rep3" = "firebrick")) + 
  scale_size_area(limits = c(0,1), breaks = c(0,0.001,0.01,0.1,1)) +
  facet_wrap(~sample, nrow = 1) + 
  #facet_grid(sample ~ rep) +
  scale_x_continuous(breaks = c(-10:0)) + 
  scale_y_continuous(breaks = c(0:10)) + 
  cowplot::theme_cowplot() +
  cowplot::background_grid() +
  ggtitle("Nibbling and tailing AFB2 5'CF") +
  xlab("nibbling") +
  ylab("tailing") 

#ggsave(blake_plot, file = "/binf-isilon/PBgrp/xpj980/datascratch/3RACE/RACEseq3/figures/BlakeMeyersrep.svg", units = "cm", width = 32, height = 9, dpi = 300)

to_plot2 <- rem_seq_errors %>% 
  group_by(sample, trim, tail_len) %>% 
  filter(tail_len <10) %>%
  summarize(count=n()) %>% 
  ungroup() %>% 
  group_by(sample) %>% 
  mutate(ncount=count/sum(count)) 

to_plot2$sample = factor(to_plot2$sample, levels = c("WT", "urt1", "GKHeso1", "GKHeso1urt1"))

blake_plot_avg <- ggplot(to_plot2, aes(x = trim, y = tail_len)) + 
  geom_point(aes(color = sample, size = ncount)) +
  scale_colour_manual(values = rev(wes_palette("Darjeeling2", n = 5))) +
  #scale_colour_manual(values=c("rep1"="chocolate3","rep2"="cyan4", "rep3" = "firebrick")) + 
  scale_size_area(limits = c(0,1), breaks = c(0,0.01,0.1,0.4,1)) +
  facet_wrap(~sample, nrow = 1) +  
  scale_x_continuous(breaks = c(-10:0)) + 
  scale_y_continuous(breaks = c(0:10)) + 
  cowplot::theme_cowplot() +
  cowplot::background_grid() +
  ggtitle("Nibbling and tailing AFB2 5'CF") +
  xlab("nibbling") +
  ylab("tailing") 
#ggsave(blake_plot_avg, file = "/binf-isilon/PBgrp/xpj980/datascratch/3RACE/RACEseq3/figures/BlakeMeyersavg_DEDUP.svg", units = "cm", width = 32, height = 9, dpi = 300)

# keep reps seperated 

to_plot2 <- rem_seq_errors %>% 
  group_by(sample, rep, trim, tail_len) %>% 
  filter(tail_len <10) %>%
  summarize(count=n()) %>% 
  ungroup() %>% 
  group_by(sample, rep) %>% 
  mutate(ncount=count/sum(count)) 

to_plot2$sample = factor(to_plot2$sample, levels = c("WT", "urt1", "GKHeso1", "GKHeso1urt1"))

blake_plot_avg <- ggplot(to_plot2, aes(x = trim, y = tail_len)) + 
  geom_point(aes(color = sample, size = ncount)) +
  scale_colour_manual(values = rev(wes_palette("Darjeeling2", n = 5))) +
  #scale_colour_manual(values=c("rep1"="chocolate3","rep2"="cyan4", "rep3" = "firebrick")) + 
  scale_size_area(limits = c(0,1), breaks = c(0,0.01,0.1,0.4,1)) +
  #facet_wrap(~sample, nrow = 1) +  
  facet_grid(sample ~ rep) +
  scale_x_continuous(breaks = c(-10:0)) + 
  scale_y_continuous(breaks = c(0:10)) + 
  cowplot::theme_cowplot() +
  cowplot::background_grid() +
  ggtitle("Nibbling and tailing AFB2 5'CF") +
  xlab("nibbling") +
  ylab("tailing") 
#ggsave(blake_plot_avg, file = "/binf-isilon/PBgrp/xpj980/datascratch/3RACE/RACEseq3/figures/BlakeMeyersavg_DEDUP_reps.svg", units = "cm", width = 32, height = 9, dpi = 300)

how_many <- rem_seq_errors %>% 
  group_by(sample) %>% 
  summarise(n = n())

how_many_reps <- rem_seq_errors %>% 
  group_by(sample, rep) %>% 
  summarise(n = n())

## only U-tailing 

Utail <- c("T", "TT", "TTT", "TTTT", "TTTTT", "TTTTTT", "TTTTTTT", "TTTTTTTT", "TTTTTTTTT", "TTTTTTTTTT")

to_plot3 <- rem_seq_errors %>% 
  filter(tail_len < 10) %>% 
  filter(tail %in% c("", Utail)) %>%
  group_by(sample, trim, tail_len) %>%
  summarize(count=n()) %>% 
  ungroup() %>% 
  group_by(sample) %>% 
  mutate(ncount=count/sum(count)) 

to_plot3$sample = factor(to_plot3$sample, levels = c("WT", "urt1", "GKHeso1", "GKHeso1urt1"))

blake_plot_onlyU <- ggplot(to_plot3, aes(x = trim, y = tail_len)) + 
  geom_point(aes(color = sample, size = ncount)) +
  scale_colour_manual(values = rev(wes_palette("Darjeeling2", n = 5))) +
  #scale_colour_manual(values=c("rep1"="chocolate3","rep2"="cyan4", "rep3" = "firebrick")) + 
  scale_size_area(limits = c(0,1), breaks = c(0,0.2,0.4,0.6,0.8,1)) +
  facet_wrap(~sample, nrow = 1) +  
  scale_x_continuous(breaks = c(-10:0)) + 
  scale_y_continuous(breaks = c(0:10)) + 
  cowplot::theme_cowplot() +
  cowplot::background_grid() +
  ggtitle("Nibbling and U-tailing AFB2 5'CF") +
  xlab("nibbling") +
  ylab("tailing") 
#ggsave(blake_plot_onlyU, file = "/binf-isilon/PBgrp/xpj980/datascratch/3RACE/RACEseq3/figures/BlakeMeyersavgUtail_DEDUP.svg", units = "cm", width = 32, height = 9, dpi = 300)

how_many <- rem_seq_errors %>% 
  filter(tail_len < 10) %>% 
  filter(tail %in% c("", Utail)) %>%
  group_by(sample) %>% 
  summarise(n = n()) 

how_many_reps <- rem_seq_errors %>% 
  filter(tail_len < 10) %>% 
  filter(tail %in% c("", Utail)) %>%
  group_by(sample, rep) %>% 
  summarise(n = n()) 

# keep replicates seperate and only look at U tails 
to_plot3 <- rem_seq_errors %>% 
  filter(tail_len < 10) %>% 
  filter(tail %in% c("", Utail)) %>%
  group_by(sample, rep, trim, tail_len) %>%
  summarize(count=n()) %>% 
  ungroup() %>% 
  group_by(sample, rep) %>% 
  mutate(ncount=count/sum(count)) 

to_plot3$sample = factor(to_plot3$sample, levels = c("WT", "urt1", "GKHeso1", "GKHeso1urt1"))

blake_plot_onlyU <- ggplot(to_plot3, aes(x = trim, y = tail_len)) + 
  geom_point(aes(color = sample, size = ncount)) +
  scale_colour_manual(values = rev(wes_palette("Darjeeling2", n = 5))) +
  #scale_colour_manual(values=c("rep1"="chocolate3","rep2"="cyan4", "rep3" = "firebrick")) + 
  scale_size_area(limits = c(0,1), breaks = c(0,0.2,0.4,0.6,0.8,1)) +
  facet_grid(rep~sample) +  
  scale_x_continuous(breaks = c(-10:0)) + 
  scale_y_continuous(breaks = c(0:10)) + 
  cowplot::theme_cowplot() +
  cowplot::background_grid() +
  ggtitle("Nibbling and U-tailing AFB2 5'CF") +
  xlab("nibbling") +
  ylab("tailing") 
#ggsave(blake_plot_onlyU, filename = "/binf-isilon/PBgrp/xpj980/datascratch/3RACE/RACEseq3/figures/Meyerplot_reps_DEDUP.svg", units = "cm", width = 32, height = 9, dpi = 300)

#### 5CF GATHERED ####
### make same plot but with reps from both libs 

gathered_toplot <- gathered %>% 
  filter(sample != "rdr6") %>%
  mutate(sample_renamed = case_when(sample == "heso1" ~ "GKHeso1", 
                                    T ~ sample)) 

gathered_toplot2 <- gathered_toplot %>% 
  mutate(sample_renamed2 = case_when(sample_renamed == "heso1urt1" ~ "GKHeso1urt1", 
                                     T ~ sample_renamed)) 
#write.xlsx(x = gathered_toplot, file = "/binf-isilon/PBgrp/xpj980/R_scripts_collected/Figure3B.xlsx")

to_plot3 <- gathered_toplot2 %>% 
  filter(tail_len < 10) %>% 
  filter(trim > -6) %>%
  filter(tail %in% c("", Utail)) %>%
  group_by(sample_renamed2, rep, trim, tail_len, lib) %>%
  summarize(count=n()) %>% 
  ungroup() %>% 
  group_by(sample_renamed2, rep, lib) %>% 
  mutate(ncount=count/sum(count)) 

to_plot3$sample_renamed2 = factor(to_plot3$sample_renamed2, levels = c("WT", "urt1", "GKHeso1", "GKHeso1urt1"))

##### Figure 3B ######
blake_plot_onlyU <- ggplot(to_plot3, aes(x = trim, y = tail_len)) + 
  geom_point(aes(color = sample_renamed2, size = ncount)) +
  scale_colour_manual(values = rev(wes_palette("Darjeeling2", n = 5))) +
  #scale_colour_manual(values=c("rep1"="chocolate3","rep2"="cyan4", "rep3" = "firebrick")) + 
  scale_size_area(limits = c(0,1), breaks = c(0,0.2,0.4,0.6,0.8,1)) +
  facet_grid(rep~sample_renamed2) +  
  scale_x_continuous(breaks = c(-5:0)) + 
  scale_y_continuous(breaks = c(0:5)) + 
  cowplot::theme_cowplot() +
  cowplot::background_grid() +
  ggtitle("Nibbling and U-tailing AFB2 5'CF") +
  xlab("nibbling") +
  ylab("tailing") 
#ggsave(blake_plot_onlyU, filename = "/binf-isilon/PBgrp/xpj980/datascratch/3RACE/RACEseq3/figures/Meyerplot_reps_DEDUP_BOTHLIBS_onlyU_scaled.svg", dpi = 300, units = "cm", width = 32, height = 32)

## a tail should at least have 50% of same base

test <- gathered_toplot2 %>% 
  mutate(frac.t = str_count(gathered_toplot2$tail, "T")/tail_len, 
         frac.c = str_count(gathered_toplot2$tail, "C")/tail_len,
         frac.g = str_count(gathered_toplot2$tail, "G")/tail_len,
         frac.a = str_count(gathered_toplot2$tail, "A")/tail_len) %>% 
  filter(., tail == "" | frac.t > 0.5 | frac.c > 0.5 | frac.g > 0.5 | frac.a > 0.5)

how_many <- test %>% 
  group_by(sample_renamed2) %>% 
  summarise(n = n())

how_many_reps <- test %>% 
  group_by(sample_renamed2, rep) %>% 
  summarise(n = n())

test_to_plot <- test %>% 
  group_by(sample_renamed2, trim, tail_len) %>% 
  filter(tail_len <10) %>%
  summarize(count=n()) %>% 
  ungroup() %>% 
  group_by(sample_renamed2) %>% 
  mutate(ncount=count/sum(count)) 

test_to_plot$sample_renamed2 = factor(test_to_plot$sample_renamed2, levels = c("WT", "urt1", "GKHeso1", "GKHeso1urt1"))

blake_plot_avg_test <- ggplot(test_to_plot, aes(x = trim, y = tail_len)) + 
  geom_point(aes(color = sample_renamed2, size = ncount)) +
  scale_colour_manual(values = rev(wes_palette("Darjeeling2", n = 5))) +
  #scale_colour_manual(values=c("rep1"="chocolate3","rep2"="cyan4", "rep3" = "firebrick")) + 
  scale_size_area(limits = c(0,1), breaks = c(0,0.2,0.4,0.6,0.8,1)) +
  facet_wrap(~sample_renamed2, nrow = 1) +  
  scale_x_continuous(breaks = c(-10:0)) + 
  scale_y_continuous(breaks = c(0:10)) + 
  cowplot::theme_cowplot() +
  cowplot::background_grid() +
  ggtitle("Nibbling and tailing AFB2 5'CF") +
  xlab("nibbling") +
  ylab("tailing") 
#ggsave(blake_plot_avg_test, file = "/binf-isilon/PBgrp/xpj980/datascratch/3RACE/RACEseq3/figures/BlakeMeyersavgtailreq_DEDUP_BOTHLIBS.svg", units = "cm", width = 32, height = 9, dpi = 300)

### Blake Meyers plot zoomed out ### 
test_to_plot <- test %>% 
  group_by(sample_renamed2, trim, tail_len) %>% 
  summarize(count=n()) %>% 
  ungroup() %>% 
  group_by(sample_renamed2) %>% 
  mutate(ncount=count/sum(count)) 

test_to_plot$sample_renamed2 = factor(test_to_plot$sample_renamed2, levels = c("WT", "urt1", "GKHeso1", "GKHeso1urt1"))

blake_plot_avg_test <- ggplot(test_to_plot, aes(x = trim, y = tail_len)) + 
  geom_point(aes(color = sample_renamed2, size = ncount)) +
  scale_colour_manual(values = rev(wes_palette("Darjeeling2", n = 5))) +
  #scale_colour_manual(values=c("rep1"="chocolate3","rep2"="cyan4", "rep3" = "firebrick")) + 
  scale_size_area(limits = c(0,1), breaks = c(0,0.2,0.4,0.6,0.8,1)) +
  facet_wrap(~sample_renamed2, nrow = 1) +  
  scale_x_continuous(breaks = c(-10:0)) + 
  cowplot::theme_cowplot() +
  cowplot::background_grid() +
  ggtitle("Nibbling and tailing AFB2 5'CF") +
  xlab("nibbling") +
  ylab("tailing") 
#ggsave(blake_plot_avg_test, file = "/binf-isilon/PBgrp/xpj980/datascratch/3RACE/05.figures/BlakeMeyersavgtailreqnozoom_DEDUP_BOTHLIBS.svg", units = "cm", width = 32, height = 9, dpi = 300)

### TAIL COMPOSITION ###

test2 <- test %>% 
  group_by(sample_renamed2, tail) %>%
  filter(tail_len < 20) %>% 
  summarise(count = n()) %>% 
  group_by(sample_renamed2) %>% 
  mutate(ncount = count/sum(count)) %>% 
  filter(tail != "") %>% 
  mutate(frac.u = str_count(tail, "T")/nchar(tail), 
         frac.c = str_count(tail, "C")/nchar(tail),
         frac.g = str_count(tail, "G")/nchar(tail),
         frac.a = str_count(tail, "A")/nchar(tail)) %>% 
  gather(base, basecount, contains("frac"))

test2$sample_renamed2 = factor(test2$sample_renamed2, levels = c("WT", "urt1", "GKHeso1", "GKHeso1urt1"))    
plot <- ggplot(test2, aes(x = factor(sample_renamed2), y = (basecount*ncount)*100, fill = base %>% str_remove("frac.") %>% toupper())) +
  geom_bar(stat="identity", width = 0.7) +
  labs(x = "genotype", y = "% of reads with a tail", fill = "base") +
  cowplot::theme_cowplot() +
  cowplot::background_grid() +
  #geom_text(aes(label = paste0(Round_off,"%")), 
  #position = position_stack(vjust = 0.5), size = 4) +
  scale_fill_manual(values = wes_palette("GrandBudapest1", n = 4)) +
  ggtitle("Proportion of reads with a tail and composition of the tails (1-20 nt long)") 
#ggsave(plot, file = "/binf-isilon/PBgrp/xpj980/datascratch/3RACE/RACEseq3/figures/readprop_DEDUP_BOTHLIBS.svg", units = "in", width = 30, height = 30, dpi = 300)

# again, plot replicates individually 
test2 <- test %>% 
  group_by(rep, sample_renamed2, tail) %>%
  filter(tail_len < 20) %>% 
  summarise(count = n()) %>% 
  group_by(rep, sample_renamed2) %>% 
  mutate(ncount = count/sum(count)) %>% 
  filter(tail != "") %>% 
  mutate(frac.u = str_count(tail, "T")/nchar(tail), 
         frac.c = str_count(tail, "C")/nchar(tail),
         frac.g = str_count(tail, "G")/nchar(tail),
         frac.a = str_count(tail, "A")/nchar(tail)) %>% 
  gather(base, basecount, contains("frac"))

test2$sample_renamed2 = factor(test2$sample_renamed2, levels = c("WT", "urt1", "GKHeso1", "GKHeso1urt1"))    
plot <- ggplot(test2, aes(x = factor(sample_renamed2), y = (basecount*ncount)*100, fill = base %>% str_remove("frac.") %>% toupper())) +
  geom_bar(stat="identity", width = 0.7) +
  labs(x = "genotype", y = "% of reads with a tail", fill = "base") +
  cowplot::theme_cowplot() +
  cowplot::background_grid() +
  facet_wrap(~rep) +
  #geom_text(aes(label = paste0(Round_off,"%")), 
  #position = position_stack(vjust = 0.5), size = 4) +
  scale_fill_manual(values = wes_palette("GrandBudapest1", n = 4)) +
  ggtitle("Proportion of reads with a tail and composition of the tails (1-20 nt long)") 
#ggsave(plot, file = "/binf-isilon/PBgrp/xpj980/datascratch/3RACE/RACEseq3/figures/test_reps_DEDUP_BOTHLIBS.svg", units = "in", width = 30, height = 30, dpi = 300)

#### Figure EV2 #####

plot <- ggplot(test2, aes(x = rep, y = (basecount*ncount)*100, fill = base %>% str_remove("frac.") %>% toupper())) +
  geom_bar(stat="identity", width = 0.7) +
  labs(x = "genotype", y = "% of reads with a tail", fill = "base") +
  cowplot::theme_cowplot() +
  cowplot::background_grid() +
  facet_wrap(~sample_renamed2) +
  #geom_text(aes(label = paste0(Round_off,"%")), 
  #position = position_stack(vjust = 0.5), size = 4) +
  scale_fill_manual(values = wes_palette("GrandBudapest1", n = 4)) +
  ggtitle("Proportion of reads with a tail and composition of the tails (1-20 nt long)") 
#ggsave(plot, file = "/binf-isilon/PBgrp/xpj980/datascratch/3RACE/RACEseq3/figures/test_reps_DEDUP_BOTHLIBS_FACETBEST.svg", units = "in", width = 30, height = 30, dpi = 300)
