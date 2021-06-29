#!/usr/bin/env Rscript --vanilla
library(tidyverse)
# Used in calculating isotypes
stack_list <- function(x) {
  if (is.null(names(x)) == T) {
    names(x) <- as.character(1:length(x))
  }
  stack(x)
}

coverage_20 <- (readr::read_tsv("SM_coverage.tsv") %>%
                  dplyr::filter(coverage > 20))$strain


WI_strain_info <- readr::read_tsv("https://docs.google.com/spreadsheets/d/1V6YHzblaDph01sFDI8YK_fP0H7sVebHQTXypGdiQIjI/pub?output=tsv") %>%
  dplyr::select(strain, isotype) %>%
  dplyr::filter(!is.na(isotype)) %>%
  dplyr::filter(strain %in% coverage_20)

# Get number of Sites in VCF

nsites = system("cat filtered.stats.txt  | grep 'SN.*number of SNPs' | cut -f 4", intern = T)
nsites = as.numeric(nsites)


gtcheck <- readr::read_tsv("gtcheck.tsv") %>%
  dplyr::filter((i %in% coverage_20) & (j %in% coverage_20)) %>%
  dplyr::mutate(concordance = (nsites - discordance) / nsites) %>%
  dplyr::mutate(isotype = concordance > 0.9993) 
  
#
# Handle single strains
#

single_strains <- gtcheck %>%
  dplyr::select(i, j, isotype) %>%
  tidyr::gather(col, strain, -isotype) %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(single_strain = (sum(isotype) == 0)) %>%
  dplyr::distinct(.keep_all = T) %>%
  dplyr::select(strain, single_strain) %>% dplyr::filter(single_strain)
# Add LSJ1
single_strains <- dplyr::bind_rows(list(strain = "LSJ1", single_strain = T), single_strains)

single_strains <- as.data.frame(list(i = single_strains$strain,
                                     j = single_strains$strain,
                                     discordance = 0,
                                     concordance = 1,
                                     isotype = TRUE))


gtcheck <- dplyr::bind_rows(gtcheck, single_strains)

iso_gtcheck <- gtcheck %>% dplyr::filter(isotype == T)

strain_list <- sort(unique(c(gtcheck$i, gtcheck$j)))

isotype_groups <- stack_list(unique(lapply(strain_list, function(x) {
  grouped_strains <- dplyr::filter(iso_gtcheck, (i == x | j == x)) %>%
    dplyr::select(i, j)
  sort(unique(unlist(grouped_strains)))
}))) %>%
  dplyr::rename(strain = values, group = ind) %>%
  dplyr::mutate(group = ifelse(strain == "LSJ1", 0, group)) %>%
  dplyr::group_by(group) %>%
  tidyr::nest(strain) %>%
  dplyr::distinct(data, .keep_all = T) %>%
  tidyr::unnest()

isotype_groups <- dplyr::left_join(isotype_groups, WI_strain_info) %>%
  dplyr::left_join(readr::read_tsv("SM_coverage.tsv")) %>%
  dplyr::mutate(group = as.integer(group)) %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(unique_isotypes_per_group = length(unique(purrr::discard(isotype, is.na )))) %>%
  dplyr::ungroup(isotype_groups) %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(strain_conflicts = length(strain) != 1) 

ggplot(gtcheck) +
  geom_histogram(aes(x=concordance, fill = isotype), binwidth = 0.00025) +
  scale_fill_manual(values = c("#808080", "#0080FF"))

ggsave("concordance.pdf")
ggsave("concordance.png")

ggplot(gtcheck) +
  geom_histogram(aes(x=concordance, fill = isotype), binwidth = 0.000025) +
  scale_fill_manual(values = c("#808080", "#0080FF")) +
  scale_x_continuous(limits = c(0.99, 1.0)) +
  labs(x = "Concordance", y = "Number of Comparisons") +
  geom_vline(aes(xintercept = 0.9993), color = "red") +
  theme(axis.title = ggplot2::element_text(size=14, face="bold", color="black", vjust=5))

ggsave("xconcordance.pdf")
ggsave("xconcordance.png")

# Save text files
readr::write_tsv(gtcheck, "gtcheck.tsv")
readr::write_tsv(isotype_groups, "isotype_groups.tsv")

