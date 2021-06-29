library(tidyverse)
try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))
df <- readr::read_tsv("tajima.bed.gz", col_names = c("CHROM", "START", "END","N_SITES", "N_SNPs", "TajimaD"))

df <- df %>% dplyr::mutate(POS = ((END - START)/2) + START)

genetic_scale <- function(n) {
  paste0(n/1E6, "Mb")
}

ggplot(df %>% dplyr::filter(CHROM != "MtDNA")) +
  geom_point(aes(y = TajimaD, x = POS)) +
  facet_grid(.~CHROM, scales="free_x") +
  theme_bw() +
  scale_x_continuous(labels = genetic_scale)

  
ggsave("tajima_d.png", width=15, height = 10)
