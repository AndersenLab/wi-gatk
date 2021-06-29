library(tidyverse)


df_use <- tseries::read.matrix(gzfile("impute_gts.tsv.gz"))
storage.mode(df_use) <- "logical"


n_samples <- ncol(df_use)
results <- sapply(1:10, function(x) {
  cumsum(
    (tidyr::complete(
      (plyr::count(
        apply(
          df_use[,sample(n_samples, n_samples)],
          1,
          function(x) {
            min(which(x))
          }
        )
      )
      ),
      x = 1:n_samples,
      fill = list(freq = as.integer(0))
    ))$freq
  )
})

colnames(results) <- 1:ncol(results)
results <- tbl_df(results) %>% 
  tidyr::gather(rep, sites) %>%
  dplyr::group_by(rep) %>%
  dplyr::mutate(isotypes = dplyr::row_number(rep)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(isotypes)

summarized_results <- results %>% 
  dplyr::group_by(isotypes) %>%
  dplyr::summarize(mean_sites = mean(sites))


summarized_results %>% readr::write_tsv("variant_accumulation.tsv")


ggplot(results) +
  geom_step(aes(y = sites, x = isotypes, color=rep), size = 0.5, alpha=0.5) +
  geom_step(aes(y = mean_sites, x = isotypes), size=2, data=summarized_results) +
  scale_y_continuous(sec.axis = sec_axis(~./nrow(df_use)*100, name="SNVs (%)")) +
  labs(x = "Isotypes", y = "SNVs") +
  theme_bw() +
  theme(legend.position = 'none',
        plot.margin = unit(c(1.0,0.5,0.5,1.0),"cm"),
        axis.text = element_text(size=12, face='bold'))
  
ggsave("variant_accumulation.pdf", width=10, height=10)
