# This script takes parsed variant summary json files and converts them to TSVs
library(tidyverse)

flist <- list.files(pattern="*summary.json")

for(f in flist) {
    df <- jsonlite::read_json(f,  simplifyDataFrame = TRUE)
    df <- df[,names(df) != ""]
    df %>% 
    replace(is.na(.), 0) %>%
    dplyr::tbl_df()
    df %>%
    readr::write_tsv(gsub("json", "tsv", f))
}
