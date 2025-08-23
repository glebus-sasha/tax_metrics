#!/usr/bin/env Rscript
 
suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
})

# bracken_file <- 'data/sample_0_bracken_result.txt'
# output_file <- 'sample_0_alpha.csv'

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript alpha_div.R <input_file.csv> <output_file.csv>")
}

bracken_file <- args[1]
output_file <- args[2]

counts <- read_table(bracken_file) %>%
  pull(new_est_reads) %>% 
  as.numeric() %>%
  discard(is.na)

# Считаем все нужные метрики альфа-разнообразия
alpha_df <- tibble(
  observed = specnumber(counts),
  shannon  = diversity(counts, index = "shannon"),
  simpson  = diversity(counts, index = "simpson"),
  invsimpson = diversity(counts, index = "invsimpson")
)

# Сохраняем
write_csv2(alpha_df, output_file)

cat("Written to:", output_file, "\n")
print(alpha_df)
