#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

bracken_file <- 'data/sample_0_bracken_result.txt'
output_file <- 'sample_0_bracken.csv'

# args <- commandArgs(trailingOnly = TRUE)
# 
# if (length(args) != 2) {
#   stop("Usage: Rscript convert_bracken.R <input_file.tsv> <output_file.csv>")
# }
# 
# bracken_file <- args[1]
# output_file <- args[2]

bracken <- read.delim(bracken_file, header = TRUE, stringsAsFactors = FALSE) %>%
  rename(
    tax_id = taxonomy_id,
    abundance = fraction_total_reads
  ) %>%
  select(tax_id, abundance) %>% 
  mutate(abundance = abundance / sum(abundance) * 100) %>% 
  filter(abundance > 0.000001)

write_csv2(bracken, output_file)

cat("Written to:", output_file, "\n")
