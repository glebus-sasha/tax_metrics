#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

# truth_file <- 'data/short_reads_mapping.tsv'
# output_file <- 'short_reads_truth.csv'

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript convert_truth.R <input_file.tsv> <output_file.csv>")
}

truth_file <- args[1]
output_file <- args[2]

truth_abundance <- read.table(truth_file, comment.char = "@", header = TRUE) %>% 
  rename(anonymous_read_id = X.anonymous_read_id) %>%
  group_by(tax_id) %>%
  summarise(reads = n(), .groups = "drop") %>%
  mutate(abundance = reads / sum(reads)) %>%
  select(tax_id, abundance) %>% 
  mutate(abundance = abundance / sum(abundance) * 100)

write_csv2(truth_abundance, output_file)

cat("Written to:", output_file, "\n")
