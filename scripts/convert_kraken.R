#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

# kraken2_result <- 'data/short_kraken2_result.txt'
# output_file <- 'short_reads_kraken.csv'

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript convert_kraken.R <input_file.tsv> <output_name.csv>")
}

kraken2_result <- args[1]
output_file <- args[2]

kraken_abundance <- read_delim(kraken2_result, 
                               delim = "\t", escape_double = FALSE, col_names = FALSE,
                               trim_ws = TRUE) %>%
  set_names(c("classification", "read_id", "tax_id", "kmer_hits", "kmer_total")) %>% 
  extract(tax_id, into = c("name", "tax_id"), regex = "^(.*) \\(taxid ([0-9]+)\\)$") %>%
  mutate(tax_id = as.integer(tax_id)) %>%
  group_by(tax_id) %>%
  summarise(reads = n(), .groups = "drop") %>%
  mutate(abundance = reads / sum(reads)) %>% 
  select(-reads) %>% 
  mutate(abundance = abundance / sum(abundance) * 100) %>% 
  filter(abundance > 0.000001)

write_csv2(kraken_abundance, output_file)

cat("Written to:", output_file, "\n")
