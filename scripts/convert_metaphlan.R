#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

# metaphlan_file <- 'data/metaphlan.txt'
# output_file <- 'short_reads_metaphlan.csv'

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript convert_metaphlan.R <input_file.tsv> <output_file.csv>")
}

metaphlan_file <- args[1]
output_file <- args[2]

metaphlan <- read.delim(metaphlan_file, 
                        comment.char = "#", 
                        header = TRUE, 
                        stringsAsFactors = FALSE,   
                        col.names = c("clade_name", "tax_id", "abundance", "additional_species")) %>%
  mutate(
    # Извлечь последний tax_id, разделённый "|"
    tax_id = str_extract(tax_id, "[^|]+$"),
    tax_id = as.integer(tax_id)
  ) %>%
  select(tax_id, abundance) %>% 
  mutate(abundance = abundance / sum(abundance) * 100) %>% 
  filter(abundance > 0.000001)

write_csv2(metaphlan, output_file)

cat("Written to:", output_file, "\n")
