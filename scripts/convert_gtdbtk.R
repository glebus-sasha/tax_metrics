#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

# gtdbtk_file <- 'data/anonymous_reads.summary.tsv'
# output_file <- 'short_reads_gtdbtk.csv'
# ar122_file <- 'data/ar122_metadata_r202.tsv'
# bac120_file <- 'data/bac120_metadata_r202.tsv'

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop("Usage: Rscript convert_gtdbtk.R <summary.tsv> <ar122_metadata.tsv> <bac120_metadata.tsv> <output.csv>")
}

gtdbtk_file <- args[1]
ar122_file <- args[2]
bac120_file <- args[3]
output_file <- args[4]

# 1. Загрузка и объединение метаданных
ar122 <- read_tsv(ar122_file, show_col_types = FALSE) %>%
  select(gtdb_taxonomy, tax_id = ncbi_taxid)

bac120 <- read_tsv(bac120_file, show_col_types = FALSE) %>%
  select(gtdb_taxonomy, tax_id = ncbi_taxid)

metadata <- bind_rows(ar122, bac120) %>%
  distinct(gtdb_taxonomy, .keep_all = TRUE)

# 2. Загрузка summary.tsv
summary_df <- read.delim(gtdbtk_file, sep = "\t", stringsAsFactors = FALSE) %>% 
  filter(classification != 'classification')

# 3. Join + формирование универсального формата
result <- summary_df %>%
  left_join(metadata, by = c("classification" = "gtdb_taxonomy")) %>%
  mutate(
    tax_id = as.integer(tax_id),
    abundance = 1 / n()
  ) %>%
  select(tax_id, abundance) %>% 
  mutate(abundance = abundance / sum(abundance) * 100)

# 4. Сохранение
write_csv2(result, output_file)
cat("Written to:", output_file, "\n")