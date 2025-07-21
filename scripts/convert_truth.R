#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

# truth_file <- 'data/taxonomic_profile_0.txt'
# output_file <- 'sample_0_truth.csv'

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript convert_truth.R <input_file.tsv> <output_file.csv>")
}

truth_file <- args[1]
output_file <- args[2]

truth_abundance <- readLines(truth_file) %>% 
  discard(~ str_starts(.x, "@") & !str_starts(.x, "@@")) %>%  # сохраняем строку с заголовками
  modify_at(1, ~ str_remove(.x, "^@@")) %>%                   # убираем @@ в первой строке
  {\(x) read_tsv(I(x), na = c("", "NA")) }() %>%              # читаем как TSV
  mutate(
    TAXID = as.integer(TAXID),
    PERCENTAGE = as.numeric(PERCENTAGE)
  ) %>% 
  select(tax_id = TAXID, abundance = PERCENTAGE) %>% 
  mutate(abundance = abundance / sum(abundance) * 100) %>% 
  filter(abundance > 0.000001) %>% 
  arrange(-abundance)

write_csv2(truth_abundance, output_file)

cat("Written to:", output_file, "\n")
