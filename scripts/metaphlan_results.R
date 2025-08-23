#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

# metaphlan_file <- 'data/sample_0.txt'
# output_file <- 'taxonomy.csv'
 
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript metaphlan_results.R <metaphlan_file.tsv> <output_file.csv>")
}

metaphlan_file <- args[1]
output_file <- args[2]


metaphlan <- read_tsv(
  metaphlan_file,
  comment = "#",
  col_names = c("clade_name", "NCBI_tax_id", "relative_abundance", "additional_species"),
  col_types = cols(.default = "c")
) %>%
  mutate(relative_abundance = as.numeric(relative_abundance)) %>%
  filter(str_detect(clade_name, "t__")) %>%
  select(clade_name, relative_abundance) %>%
  arrange(-relative_abundance) %>%
  separate(
    clade_name,
    into = c("kingdom", "phylum", "class", "order", "family", "genus", "species", "sgb"),
    sep = "\\|",
    fill = "right"
  ) %>%
  mutate(across(
    .cols = kingdom:sgb,
    .fns = ~ str_remove(.x, "^[a-z]__")  # удаляет префиксы вроде k__, p__, t__ и т.д.
  )) %>%
  rename(
    Царство = kingdom,
    Тип = phylum,
    Класс = class,
    Порядок = order,
    Семейство = family,
    Род = genus,
    Вид = species,
    SGB = sgb,
    `Относительное содержание` = relative_abundance
  )

write_csv2(metaphlan, output_file)

cat("Written to:", output_file, "\n")