#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

# metaphlan_file <- 'raw/250822_RnD-L_250822_21_Metagenom1_n17K_L00_profile.txt'
# output_file <- 'results/250822_RnD-L_250822_21_Metagenom1_n17K_L00_taxonomy.csv'

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript metaphlan_results.R <metaphlan_file.tsv> <output_file.csv>")
}

metaphlan_file <- args[1]
output_file <- args[2]

metaphlan <- read_tsv(
  metaphlan_file,
  comment = "#",
  col_names = c("clade_name", "clade_taxid", "relative_abundance", "coverage", "estimated_number_of_reads_from_the_clade"),
  col_types = cols(.default = "c")
) %>%
  mutate(relative_abundance = as.numeric(relative_abundance)) %>%
  filter(str_detect(clade_name, "t__")) %>%
  select(clade_name, estimated_number_of_reads_from_the_clade, relative_abundance) %>%
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
  mutate(
    total_filtered_abundance = sum(relative_abundance),
    relative_abundance = (relative_abundance / total_filtered_abundance) * 100
  ) %>% 
  rename(
    Царство = kingdom,
    Тип = phylum,
    Класс = class,
    Порядок = order,
    Семейство = family,
    Род = genus,
    Вид = species,
    SGB = sgb,
    `Количество прочтений` = estimated_number_of_reads_from_the_clade,
    `Относительное содержание, %` = relative_abundance
  ) %>% 
  select(-SGB, -total_filtered_abundance) %>% 
  filter(`Относительное содержание, %` > 0.01) %>% 
  arrange(desc(`Относительное содержание, %`))

write_excel_csv2(metaphlan, output_file)

cat("Written to:", output_file, "\n")