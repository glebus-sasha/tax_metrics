#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

# input_file <- 'raw/sample_0_taxonomy.csv'
# output_file <- 'results/sample_0_taxonomy.csv'

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript bracken_results.R <kreport2mpa_file.tsv> <output_file.csv>")
}

input_file <- args[1]
output_file <- args[2]

# Чтение данных
kreport_data <- read_tsv(
  input_file,
  col_names = c("taxonomy_path", "reads_count"),
  col_types = cols(.default = "c")
) %>%
  mutate(reads_count = as.numeric(reads_count)) %>%
  filter(!is.na(reads_count), reads_count > 0)

# Функция для извлечения последнего таксономического уровня
extract_last_taxon <- function(path) {
  taxa <- unlist(str_split(path, "\\|"))
  last_taxon <- tail(taxa, 1)
  return(last_taxon)
}

# Функция для извлечения таксономического уровня по префиксу
extract_taxon_level <- function(path, level_prefix) {
  pattern <- paste0(level_prefix, "__([^|]+)")
  match <- str_extract(path, pattern)
  if (!is.na(match)) {
    return(str_remove(match, paste0(level_prefix, "__")))
  } else {
    return(NA_character_)
  }
}

# Обработка данных
processed_data <- kreport_data %>%
  # Извлекаем последний таксон (самый специфичный уровень)
  mutate(last_taxon = map_chr(taxonomy_path, extract_last_taxon)) %>%
  # Извлекаем отдельные таксономические уровни
  mutate(
    kingdom = map_chr(taxonomy_path, ~extract_taxon_level(.x, "k")),
    phylum = map_chr(taxonomy_path, ~extract_taxon_level(.x, "p")),
    class = map_chr(taxonomy_path, ~extract_taxon_level(.x, "c")),
    order = map_chr(taxonomy_path, ~extract_taxon_level(.x, "o")),
    family = map_chr(taxonomy_path, ~extract_taxon_level(.x, "f")),
    genus = map_chr(taxonomy_path, ~extract_taxon_level(.x, "g")),
    species = map_chr(taxonomy_path, ~extract_taxon_level(.x, "s"))
  ) %>%
  filter(!is.na(species)) %>%
  filter(species != "Homo_sapiens") %>% 
  mutate(total_reads = sum(reads_count)) %>%
  mutate(relative_abundance = (reads_count / total_reads) * 100) %>%
  select(
    Царство = kingdom,
    Тип = phylum,
    Класс = class,
    Порядок = order,
    Семейство = family,
    Род = genus,
    Вид = species,
    `Количество прочтений` = reads_count,
    `Относительное содержание` = relative_abundance
  ) %>%
  mutate(across(everything(), ~replace_na(.x, ""))) %>%
  filter(`Относительное содержание` > 0.01) %>% 
  arrange(desc(`Относительное содержание`))

# Сохраняем результат
write_excel_csv2(processed_data, output_file)

cat("Written to:", output_file, "\n")
cat("Total rows processed:", nrow(processed_data), "\n")