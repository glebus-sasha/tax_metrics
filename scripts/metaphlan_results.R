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
  mutate(estimated_number_of_reads_from_the_clade = as.numeric(estimated_number_of_reads_from_the_clade)) %>%
  filter(str_detect(clade_name, "t__")) %>%
  select(clade_taxid, clade_name, estimated_number_of_reads_from_the_clade, relative_abundance) %>%
  # Разделяем таксономию и taxid по уровням
  separate(
    clade_name,
    into = c("kingdom", "phylum", "class", "order", "family", "genus", "species", "strain"),
    sep = "\\|",
    fill = "right"
  ) %>%
  separate(
    clade_taxid,
    into = c("kingdom_taxid", "phylum_taxid", "class_taxid", "order_taxid", "family_taxid", "genus_taxid", "species_taxid", "strain_taxid"),
    sep = "\\|",
    fill = "right"
  ) %>%
  # Удаляем префиксы из названий таксонов
  mutate(across(
    .cols = kingdom:strain,
    .fns = ~ str_remove(.x, "^[a-z]__")
  )) %>%
  # Заменяем пустые значения на NA
  mutate(across(
    .cols = kingdom:strain,
    .fns = ~ ifelse(.x == "" | is.na(.x), NA, .x)
  )) %>%
  mutate(across(
    .cols = kingdom_taxid:strain_taxid,
    .fns = ~ ifelse(.x == "" | is.na(.x), NA, .x)
  )) %>%
  # Пересчитываем относительную abundance
  mutate(
    total_filtered_abundance = sum(estimated_number_of_reads_from_the_clade),
    relative_abundance = (estimated_number_of_reads_from_the_clade / total_filtered_abundance) * 100,
    relative_abundance = round(relative_abundance, 3)
  ) %>% 
  # Переименовываем столбцы
  rename(
    Царство = kingdom,
    Тип = phylum,
    Класс = class,
    Порядок = order,
    Семейство = family,
    Род = genus,
    Вид = species,
    `Количество прочтений` = estimated_number_of_reads_from_the_clade,
    `Относительное содержание %` = relative_abundance
  ) %>% 
  # Добавляем taxid для каждого уровня с русскими названиями
  rename(
    `TaxID царства` = kingdom_taxid,
    `TaxID типа` = phylum_taxid,
    `TaxID класса` = class_taxid,
    `TaxID порядка` = order_taxid,
    `TaxID семейства` = family_taxid,
    `TaxID рода` = genus_taxid,
    `TaxID вида` = species_taxid
  ) %>% 
  # Убираем ненужные столбцы (штамм/MAG и их taxid)
  select(-strain, -strain_taxid, -total_filtered_abundance) %>% 
  # Переупорядочиваем столбцы для удобства
  select(
   `TaxID царства`, `TaxID типа`, `TaxID класса`, `TaxID порядка`, 
    `TaxID семейства`, `TaxID рода`, `TaxID вида`,
    Царство, Тип, Класс, Порядок, Семейство, Род, Вид,
    `Количество прочтений`, `Относительное содержание %`
  ) %>% 
  arrange(desc(`Относительное содержание %`))

write_excel_csv2(metaphlan, output_file)

cat("Written to:", output_file, "\n")