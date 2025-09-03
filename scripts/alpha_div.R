#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(jsonlite)
})

# taxonomy_file <- 'results/sample_0_taxonomy.csv'
# output_file <- 'results/sample_0_alpha.json'

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript alpha_div.R <input_file.csv> <output_file.csv>")
}

bracken_file <- args[1]
output_file <- args[2]

taxonomy_df <- read_csv2(taxonomy_file) 

counts <- taxonomy_df %>%
  pull(`Количество прочтений`) %>% 
  as.numeric() %>%
  discard(is.na)

genus_count <- taxonomy_df %>%
  distinct(Род) %>%
  filter(!is.na(Род), Род != "") %>%
  nrow()

# Расчет соотношений бактерий
# 1. Соотношение Bacillota/Bacteroidota (типы)
firmicutes_count <- taxonomy_df %>%
  filter(Тип == "Bacillota" | Тип == "Firmicutes") %>%
  pull(`Количество прочтений`) %>%
  sum(na.rm = TRUE)

bacteroidota_count <- taxonomy_df %>%
  filter(Тип == "Bacteroidota" | Тип == "Bacteroidetes") %>%
  pull(`Количество прочтений`) %>%
  sum(na.rm = TRUE)

firmicutes_bacteroidota_ratio <- ifelse(bacteroidota_count > 0, 
                                        firmicutes_count / bacteroidota_count, 
                                        NA)

# 2. Соотношение Bacteroides fragilis/Faecalibacterium prausnitzii (виды)
b_fragilis_count <- taxonomy_df %>%
  filter(Вид == "Bacteroides_fragilis" | Вид == "Bacteroides fragilis") %>%
  pull(`Количество прочтений`) %>%
  sum(na.rm = TRUE)

f_prausnitzii_count <- taxonomy_df %>%
  filter(Вид == "Faecalibacterium_prausnitzii" | Вид == "Faecalibacterium prausnitzii") %>%
  pull(`Количество прочтений`) %>%
  sum(na.rm = TRUE)

b_fragilis_f_prausnitzii_ratio <- ifelse(f_prausnitzii_count > 0, 
                                         b_fragilis_count / f_prausnitzii_count, 
                                         NA)

# 3. Расчет доминирующих родов для определения энтеротипа
dominant_genera <- taxonomy_df %>%
  filter(!is.na(Род), Род != "") %>%
  group_by(Род) %>%
  summarise(total_count = sum(`Количество прочтений`, na.rm = TRUE)) %>%
  arrange(desc(total_count)) %>%
  mutate(relative_abundance = total_count / sum(total_count))

# Определение энтеротипа
bacteroides_abundance <- dominant_genera %>%
  filter(Род == "Bacteroides") %>%
  pull(relative_abundance) %>%
  ifelse(length(.) == 0, 0, .)

prevotella_abundance <- dominant_genera %>%
  filter(Род == "Prevotella") %>%
  pull(relative_abundance) %>%
  ifelse(length(.) == 0, 0, .)

ruminococcus_abundance <- dominant_genera %>%
  filter(Род == "Ruminococcus") %>%
  pull(relative_abundance) %>%
  ifelse(length(.) == 0, 0, .)

# Правило определения энтеротипа
entero_type <- case_when(
  bacteroides_abundance > prevotella_abundance & bacteroides_abundance > ruminococcus_abundance ~ "Bacteroides",
  prevotella_abundance > bacteroides_abundance & prevotella_abundance > ruminococcus_abundance ~ "Prevotella",
  ruminococcus_abundance > bacteroides_abundance & ruminococcus_abundance > prevotella_abundance ~ "Ruminococcus",
  TRUE ~ "Mixed"
)

# 4. Дополнительные важные показатели
# Общее количество прочтений
total_reads <- sum(counts, na.rm = TRUE)

# Процент основных типов
percent_firmicutes <- (firmicutes_count / total_reads) * 100
percent_bacteroidota <- (bacteroidota_count / total_reads) * 100
percent_actinobacteria <- taxonomy_df %>%
  filter(Тип == "Actinobacteria" | Тип == "Actinomycetota") %>%
  pull(`Количество прочтений`) %>%
  sum(na.rm = TRUE) / total_reads * 100

percent_proteobacteria <- taxonomy_df %>%
  filter(Тип == "Proteobacteria" | Тип == "Pseudomonadota") %>%
  pull(`Количество прочтений`) %>%
  sum(na.rm = TRUE) / total_reads * 100

# Считаем все нужные метрики альфа-разнообразия
alpha_df <- tibble(
  # Базовая информация
  total_reads = total_reads,
  
  # Индексы биоразнообразия
  num_species = specnumber(counts),
  num_genus = genus_count,
  shannon  = diversity(counts, index = "shannon"),
  pielou   = ifelse(num_species > 1, shannon / log(num_species), NA), 
  simpson  = diversity(counts, index = "simpson"),
  invsimpson = diversity(counts, index = "invsimpson"),
  
  # Соотношения бактерий
  firmicutes_bacteroidota_ratio = firmicutes_bacteroidota_ratio,
  b_fragilis_f_prausnitzii_ratio = b_fragilis_f_prausnitzii_ratio,
  
  # Энтеротип
  entero_type = entero_type,
  bacteroides_abundance = bacteroides_abundance,
  prevotella_abundance = prevotella_abundance,
  ruminococcus_abundance = ruminococcus_abundance,
  
  # Распределение по типам (%)
  percent_firmicutes = percent_firmicutes,
  percent_bacteroidota = percent_bacteroidota,
  percent_actinobacteria = percent_actinobacteria,
  percent_proteobacteria = percent_proteobacteria,
  
  # Дополнительная информация для проверки
  firmicutes_count = firmicutes_count,
  bacteroidota_count = bacteroidota_count,
  b_fragilis_count = b_fragilis_count,
  f_prausnitzii_count = f_prausnitzii_count
)

# Сохраняем
write_json(alpha_df, output_file, pretty = TRUE)

cat("Written to:", output_file, "\n")
cat("Энтеротип:", entero_type, "\n")
cat("Соотношение F/B:", round(firmicutes_bacteroidota_ratio, 2), "\n")
cat("Соотношение B.frag/F.prau:", round(b_fragilis_f_prausnitzii_ratio, 2), "\n")
cat("Индекс Шеннона:", round(alpha_df$shannon, 2), "\n")
cat("Количество родов:", genus_count, "\n")

print(alpha_df)