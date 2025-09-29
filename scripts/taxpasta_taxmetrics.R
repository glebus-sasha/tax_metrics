#!/usr/bin/env Rscript

input_folder <- "raw/taxpasta/"
output_folder <- "results"

suppressPackageStartupMessages({
  library(tidyverse)
  library(plotly)
})
# 
# args <- commandArgs(trailingOnly = TRUE)
# 
# if (length(args) != 2) {
#   stop("Usage: Rscript taxpasta_taxmetrics.R <input_folder> <output_folder>")
# }
# 
# input_folder <- args[1]
# output_folder <- args[2]

if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# --- Функции метрик ---
calc_binary_metrics <- function(truth_set, pred_set) {
  tp <- length(intersect(truth_set, pred_set))
  fp <- length(setdiff(pred_set, truth_set))
  fn <- length(setdiff(truth_set, pred_set))
  
  precision <- ifelse((tp + fp) > 0, tp / (tp + fp), NA_real_)
  recall    <- ifelse((tp + fn) > 0, tp / (tp + fn), NA_real_)
  f1        <- ifelse(!is.na(precision) & !is.na(recall) & (precision + recall) > 0,
                      2 * precision * recall / (precision + recall),
                      NA_real_)
  
  tibble(precision = precision, recall = recall, f1 = f1)
}

calc_abundance_metrics <- function(truth_ab, pred_ab) {
  taxa_all <- union(names(truth_ab), names(pred_ab))
  truth_vec <- truth_ab[taxa_all]; truth_vec[is.na(truth_vec)] <- 0
  pred_vec  <- pred_ab[taxa_all];  pred_vec[is.na(pred_vec)] <- 0
  
  # нормализация в %
  truth_vec <- as.numeric(truth_vec)
  truth_vec <- 100 * truth_vec / sum(truth_vec)
  pred_vec  <- as.numeric(pred_vec)
  pred_vec  <- 100 * pred_vec / sum(pred_vec)
  
  l1 <- sum(abs(truth_vec - pred_vec))
  pearson <- suppressWarnings(cor(truth_vec, pred_vec, method = "pearson"))
  spearman <- suppressWarnings(cor(truth_vec, pred_vec, method = "spearman"))
  
  tibble(l1_distance = l1,
         pearson_corr = pearson,
         spearman_corr = spearman)
}

# --- Функция для списка совпадений/отсутствующих/лишних ---
get_taxon_lists <- function(truth_set, pred_set) {
  matches <- intersect(truth_set, pred_set)
  missing <- setdiff(truth_set, pred_set)
  extra   <- setdiff(pred_set, truth_set)
  
  tibble(
    matches = paste(matches, collapse = "; "),
    missing = paste(missing, collapse = "; "),
    extra   = paste(extra, collapse = "; "),
    n_matches = length(matches),
    n_missing = length(missing),
    n_extra   = length(extra),
    n_truth   = length(truth_set),
    n_pred    = length(pred_set)
  )
}

# --- Читаем эталон truth.tsv ---
truth_file <- file.path(input_folder, "truth.tsv")
if (!file.exists(truth_file)) stop("Файл truth.tsv не найден в input_folder")

truth <- read_tsv(truth_file, col_types = cols(.default = "c")) %>%
  mutate(across(-c(taxonomy_id, name), ~as.numeric(str_replace(.x, ",", "."))))

# --- Получаем список файлов инструментов ---
files <- list.files(input_folder, pattern = "\\.tsv$", full.names = TRUE)
files <- setdiff(files, truth_file)

all_metrics <- list()
all_taxon_lists <- list()

# --- Обрабатываем каждый файл ---
for (file in files) {
  tool_name <- tools::file_path_sans_ext(basename(file))
  message("Processing: ", tool_name)
  
  df <- read_tsv(file, col_types = cols(.default = "c")) %>%
    mutate(across(-c(taxonomy_id, name), ~as.numeric(str_replace(.x, ",", "."))))
  
  df_sample_cols <- setdiff(names(df), c("taxonomy_id", "name"))
  
  for (df_col in df_sample_cols) {
    # Чистим название sample
    sample_clean <- str_replace(df_col, "_db[0-9]+.*$|\\.bracken$|\\.metaphlan_profile$|\\.metaphlan$", "")
    
    if (!sample_clean %in% names(truth)) {
      message("Skipping sample ", df_col, " — not found in truth")
      next
    }
    
    pred_sample <- df %>%
      select(taxonomy_id, abundance = all_of(df_col)) %>%
      filter(!is.na(abundance), abundance > 0)
    
    truth_sample <- truth %>%
      select(taxonomy_id, abundance = all_of(sample_clean)) %>%
      filter(!is.na(abundance), abundance > 0)
    
    # --- Метрики ---
    bin <- calc_binary_metrics(
      truth_set = truth_sample$taxonomy_id,
      pred_set  = pred_sample$taxonomy_id
    )
    
    abund <- calc_abundance_metrics(
      truth_ab = setNames(truth_sample$abundance, truth_sample$taxonomy_id),
      pred_ab  = setNames(pred_sample$abundance, pred_sample$taxonomy_id)
    )
    
    all_metrics[[length(all_metrics) + 1]] <- bind_cols(
      tibble(tool = tool_name, sample = sample_clean),
      bin,
      abund
    )
    
    # --- Списки совпадений/отсутствующих/лишних ---
    all_taxon_lists[[length(all_taxon_lists) + 1]] <- bind_cols(
      tibble(tool = tool_name, sample = sample_clean),
      get_taxon_lists(truth_sample$taxonomy_id, pred_sample$taxonomy_id)
    )
  }
}

# --- Сохраняем CSV метрик ---
metrics_df <- bind_rows(all_metrics) %>% filter(tool != "truth")
output_file <- file.path(output_folder, "all_taxmetrics.csv")
write_excel_csv2(metrics_df, output_file)
message("All metrics written to: ", output_file)

# --- Сохраняем CSV со списками таксонов ---
taxon_lists_df <- bind_rows(all_taxon_lists)
output_taxon_file <- file.path(output_folder, "taxon_lists.csv")
write_excel_csv2(taxon_lists_df, output_taxon_file)
message("Taxon lists written to: ", output_taxon_file)

# --- Виолин-плоты PNG ---
metrics_long <- metrics_df %>%
  pivot_longer(cols = c(precision, recall, f1, l1_distance, pearson_corr, spearman_corr),
               names_to = "metric", values_to = "value")

plot_png <- file.path(output_folder, "taxmetrics_violin.png")
png(plot_png, width = 1400, height = 900, res = 150)
p <- ggplot(metrics_long, aes(x = tool, y = value, fill = tool)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.7) +
  facet_wrap(~ metric, scales = "free_y") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(title = "Сравнение методов таксономии по метрикам", y = "Значение", x = "Метод")
print(p)
dev.off()
message("Violin plot PNG saved to: ", plot_png)

# --- Интерактивный plotly ---
plot_html <- file.path(output_folder, "taxmetrics_violin.html")
p_plotly <- ggplotly(p)
htmlwidgets::saveWidget(p_plotly, plot_html)
message("Interactive violin plot saved to: ", plot_html)
