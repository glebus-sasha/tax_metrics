#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(plotly)
  library(htmlwidgets)
  library(taxize)
  library(writexl)
})


# csv_dir <- "."
# truth_name <- "truth"
# output_prefix <- "sample"

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Usage: Rscript compare_all_vs_truth.R <csv_dir> <truth_name> <output_prefix>")
}

csv_dir <- args[1]
truth_name <- args[2]
output_prefix <- args[3]

# Создаем структуру папок для результатов
results_dir <- paste0(output_prefix, "_results")
plots_all_dir <- file.path(results_dir, "plots_all")
plots_by_rank_dir <- file.path(results_dir, "plots_by_rank")
tables_dir <- file.path(results_dir, "tables")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_all_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_by_rank_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

# ───────────────────────────── 1. Загрузка профилей ─────────────────────────────
load_profiles <- function(folder = ".", truth_name = "truth") {
  files <- list.files(folder, pattern = "\\.csv$", full.names = TRUE)
  
  profiles <- files %>%
    set_names(~ basename(.x) %>% tools::file_path_sans_ext() %>% str_split("_") %>% map_chr(~ tail(.x, 1))) %>%
    map(~ read_csv2(.x) %>% mutate(abundance = 100 * abundance / sum(abundance, na.rm = TRUE)))
  
  if (!(truth_name %in% names(profiles))) {
    stop("No 'truth' profile found among files.")
  }
  
  profiles
}

# ───────────────────────────── 2. Метрики присутствия ───────────────────────────
compare_with_truth <- function(truth_df, test_df) {
  truth_taxa <- truth_df %>% filter(abundance > 0) %>% pull(tax_id) %>% unique()
  test_taxa  <- test_df %>% filter(abundance > 0) %>% pull(tax_id) %>% unique()
  
  tp <- length(intersect(truth_taxa, test_taxa))
  fn <- length(setdiff(truth_taxa, test_taxa))
  fp <- length(setdiff(test_taxa, truth_taxa))
  
  precision <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
  recall    <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  f1        <- if (!is.na(precision) && !is.na(recall) && (precision + recall) > 0) {
    2 * (precision * recall) / (precision + recall)
  } else {
    NA_real_
  }
  
  tibble(
    TP = tp,
    FP = fp,
    FN = fn,
    total_truth = length(truth_taxa),
    total_test = length(test_taxa),
    precision = round(precision, 3),
    recall = round(recall, 3),
    F1 = round(f1, 3)
  )
}

# ───────────────────────────── 3. Метрики abundance ─────────────────────────────
compare_abundance <- function(truth_df, test_df) {
  merged <- full_join(truth_df, test_df, by = "tax_id", suffix = c("_truth", "_test")) %>%
    replace_na(list(abundance_truth = 0, abundance_test = 0))
  
  l1_dist <- sum(abs(merged$abundance_truth - merged$abundance_test))
  l2_dist <- sqrt(sum((merged$abundance_truth - merged$abundance_test)^2))
  spearman <- suppressWarnings(cor(merged$abundance_truth, merged$abundance_test, method = "spearman"))
  pearson  <- suppressWarnings(cor(merged$abundance_truth, merged$abundance_test, method = "pearson"))
  
  tibble(
    l1 = round(l1_dist, 4),
    l2 = round(l2_dist, 4),
    spearman = round(spearman, 4),
    pearson = round(pearson, 4)
  )
}

# ───────────────────────────── 4. Основной запуск ───────────────────────────────
profiles <- load_profiles(csv_dir, truth_name)
truth <- profiles[[truth_name]]

abundance_combined <- map2_dfr(profiles, names(profiles), ~ mutate(.x, tool = .y))

# ──────────────────────── 5. Классификация tax_id → ранг ────────────────────────
summary_taxa <- abundance_combined %>%
  group_by(tax_id, tool) %>%
  summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
  mutate(tax_id = as.character(tax_id))

class_list <- tryCatch(
  classification(unique(summary_taxa$tax_id), db = "ncbi"),
  error = function(e) {
    warning("Ошибка получения классификации: ", e$message)
    list()
  }
)

names(class_list) <- unique(summary_taxa$tax_id)

get_lowest_info <- function(x) {
  if (is.null(x) || inherits(x, "try-error") || length(x) == 0 ||
      !is.data.frame(x) || !all(c("name", "rank") %in% colnames(x))) {
    return(tibble(scientific_name = NA_character_, rank = NA_character_))
  }
  last_row <- tail(x, 1)
  tibble(scientific_name = last_row$name, rank = last_row$rank)
}

tax_names_df <- purrr::map2_dfr(names(class_list), class_list, ~ {
  get_lowest_info(.y) %>% mutate(tax_id = .x)
})

summary_taxa_annotated <- summary_taxa %>%
  left_join(tax_names_df, by = "tax_id") %>%
  relocate(scientific_name, rank, .after = tax_id)

# ───────────────────────────── 6. Расчёт метрик ────────────────────────────────
calculate_metrics_for_df <- function(df, truth_name) {
  tools <- unique(df$tool)
  truth <- df %>% filter(tool == truth_name)
  
  map_dfr(setdiff(tools, truth_name), function(tool_name) {
    test <- df %>% filter(tool == tool_name)
    presence <- compare_with_truth(truth, test)
    abundance <- compare_abundance(truth, test)
    bind_cols(tool = tool_name, presence, abundance)
  })
}

# Все ранги
metrics_all <- calculate_metrics_for_df(summary_taxa_annotated, truth_name)

# Метрики по каждому рангу
ranks <- summary_taxa_annotated %>%
  filter(!is.na(rank)) %>%
  pull(rank) %>%
  unique()

metrics_by_rank <- set_names(ranks) %>%
  map(~ summary_taxa_annotated %>% filter(rank == .x) %>% calculate_metrics_for_df(truth_name))

# ─────────────── 7. Сохранение метрик в Excel (разные листы по рангам) ───────────────
excel_out <- file.path(tables_dir, "_metrics_by_rank.xlsx")
write_xlsx(c("all" = list(metrics_all), metrics_by_rank), path = excel_out)

# ─────────────── 8. Построение графиков для каждого ранга + общий ───────────────
plot_metric_barplot <- function(df, suffix) {
  df_long <- df %>%
    pivot_longer(cols = c(precision, recall, F1), names_to = "metric", values_to = "value")
  
  ggplot(df_long, aes(x = tool, y = value, fill = metric)) +
    geom_col(position = "dodge") +
    scale_fill_brewer(palette = "Set2") +
    labs(title = paste("Метрики точности по рангу:", suffix),
         x = "Инструмент", y = "Значение") +
    ylim(0, 1) +
    theme_bw()
}

# График для всех рангов
plot_all <- plot_metric_barplot(metrics_all, "all")
ggsave(file.path(plots_all_dir, "_metrics_all.png"), plot_all, width = 8, height = 5)
saveWidget(ggplotly(plot_all), file = file.path(plots_all_dir, "_metrics_all.html"), selfcontained = TRUE)

# Графики по рангам
for (r in ranks) {
  p <- plot_metric_barplot(metrics_by_rank[[r]], r)
  png_name <- file.path(plots_by_rank_dir, paste0("_metrics_rank_", r, ".png"))
  html_name <- file.path(plots_by_rank_dir, paste0("_metrics_rank_", r, ".html"))
  ggsave(png_name, p, width = 8, height = 5)
  saveWidget(ggplotly(p), file = html_name, selfcontained = TRUE)
}

# ─────────────── 9. Сохранение финальной таблицы с аннотацией ───────────────
write_csv2(summary_taxa_annotated, file.path(tables_dir, "_summary_detected_taxa_annotated.csv"))

# ─────────────── 10. Завершение ───────────────
cat("✔ Все результаты сохранены:\n")
cat(" - Метрики (все ранги):        ", file.path(plots_all_dir, "_metrics_all.png/html"), "\n")
cat(" - Excel-файл с метриками по рангам: ", excel_out, "\n")
cat(" - Аннотированная таблица:     ", file.path(tables_dir, "_summary_detected_taxa_annotated.csv"), "\n")
cat(" - Графики по рангам PNG/HTML сохранены в папке: ", plots_by_rank_dir, "\n")
