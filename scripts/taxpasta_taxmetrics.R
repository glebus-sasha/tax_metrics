#!/usr/bin/env Rscript

# args <- commandArgs(trailingOnly = TRUE)
# 
# if (length(args) != 2) {
#   stop("Usage: Rscript taxpasta_taxmetrics.R <input_folder> <output_folder>")
# }
# 
# input_folder <- args[1]
# output_folder <- args[2]

input_folder <- "raw/tools/"
output_folder <- "results"

suppressPackageStartupMessages({
  library(tidyverse)
  library(yardstick)
  library(writexl)
  library(plotly)
  library(openxlsx)
})

# Функция для нормализации (приведение к относительной abundance)
normalize_columns <- function(df) {
  df %>%
    mutate(across(-c(taxonomy_id, name), 
                  ~ .x / sum(.x, na.rm = TRUE) * 100))
}

# Создаем выходную папку если не существует
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

# Находим все файлы в папке (кроме truth.tsv)
all_files <- list.files(input_folder, pattern = "\\.tsv$", full.names = TRUE)
truth_file <- file.path(input_folder, "truth.tsv")
test_files <- all_files[!grepl("truth\\.tsv$", all_files)]

# Получаем имена инструментов из названий файлов
tool_names <- tools::file_path_sans_ext(basename(test_files))

# Читаем истинные данные (только виды)
truth <- read_delim(truth_file, delim = "\t", escape_double = FALSE, trim_ws = TRUE) %>% 
  normalize_columns() %>%
  rename_with(~ ifelse(.x == "taxonomy_id",
                       .x,
                       str_replace(.x, "^([^_]+_[^_]+_[^_]+_[^_]+).*", "\\1")),
              .cols = everything()) %>% 
  pivot_longer(cols = -c(taxonomy_id, name),
               names_to = "sample",
               values_to = "abundance")

# Создаем списки для хранения результатов
all_metrics <- list()
all_taxonomy <- list()
all_taxonomy_comparison <- list()  # Новый список для сравнения таксономии
comparison_data <- data.frame()

# Обрабатываем каждый файл с результатами инструментов
for (i in seq_along(test_files)) {
  tool_file <- test_files[i]
  tool_name <- tool_names[i]
  
  cat("Processing:", tool_name, "\n")
  
  # Читаем и обрабатываем данные инструмента (только виды)
  test_data <- read_delim(tool_file, delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  
  # Фильтруем только виды, если есть столбец rank
  if ("rank" %in% colnames(test_data)) {
    test_filtered <- test_data %>% 
      filter(rank == 'species') %>%
      select(-rank)
  } else {
    test_filtered <- test_data
  }
  
  test <- test_filtered %>%
    rename_with(~ ifelse(.x == "taxonomy_id",
                         .x,
                         str_replace(.x, "^([^_]+_[^_]+_[^_]+_[^_]+).*", "\\1")),
                .cols = everything()) %>% 
    normalize_columns() %>%
    pivot_longer(cols = -c(taxonomy_id, name),
                 names_to = "sample",
                 values_to = "abundance")
  
  # Сохраняем таксономические данные для Excel (только виды)
  all_taxonomy[[tool_name]] <- test_filtered
  
  # Создаем таблицу сравнения таксономии (найденные, ненайденные, ошибочные) - только виды
  # Используем только taxonomy_id для сравнения, name берем из truth для удобства
  truth_taxa <- truth %>% 
    filter(abundance > 0) %>% 
    distinct(taxonomy_id) %>%
    mutate(present_in_truth = TRUE)
  
  test_taxa <- test %>% 
    filter(abundance > 0) %>% 
    distinct(taxonomy_id) %>%
    mutate(present_in_test = TRUE)
  
  # Берем названия из истинных данных для удобства
  truth_names <- truth %>% 
    distinct(taxonomy_id, name) %>%
    rename(name_truth = name)
  
  taxonomy_comparison <- truth_taxa %>%
    full_join(test_taxa, by = "taxonomy_id") %>%
    left_join(truth_names, by = "taxonomy_id") %>%
    mutate(
      present_in_truth = replace_na(present_in_truth, FALSE),
      present_in_test = replace_na(present_in_test, FALSE),
      status = case_when(
        present_in_truth & present_in_test ~ "Correctly identified",
        present_in_truth & !present_in_test ~ "Missed",
        !present_in_truth & present_in_test ~ "False positive",
        TRUE ~ "Not present"
      )
    ) %>%
    filter(status != "Not present") %>%  # Убираем таксоны, которых нет нигде
    select(taxonomy_id, name_truth, status) %>%
    rename(name = name_truth) %>%
    arrange(status, taxonomy_id)
  
  all_taxonomy_comparison[[tool_name]] <- taxonomy_comparison
  
  # Объединяем с истинными данными - используем только taxonomy_id
  combined_data <- truth %>% 
    select(taxonomy_id, sample, abundance_truth = abundance) %>%
    full_join(test %>% 
                select(taxonomy_id, sample, abundance_test = abundance), 
              by = c("taxonomy_id", "sample")) %>% 
    mutate(
      abundance_truth = case_when(
        !is.na(abundance_truth) ~ abundance_truth,
        TRUE ~ 0
      ),
      abundance_test = case_when(
        !is.na(abundance_test) ~ abundance_test,
        TRUE ~ 0
      ),
      presence_truth = as.factor(case_when(
        abundance_truth > 0 ~ 1,
        TRUE ~ 0
      )),
      presence_test = as.factor(case_when(
        abundance_test > 0 ~ 1,
        TRUE ~ 0
      ))
    )
  
  # Вычисляем бинарные метрики
  binary_metrics <- combined_data %>% 
    group_by(sample) %>%
    summarise(
      # Базовые счетчики
      tp = sum(presence_truth == 1 & presence_test == 1),
      tn = sum(presence_truth == 0 & presence_test == 0),
      fp = sum(presence_truth == 0 & presence_test == 1),
      fn = sum(presence_truth == 1 & presence_test == 0),
      total = n(),
      
      # Accuracy (всегда можно посчитать)
      accuracy = (tp + tn) / total,
      
      # Precision (только если есть предсказанные положительные)
      precision = ifelse((tp + fp) > 0, tp / (tp + fp), NA_real_),
      
      # Recall (только если есть настоящие положительные)
      recall = ifelse((tp + fn) > 0, tp / (tp + fn), NA_real_),
      
      # F1-score (только если можно посчитать и precision и recall)
      f1 = ifelse(!is.na(precision) & !is.na(recall) & (precision + recall) > 0,
                  2 * (precision * recall) / (precision + recall), NA_real_),
      
      # Specificity (только если есть настоящие отрицательные)
      specificity = ifelse((tn + fp) > 0, tn / (tn + fp), NA_real_),
      
      .groups = 'drop'
    )
  
  # Вычисляем дистанционные метрики
  distance_metrics <- combined_data %>%
    group_by(sample) %>%
    summarise(
      # Bray-Curtis (самая популярная в метагеномике)
      bray_curtis = 1 - (2 * sum(pmin(abundance_truth, abundance_test)) / 
                           (sum(abundance_truth) + sum(abundance_test))),
      
      # Jaccard (для presence/absence с abundance)
      jaccard = sum(pmin(abundance_truth, abundance_test)) / 
        sum(pmax(abundance_truth, abundance_test)),
      
      # Cosine similarity
      cosine = sum(abundance_truth * abundance_test) / 
        (sqrt(sum(abundance_truth^2)) * sqrt(sum(abundance_test^2))),
      
      # Manhattan distance
      manhattan = sum(abs(abundance_truth - abundance_test)),
      
      # Euclidean distance
      euclidean = sqrt(sum((abundance_truth - abundance_test)^2)),
      
      .groups = 'drop'
    )
  
  # Объединяем все метрики
  tool_metrics <- binary_metrics %>%
    full_join(distance_metrics, by = "sample") %>%
    mutate(tool = tool_name) %>%
    select(tool, sample, everything())
  
  # Добавляем в общий список
  all_metrics[[tool_name]] <- tool_metrics
  
  # Собираем данные для сравнения инструментов (усредняем по образцам)
  avg_metrics <- tool_metrics %>%
    summarise(across(c(accuracy, precision, recall, f1, specificity, 
                       bray_curtis, jaccard, cosine, manhattan, euclidean), 
                     mean, na.rm = TRUE)) %>%
    mutate(tool = tool_name)
  
  comparison_data <- bind_rows(comparison_data, avg_metrics)
}

# 1. Сохраняем метрики в Excel с правильным порядком (без truth)
metrics_output_file <- file.path(output_folder, "taxpasta_metrics_comparison.xlsx")

# Создаем workbook и добавляем листы в правильном порядке
wb <- createWorkbook()

# Лист с сравнением инструментов (первый и самый важный)
addWorksheet(wb, "Tools Comparison")
writeData(wb, "Tools Comparison", comparison_data)

# Лист с истинными данными
truth_data <- read_delim(truth_file, delim = "\t", escape_double = FALSE, trim_ws = TRUE)
addWorksheet(wb, "Ground Truth")
writeData(wb, "Ground Truth", truth_data)

# Листы с метриками для каждого инструмента (только тестовые инструменты, не truth)
for (tool_name in tool_names) {
  if (tool_name != "truth") {
    addWorksheet(wb, paste0("Metrics_", tool_name))
    writeData(wb, paste0("Metrics_", tool_name), all_metrics[[tool_name]])
  }
}

saveWorkbook(wb, metrics_output_file, overwrite = TRUE)

# 2. Сохраняем таксономические данные в отдельный Excel с таблицей сравнения
taxonomy_output_file <- file.path(output_folder, "taxpasta_taxonomy_results.xlsx")

taxonomy_wb <- createWorkbook()

# Добавляем исходные данные инструментов (только виды)
for (tool_name in tool_names) {
  addWorksheet(taxonomy_wb, tool_name)
  writeData(taxonomy_wb, tool_name, all_taxonomy[[tool_name]])
}

# Добавляем таблицы сравнения таксономии
for (tool_name in tool_names) {
  if (tool_name != "truth") {
    comparison_sheet_name <- paste0("Comparison_", tool_name)
    addWorksheet(taxonomy_wb, comparison_sheet_name)
    writeData(taxonomy_wb, comparison_sheet_name, all_taxonomy_comparison[[tool_name]])
    
    # Добавляем сводную статистику
    summary_stats <- all_taxonomy_comparison[[tool_name]] %>%
      count(status) %>%
      mutate(percentage = n / sum(n) * 100)
    
    addWorksheet(taxonomy_wb, paste0("Summary_", tool_name))
    writeData(taxonomy_wb, paste0("Summary_", tool_name), summary_stats)
  }
}

saveWorkbook(taxonomy_wb, taxonomy_output_file, overwrite = TRUE)

# 3. Создаем интерактивные графики для сравнения инструментов
create_comparison_plots <- function(comparison_data, output_folder) {
  
  # Убираем truth из данных для графиков
  plot_data <- comparison_data %>% filter(tool != "truth")
  
  # График 1: Бинарные метрики классификации
  binary_plot <- plot_data %>%
    plot_ly() %>%
    add_trace(x = ~tool, y = ~accuracy, type = 'bar', name = 'Accuracy',
              marker = list(color = '#1f77b4')) %>%
    add_trace(x = ~tool, y = ~precision, type = 'bar', name = 'Precision',
              marker = list(color = '#ff7f0e')) %>%
    add_trace(x = ~tool, y = ~recall, type = 'bar', name = 'Recall',
              marker = list(color = '#2ca02c')) %>%
    add_trace(x = ~tool, y = ~f1, type = 'bar', name = 'F1-Score',
              marker = list(color = '#d62728')) %>%
    add_trace(x = ~tool, y = ~specificity, type = 'bar', name = 'Specificity',
              marker = list(color = '#9467bd')) %>%
    layout(title = 'Binary Classification Metrics by Tool',
           xaxis = list(title = 'Tool'),
           yaxis = list(title = 'Score', range = c(0, 1)),
           barmode = 'group')
  
  # График 2: Дистанционные метрики (сходство)
  similarity_plot <- plot_data %>%
    plot_ly() %>%
    add_trace(x = ~tool, y = ~jaccard, type = 'bar', name = 'Jaccard',
              marker = list(color = '#8c564b')) %>%
    add_trace(x = ~tool, y = ~cosine, type = 'bar', name = 'Cosine',
              marker = list(color = '#e377c2')) %>%
    layout(title = 'Similarity Metrics by Tool',
           xaxis = list(title = 'Tool'),
           yaxis = list(title = 'Similarity Score', range = c(0, 1)),
           barmode = 'group')
  
  # График 3: Дистанционные метрики (расстояния)
  distance_plot <- plot_data %>%
    plot_ly() %>%
    add_trace(x = ~tool, y = ~bray_curtis, type = 'bar', name = 'Bray-Curtis',
              marker = list(color = '#7f7f7f')) %>%
    add_trace(x = ~tool, y = ~manhattan, type = 'bar', name = 'Manhattan',
              marker = list(color = '#bcbd22')) %>%
    add_trace(x = ~tool, y = ~euclidean, type = 'bar', name = 'Euclidean',
              marker = list(color = '#17becf')) %>%
    layout(title = 'Distance Metrics by Tool',
           xaxis = list(title = 'Tool'),
           yaxis = list(title = 'Distance'),
           barmode = 'group')
  
  # Сохраняем все графики в один HTML файл как отдельные виджеты
  html_file <- file.path(output_folder, "tools_comparison_plots.html")
  
  # Создаем комбинированный HTML с отдельными графиками
  combined_html <- htmltools::tagList(
    htmltools::tags$h1("Comprehensive Tool Comparison Dashboard"),
    
    htmltools::tags$h2("Binary Classification Metrics"),
    binary_plot,
    
    htmltools::tags$h2("Similarity Metrics"),
    similarity_plot,
    
    htmltools::tags$h2("Distance Metrics"), 
    distance_plot
  )
  
  htmltools::save_html(combined_html, html_file)
  return(html_file)
}

# Создаем графики (только для тестовых инструментов, без truth)
html_plot_file <- create_comparison_plots(comparison_data, output_folder)

cat("Results saved to:\n")
cat("  Metrics comparison:", metrics_output_file, "\n")
cat("  Taxonomy results:", taxonomy_output_file, "\n")
cat("  Interactive plots:", html_plot_file, "\n")
cat("Tools processed:", paste(tool_names[tool_names != "truth"], collapse = ", "), "\n")