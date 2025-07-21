# Получаем аргументы из командной строки
 args <- commandArgs(trailingOnly = TRUE)
 if (length(args) != 2) {
  stop("Usage: Rscript compare_abundance.R <truth_file.tsv> <kraken_result.txt> <sid>")
 }
truth_file <- args[1]
kraken2_result <- args[2]
samlpe_name <- args[3]

truth_file <- 'data/short_reads_mapping.tsv'
kraken2_result <- 'data/short_kraken2_result.txt'

truth_abundance <- read.table(truth_file, comment.char = "@", header = T) %>% 
  rename(anonymous_read_id = X.anonymous_read_id) %>%
  group_by(tax_id) %>%
  summarise(reads = n(), .groups = "drop") %>%
  mutate(fraction = reads / sum(reads))
kraken_abundance <- read_delim(kraken2_result, 
                               delim = "\t", escape_double = FALSE, col_names = FALSE,
                               trim_ws = TRUE) %>%
  set_names(c("classification", "read_id", "tax_id", "kmer_hits", "kmer_total")) %>% 
  extract(tax_id, into = c("name", "tax_id"), regex = "^(.*) \\(taxid ([0-9]+)\\)$") %>%
  mutate(tax_id = as.integer(tax_id)) %>%
  group_by(tax_id) %>%
  summarise(reads = n(), .groups = "drop") %>%
  mutate(fraction = reads / sum(reads))

calculate_all_metrics <- function(truth_df, pred_df, tax_col = "tax_id", fraction_col = "fraction") {
  # Слепим датафреймы по tax_id и подставим 0, если taxa нет в одном из наборов
  df_compare <- full_join(
    truth_df %>% select(!!sym(tax_col), truth_fraction = !!sym(fraction_col)),
    pred_df %>% select(!!sym(tax_col), pred_fraction = !!sym(fraction_col)),
    by = tax_col
  ) %>%
    replace_na(list(truth_fraction = 0, pred_fraction = 0))
  
  # Векторы для количественных метрик
  truth_vec <- df_compare$truth_fraction
  pred_vec <- df_compare$pred_fraction
  
  # Метрики для presence/absence
  truth_presence <- df_compare[[tax_col]] %in% truth_df[[tax_col]]
  pred_presence <- df_compare[[tax_col]] %in% pred_df[[tax_col]]
  
  TP <- sum(truth_presence & pred_presence)
  FP <- sum(!truth_presence & pred_presence)
  FN <- sum(truth_presence & !pred_presence)
  TN <- sum(!truth_presence & !pred_presence)
  
  sensitivity <- TP / (TP + FN)
  precision <- TP / (TP + FP)
  specificity <- TN / (TN + FP)
  f1_score <- 2 * precision * sensitivity / (precision + sensitivity)
  
  # Jensen-Shannon divergence
  mat <- rbind(truth_vec, pred_vec)
  jsd <- distance(mat, method = "jensen-shannon")
  
  # Корреляции
  pearson_cor <- cor(truth_vec, pred_vec, method = "pearson")
  spearman_cor <- cor(truth_vec, pred_vec, method = "spearman")
  
  tibble(
    sensitivity = sensitivity,
    precision = precision,
    specificity = specificity,
    f1_score = f1_score,
    jensen_shannon = jsd,
    pearson = pearson_cor,
    spearman = spearman_cor
  )
}
metrics <- calculate_all_metrics(truth_abundance, kraken_abundance)
g <- metrics %>%
  pivot_longer(cols = everything(), names_to = "metric", values_to = "value") %>%
  ggplot(aes(x = metric, y = value, fill = metric)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = round(value, 2)), vjust = -0.3) +
  labs(title = "Метрики классификации", x = NULL, y = "Значение") +
  theme_bw()
ggsave(c(samlpe_name, 'metrics.png'), g)
write_csv2(metrics, c(samlpe_name,'_metrics.csv'))
write_csv2(truth_abundance, c(samlpe_name,'_truth_abundance.csv'))
write_csv2(kraken_abundance, c(samlpe_name,'_kraken_abundance.csv'))
