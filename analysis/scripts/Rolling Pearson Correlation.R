rm(list = ls())

library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(zoo)

# Paths
balancer_info_path <- "analysis/balancer_CBF_Pearson_correlation_analysis/csv/balancer_info.csv"
csv_folder <- "analysis/balancer_CBF_Pearson_correlation_analysis/csv"
output_barplot <- "analysis/balancer_CBF_Pearson_correlation_analysis/output/CBF_barplot"
output_corrplot <- "analysis/balancer_CBF_Pearson_correlation_analysis/output/pearson_rolling"
output_boxplot <- "analysis/balancer_CBF_Pearson_correlation_analysis/output/pearson_boxplot"

dir.create(output_barplot, recursive = TRUE, showWarnings = FALSE)
dir.create(output_corrplot, recursive = TRUE, showWarnings = FALSE)
dir.create(output_boxplot, recursive = TRUE, showWarnings = FALSE)

# Read info
balancer_info <- read_csv(balancer_info_path)

# Filter by abs_angle
filtered_files <- balancer_info %>%
  filter(abs_angle >= 0 & abs_angle < 20) %>%
  select(file_name, plane)

# Barplot
plot_bar <- function(df, file_name, plane) {
  df <- df %>%
    replace_na(list(left = 0, right = 0)) %>%
    filter(is.finite(left) & is.finite(right)) %>%
    mutate(frame = row_number())
  
  colors <- list(S = c("deepskyblue", "royalblue"), T = c("orange", "orangered"))
  c_left <- colors[[plane]][1]; c_right <- colors[[plane]][2]
  
  p <- ggplot(df, aes(x = frame)) +
    geom_bar(aes(y = left), stat = "identity", fill = c_left, alpha = 0.7, na.rm = TRUE) +
    geom_bar(aes(y = -right), stat = "identity", fill = c_right, alpha = 0.7, na.rm = TRUE) +
    labs(x = "Time (sec)", y = "CBF (Hz)") +
    ylim(-25, 25) + xlim(0, 120) +
    theme_minimal() +
    theme(plot.title = element_blank())
  
  ggsave(file.path(output_barplot, paste0("CBF_barplot_", file_name, ".png")),
         plot = p, width = 8, height = 4.5, dpi = 300)
}

# Rolling correlation
calc_rolling_corr <- function(df, window) {
  rollapply(1:(nrow(df) - window + 1), 1, function(i) {
    x <- df$left[i:(i + window - 1)]
    y <- df$right[i:(i + window - 1)]
    clean <- na.omit(data.frame(x, y))
    if (nrow(clean) < 2 || sd(clean$x) == 0 || sd(clean$y) == 0) return(NA)
    cor(clean$x, clean$y, use = "complete.obs")
  }, by.column = FALSE)
}

# Plot rolling correlation
plot_rolling <- function(corr, file_name, plane, window) {
  plot_data <- data.frame(time_index = 1:length(corr), rolling_corr = corr)
  color <- ifelse(plane == "S", "#0F52BA", "#F28500")
  
  p <- ggplot(plot_data, aes(x = time_index, y = rolling_corr)) +
    geom_line(color = color, size = 1) +
    geom_point(color = color, size = 2) +
    labs(x = "Time Index", y = "Rolling Pearson Correlation") +
    theme_minimal() +
    theme(plot.title = element_blank()) +
    ylim(-1, 1)
  
  sub_dir <- file.path(output_corrplot, paste0("window_", window))
  if (!dir.exists(sub_dir)) dir.create(sub_dir, recursive = TRUE)
  
  ggsave(file.path(sub_dir, paste0("rolling_corr_", file_name, "_win_", window, ".png")),
         plot = p, width = 8, height = 5, dpi = 300)
}

# Boxplot
plot_boxplot <- function(results, window) {
  s_vals <- results %>% filter(plane == "S") %>% pull(mean_correlation)
  t_vals <- results %>% filter(plane == "T") %>% pull(mean_correlation)
  
  shapiro_s <- shapiro.test(s_vals)
  shapiro_t <- shapiro.test(t_vals)
  if (shapiro_s$p.value > 0.05 & shapiro_t$p.value > 0.05) {
    test <- t.test(s_vals, t_vals); test_name <- "t-test"
  } else {
    test <- wilcox.test(s_vals, t_vals); test_name <- "Wilcoxon test"
  }
  
  p <- ggplot(results, aes(x = plane, y = mean_correlation, color = plane)) +
    geom_boxplot(fill = "lightgray", alpha = 0.5, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 3, alpha = 0.6) +
    scale_color_manual(values = c("S" = "#0F52BA", "T" = "#F28500")) +
    labs(x = "Plane", y = paste("Mean Rolling Pearson Correlation (Window =", window, ")"),
         title = paste("Angle 0–20,", test_name, "p =", round(test$p.value, 3))) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylim(-1, 1)
  
  ggsave(file.path(output_boxplot, paste0("boxplot_correlation_angle_0-20_win_", window, ".png")),
         plot = p, width = 8, height = 5, dpi = 300)
  print(test)
}

# Run
window_sizes <- c(10, 20, 30)
for (window in window_sizes) {
  results <- data.frame(file_name = character(), plane = character(), mean_correlation = numeric())
  for (i in 1:nrow(filtered_files)) {
    file_name <- filtered_files$file_name[i]
    plane_info <- filtered_files$plane[i]
    file_path <- file.path(csv_folder, file_name)
    if (!file.exists(file_path)) next
    
    df <- read_csv(file_path)
    if (nrow(df) < window) next
    
    plot_bar(df, file_name, plane_info)
    
    corr <- calc_rolling_corr(df, window)
    plot_rolling(corr, file_name, plane_info, window)
    
    mean_corr <- mean(corr, na.rm = TRUE)
    results <- rbind(results, data.frame(file_name = file_name, plane = plane_info, mean_correlation = mean_corr))
  }
  plot_boxplot(results, window)
}
