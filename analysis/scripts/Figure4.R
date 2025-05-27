# Code to generate Figure 4 of the Jokura et al 2024 Ctenophore apical organ connectome paper

# source packages and functions ------------------------------------------------
source("analysis/scripts/packages_and_functions.R")


# read arrest re-beat data------------------------------------------------------

df_time_dif <- read_csv("analysis/data/arrest_rebeat_time_differences.csv")

# plot arrest re-beat graph ----------------------------------------------------

plot_arrest_rebeat <- 
  ggplot(df_time_dif) +
  aes(x = response, y = time, fill = response) +  
  geom_boxplot(aes(color = plane), 
               size = 0.5, 
               outlier.shape = NA) +
  geom_beeswarm(aes(color = plane),    
                size = 1.75,
                cex = 3.2,
                alpha = 0.5,
                show.legend = FALSE) +
  labs(y = "Time difference (sec)") +
  xlab("") +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.3, 0.75),
    strip.text = element_text(size = 11),  
    strip.placement = "outside",  
    strip.background = element_blank(),  
    axis.text.x = element_blank(),  
    axis.ticks.x = element_blank()  
  ) +
  scale_fill_manual(
    values = c("re-beat" = "gray70", 
               "arrest" = "white")) +
  scale_color_manual(
    values = c("sagittal" = "#0097b7", 
               "tentacular" = "#f39500"),
    guide = "none") +
  facet_grid(~ plane, scales = "free_x", space = "free_x", switch = "x")

plot_arrest_rebeat

ggsave(
  filename = "manuscript/pictures/arrest_rebeat_graph.png",
  plot = plot_arrest_rebeat,
  width = 900,
  height = 1050,
  units = "px",
  dpi = 300, 
  bg = "white"
)

# save to source data (plot arrest re-beat graph) -------------------------------------------------------

df_time_dif %>%
  write_csv("manuscript/source_data/Figure4_source_data1.csv")


# CBF barplot sagittal(S) vs tentacular(T) -------------------------------------------

# Lists files to be processed and Plane information
files_info <- list(
  list(
    path = "analysis/data/balancer_CBF_Pearson_correlation_analysis/csv/23_08_04_WT1_030_Rec_at_100.73fps_dark_f.csv",
    plane = "S"
  ),
  list(
    path = "analysis/data/balancer_CBF_Pearson_correlation_analysis/csv/23_08_01_WT1_007_Rec_at_100.20fps_dark.csv",
    plane = "T"
  )
)

# destination folder
output_dir <- "manuscript/pictures"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Color Setting
colors <- list(S = c("#0e6b76", "#77c2cc"), T = c("#a55e00", "#f5be6b"))

# Loop processing for each file
for (info in files_info) {
  
  file_path <- info$path
  plane <- info$plane
  file_name <- tools::file_path_sans_ext(basename(file_path))
  
  # File loading and preprocessing
  df <- read_csv(file_path, show_col_types = FALSE) %>%
    filter(!is.na(left) & !is.na(right)) %>%
    filter(left != 0 & right != 0) %>%
    filter(is.finite(left) & is.finite(right)) %>%
    mutate(frame = row_number())
  
  # color coding
  c_left <- colors[[plane]][1]
  c_right <- colors[[plane]][2]
  
  # Barplot
  cbf <- ggplot(df, aes(x = frame)) +
    geom_bar(aes(y = left), stat = "identity", fill = c_left, alpha = 0.7, na.rm = TRUE) +
    geom_bar(aes(y = -right), stat = "identity", fill = c_right, alpha = 0.7, na.rm = TRUE) +
    labs(x = "Time (sec)", y = "CBF (Hz)") +
    theme_minimal() +
    theme(
      plot.title = element_blank(),
      text = element_text(size = 14),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12)
    )
  
  # save
  ggsave(
    filename = file.path(output_dir, paste0("CBF_barplot_", file_name, ".png")),
    plot = cbf,
    width = 4,
    height = 3,
    dpi = 300,
    bg = "white"
  )
}

panel_cbf_bar_sag <- ggdraw() + draw_image(readPNG("manuscript/pictures/CBF_barplot_23_08_04_WT1_030_Rec_at_100.73fps_dark_f.png"))
panel_cbf_bar_tent <- ggdraw() + draw_image(readPNG("manuscript/pictures/CBF_barplot_23_08_01_WT1_007_Rec_at_100.20fps_dark.png"))

panel_cbf_bars_combined <- panel_cbf_bar_sag + panel_cbf_bar_tent +
  plot_layout(ncol = 2, widths = c(1, 1))

ggsave("manuscript/pictures/panel_cbf_bars_combined.png",
       plot = panel_cbf_bars_combined,
       width = 8, height = 3, dpi = 300, bg = "white")

# save to source data (CBF barplot) -------------------------------------------------------

df_s <- read_csv("analysis/data/balancer_CBF_Pearson_correlation_analysis/csv/23_08_04_WT1_030_Rec_at_100.73fps_dark_f.csv", show_col_types = FALSE) %>%
  filter(!is.na(left) & !is.na(right)) %>%
  filter(left != 0 & right != 0) %>%
  filter(is.finite(left) & is.finite(right)) %>%
  mutate(frame = row_number())

df_s %>%
  write_csv("manuscript/source_data/Figure4_source_data2.csv")


df_t <- read_csv("analysis/data/balancer_CBF_Pearson_correlation_analysis/csv/23_08_01_WT1_007_Rec_at_100.20fps_dark.csv", show_col_types = FALSE) %>%
  filter(!is.na(left) & !is.na(right)) %>%
  filter(left != 0 & right != 0) %>%
  filter(is.finite(left) & is.finite(right)) %>%
  mutate(frame = row_number())

df_t %>%
  write_csv("manuscript/source_data/Figure4_source_data3.csv")

# boxplot of pearson comparison S vs T (windows = 20)---------------------------

# file and output path
input_csv <- "analysis/data/balancer_CBF_Pearson_correlation_analysis/output/csv/correlation_results.csv"

# data loading
df <- read_csv(input_csv, show_col_types = FALSE)

# Get S and T values
s_vals <- df %>% filter(plane == "S") %>% pull(mean_correlation)
t_vals <- df %>% filter(plane == "T") %>% pull(mean_correlation)

# Automatic t-test or Wilcoxon test
shapiro_s <- shapiro.test(s_vals)
shapiro_t <- shapiro.test(t_vals)
if (shapiro_s$p.value > 0.05 & shapiro_t$p.value > 0.05) {
  test <- t.test(s_vals, t_vals)
  test_name <- "t-test, "
} else {
  test <- wilcox.test(s_vals, t_vals)
  test_name <- "Wilcoxon test, "
}

# Boxplot 
pearson <- ggplot(df, aes(x = plane, y = mean_correlation, color = plane)) +
  geom_boxplot(fill = "lightgray", alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.6) +
  scale_color_manual(values = c("S" = "#148daa", "T" = "#e89014")) +
  scale_x_discrete(labels = c("S" = "sagittal", "T" = "tentacular")) +
  labs(
    x = NULL,
    y = "Mean Rolling Pearson Correlation",
    title = paste(test_name, "**p =", round(test$p.value, 3))
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),
    axis.text = element_text(size = 14),
    legend.position = "none"
  ) +
  ylim(-1, 1)

pearson

# Save
ggsave(filename = "manuscript/pictures/boxplot_correlation_angle_0-20_win_20.png", 
       plot = pearson,
       width = 1050,
       height = 1150,
       units = "px",
       dpi = 300, 
       bg = "white")

# save to source data (boxplot of pearson comparison) -------------------------------------------------------

df %>%
  write_csv("manuscript/source_data/Figure4_source_data4.csv")


# assemble figure -------------------------------------------------------------

panel_ms <- ggdraw() + draw_image(readPNG("manuscript/pictures/tilt_microscope.png"))
panel_balancer <- ggdraw() + draw_image(readPNG("manuscript/pictures/balancer_closeup.png")) +
  draw_label("sagittal plane", x = 0.25, y = 0.95, size = 10, hjust = 0.5) +
  draw_label("tentacular plane", x = 0.75, y = 0.95, size = 10, hjust = 0.5)

panel_kymograph <- ggdraw() + draw_image(readPNG("manuscript/pictures/balacer_arrest_rebeat_kymograph.png"))
panel_graph <- ggdraw() + draw_image(readPNG("manuscript/pictures/arrest_rebeat_graph.png"))

panel_cbf_bars <- ggdraw() + draw_image(readPNG("manuscript/pictures/panel_cbf_bars_combined.png")) +
  draw_label("sagittal plane", x = 0.25, y = 0.975, size = 10, hjust = 0.5) +
  draw_label("tentacular plane", x = 0.75, y = 0.975, size = 10, hjust = 0.5) 

panel_cor_graph <- ggdraw() + draw_image(readPNG("manuscript/pictures/boxplot_correlation_angle_0-20_win_20.png"))

layout <- "
A#BBB
#####
CCC#D
#####
EEE#F
"

Figure4 <- panel_ms + panel_balancer + 
  panel_kymograph + panel_graph + 
  panel_cbf_bars + panel_cor_graph +
  plot_layout(design = layout,
              heights = c(1, 0.05, 1.2, 0.2, 1),
              widths = c(1, 0.05, 1, 0.05, 1)) + 
  patchwork::plot_annotation(tag_levels = "A") &  
  ggplot2::theme(plot.tag = element_text(size = 12, 
                                         face='plain', color='black'))

ggsave("manuscript/figures/Figure4.png", limitsize = FALSE, 
       units = c("px"), Figure4, width = 2100, height = 2200, bg='white')  


ggsave("manuscript/figures/Figure4.pdf", limitsize = FALSE, 
       units = c("px"), Figure4, width = 2100, height = 2200) 


