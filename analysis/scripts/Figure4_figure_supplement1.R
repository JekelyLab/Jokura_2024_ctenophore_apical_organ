# Code to generate Figure 4 Supplement 1 of the Jokura et al 2024 Ctenophore apical organ connectome paper

# source packages and functions ------------------------------------------------
source("analysis/scripts/packages_and_functions.R")


# CBF barplot sagittal(S) vs tentacular(T) -------------------------------------------

# Lists files to be processed and Plane information
files_info <- list(
  list(
    path = "analysis/data/balancer_CBF_Pearson_correlation_analysis/csv/23_06_29_WT10_027_Rec_at_100.75fps_light.csv",
    plane = "S"
  ),
  list(
    path = "analysis/data/balancer_CBF_Pearson_correlation_analysis/csv/23_07_25_WT2_029_Rec_at_100.07fps_dark.csv",
    plane = "S"
  ),
  list(
    path = "analysis/data/balancer_CBF_Pearson_correlation_analysis/csv/23_08_09_WT1_001_Rec_at_100.24fps_dark.csv",
    plane = "S"
  ),
  list(
    path = "analysis/data/balancer_CBF_Pearson_correlation_analysis/csv/23_08_01_WT1_006_Rec_at_100.20fps_dark.csv",
    plane = "T"
  ),
  list(
    path = "analysis/data/balancer_CBF_Pearson_correlation_analysis/csv/23_08_01_WT1_008_Rec_at_100.20fps_light.csv",
    plane = "T"
  ),
  list(
    path = "analysis/data/balancer_CBF_Pearson_correlation_analysis/csv/23_08_01_WT1_009_Rec_at_100.20fps_light.csv",
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
    filename = file.path(output_dir, paste0("CBF_barplot_sup_", file_name, ".png")),
    plot = cbf,
    width = 4,
    height = 3,
    dpi = 300,
    bg = "white"
  )
}

# save to source data (CBF barplot) -------------------------------------------------------

df_s1 <- read_csv("analysis/data/balancer_CBF_Pearson_correlation_analysis/csv/23_06_29_WT10_027_Rec_at_100.75fps_light.csv", show_col_types = FALSE) %>%
  filter(!is.na(left) & !is.na(right)) %>%
  filter(left != 0 & right != 0) %>%
  filter(is.finite(left) & is.finite(right)) %>%
  mutate(frame = row_number())

df_s1 %>%
  write_csv("manuscript/source_data/Figure4_Supplement1_source_data1.csv")

df_s2 <- read_csv("analysis/data/balancer_CBF_Pearson_correlation_analysis/csv/23_07_25_WT2_029_Rec_at_100.07fps_dark.csv", show_col_types = FALSE) %>%
  filter(!is.na(left) & !is.na(right)) %>%
  filter(left != 0 & right != 0) %>%
  filter(is.finite(left) & is.finite(right)) %>%
  mutate(frame = row_number())

df_s2 %>%
  write_csv("manuscript/source_data/Figure4_Supplement1_source_data2.csv")

df_s3 <- read_csv("analysis/data/balancer_CBF_Pearson_correlation_analysis/csv/23_08_09_WT1_001_Rec_at_100.24fps_dark.csv", show_col_types = FALSE) %>%
  filter(!is.na(left) & !is.na(right)) %>%
  filter(left != 0 & right != 0) %>%
  filter(is.finite(left) & is.finite(right)) %>%
  mutate(frame = row_number())

df_s3 %>%
  write_csv("manuscript/source_data/Figure4_Supplement1_source_data3.csv")

df_t1 <- read_csv("analysis/data/balancer_CBF_Pearson_correlation_analysis/csv/23_08_01_WT1_006_Rec_at_100.20fps_dark.csv", show_col_types = FALSE) %>%
  filter(!is.na(left) & !is.na(right)) %>%
  filter(left != 0 & right != 0) %>%
  filter(is.finite(left) & is.finite(right)) %>%
  mutate(frame = row_number())

df_t1 %>%
  write_csv("manuscript/source_data/Figure4_Supplement1_source_data4.csv")

df_t2 <- read_csv("analysis/data/balancer_CBF_Pearson_correlation_analysis/csv/23_08_01_WT1_008_Rec_at_100.20fps_light.csv", show_col_types = FALSE) %>%
  filter(!is.na(left) & !is.na(right)) %>%
  filter(left != 0 & right != 0) %>%
  filter(is.finite(left) & is.finite(right)) %>%
  mutate(frame = row_number())

df_t2 %>%
  write_csv("manuscript/source_data/Figure4_Supplement1_source_data5.csv")

df_t3 <- read_csv("analysis/data/balancer_CBF_Pearson_correlation_analysis/csv/23_08_01_WT1_009_Rec_at_100.20fps_light.csv", show_col_types = FALSE) %>%
  filter(!is.na(left) & !is.na(right)) %>%
  filter(left != 0 & right != 0) %>%
  filter(is.finite(left) & is.finite(right)) %>%
  mutate(frame = row_number())

df_t3 %>%
  write_csv("manuscript/source_data/Figure4_Supplement1_source_data6.csv")


# assemble figure -------------------------------------------------------------

sup_panel_cbf_sag1 <- ggdraw() + draw_image(readPNG("manuscript/pictures/CBF_barplot_sup_23_06_29_WT10_027_Rec_at_100.75fps_light.png"))
sup_panel_cbf_sag2 <- ggdraw() + draw_image(readPNG("manuscript/pictures/CBF_barplot_sup_23_07_25_WT2_029_Rec_at_100.07fps_dark.png"))
sup_panel_cbf_sag3 <- ggdraw() + draw_image(readPNG("manuscript/pictures/CBF_barplot_sup_23_08_09_WT1_001_Rec_at_100.24fps_dark.png"))

sup_panel_cbf_tent1 <- ggdraw() + draw_image(readPNG("manuscript/pictures/CBF_barplot_sup_23_08_01_WT1_006_Rec_at_100.20fps_dark.png"))
sup_panel_cbf_tent2 <- ggdraw() + draw_image(readPNG("manuscript/pictures/CBF_barplot_sup_23_08_01_WT1_008_Rec_at_100.20fps_light.png"))
sup_panel_cbf_tent3 <- ggdraw() + draw_image(readPNG("manuscript/pictures/CBF_barplot_sup_23_08_01_WT1_009_Rec_at_100.20fps_light.png"))




layout <- "
A#B#C
#####
D#E#F
"

Fig4_Sup1 <- sup_panel_cbf_sag1 + sup_panel_cbf_sag2 + sup_panel_cbf_sag3 +
  sup_panel_cbf_tent1 + sup_panel_cbf_tent2 + sup_panel_cbf_tent3 +
  plot_layout(design = layout,
              heights = c(1, 0.05, 1),
              widths = c(1, 0.05, 1, 0.05, 1)) + 
  patchwork::plot_annotation(tag_levels = "A") &  
  ggplot2::theme(plot.tag = element_text(size = 12, 
                                         face='plain', color='black'))

ggsave("manuscript/figures/Figure4_Supplement1.png", limitsize = FALSE, 
       units = c("px"), Fig4_Sup1, width = 2400, height = 1200, bg='white')  


ggsave("manuscript/figures/Figure4_Supplement1.pdf", limitsize = FALSE, 
       units = c("px"), Fig4_Sup1, width = 2400, height = 1200) 

