# Code to generate Figure 4 of the Jokura et al 2024 Ctenophore apical organ connectome paper

# source packages and functions ------------------------------------------------
source("analysis/scripts/packages_and_functions.R")


# read data---------------------------------------------------------------------

df_time_dif <- read_csv("analysis/data/arrest_rebeat_time_differences.csv")


# plot arrest re-beat graph ----------------------------------------------------

plot_arrest_rebeat <- 
  ggplot(df_time_dif) +
  aes(x = response, y = time, fill = response, color = response) +  
  geom_boxplot(aes(color = plane), size = 0.5, outlier.shape = NA) +
  geom_beeswarm(aes(color = time),
                size = 2,
                cex = 3.2,
                alpha = 0.3,
                color = "gray5") +
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
    values = c("re-beat" = "gray", 
               "arrest" = "white")) +
  scale_color_manual(
    values = c("sagittal" = "#0097b7", "tentacular" = "#f39500")) +
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


# assemble figure -------------------------------------------------------------

panel_ms <- ggdraw() + draw_image(readPNG("manuscript/pictures/tilt_microscope.png"))
panel_balancer <- ggdraw() + draw_image(readPNG("manuscript/pictures/balancer_closeup.png")) +
  draw_label("sagittal", x = 0.25, y = 0.95, size = 9, hjust = 0.5) +
  draw_label("tentacular", x = 0.75, y = 0.95, size = 9, hjust = 0.5)

panel_kymograph <- ggdraw() + draw_image(readPNG("manuscript/pictures/balacer_arrest_rebeat_kymograph.png"))
panel_graph <- ggdraw() + draw_image(readPNG("manuscript/pictures/arrest_rebeat_graph.png"))

panel_cbf_bar_sag <- ggdraw() + draw_image(readPNG("analysis/data/balancer_CBF_Pearson_correlation_analysis/output/CBF_barplot/CBF_barplot_23_08_04_WT1_030_Rec_at_100.73fps_dark.csv.png")) +
  draw_label("sagittal", x = 0.5, y = 0.95, size = 9, hjust = 0.5)
panel_cbf_bar_tent <- ggdraw() + draw_image(readPNG("analysis/data/balancer_CBF_Pearson_correlation_analysis/output/CBF_barplot/CBF_barplot_23_08_01_WT1_007_Rec_at_100.20fps_dark.csv.png")) +
  draw_label("tentacular", x = 0.5, y = 0.95, size = 9, hjust = 0.5)
panel_cor_graph <- ggdraw() + draw_image(readPNG("analysis/data/balancer_CBF_Pearson_correlation_analysis/output/pearson_boxplot/boxplot_correlation_angle_0-20_win_20.png"))

panel_comparison <- ggdraw() + draw_image(readPNG("manuscript/pictures/map_comparison.png"))

layout <- "
ABB
###
CCD
###
EFG
"

Figure4 <- panel_ms + panel_balancer + 
  panel_kymograph + panel_graph + 
  panel_cbf_bar_sag + panel_cbf_bar_tent + panel_cor_graph +
  plot_layout(design = layout,
              heights = c(1, 0.1, 1, 0.1, 1),
              widths = c(1, 1, 1)) + 
  patchwork::plot_annotation(tag_levels = "A") &  
  ggplot2::theme(plot.tag = element_text(size = 12, 
                                         face='plain', color='black'))

ggsave("manuscript/figures/Figure4.png", limitsize = FALSE, 
       units = c("px"), Figure4, width = 2000, height = 2000, bg='white')  


ggsave("manuscript/figures/Figure4.pdf", limitsize = FALSE, 
       units = c("px"), Figure4, width = 2000, height = 2000) 


