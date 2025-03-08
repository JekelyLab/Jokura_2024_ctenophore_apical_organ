# Code to generate Figure 4 of the Jokura et al 2024 Ctenophore apical organ connectome paper

# source packages and functions ------------------------------------------------
source("analysis/scripts/packages_and_functions.R")


# read data---------------------------------------------------------------------

df_time_dif <- read_csv("analysis/data/arrest_rebeat_time_differences.csv")


# plot arrest re-beat graph ----------------------------------------------------

plot_arrest_rebeat <- 
  ggplot(df_time_dif) +
  aes(x = interaction(response, plane), y = time, fill = response) +
  geom_boxplot(size = 0.5, outlier.shape = NA) +
  geom_beeswarm(aes(color = time),
                size = 2,
                cex = 3.2,
                alpha = 0.3,
                color = "gray5") +
  labs(y = "time difference (sec)") +
  xlab("") +
  theme_minimal() +
  theme(#axis.text.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.3, 0.75)
        ) +
  scale_fill_manual(
    values = c("re-beat" = "#28A8FF", 
               "arrest" = "#f56EBA")) +
  scale_x_discrete(labels = c("arrest.sagittal" = "sagittal",
                              "re-beat.sagittal" = "",
                              "arrest.tentacular" = "tentacular",
                              "re-beat.tentacular" = ""))

plot_arrest_rebeat

ggsave(
  filename = "manuscript/pictures/arrest_rebeat_graph.png",
  plot = plot_arrest_rebeat,
  width = 900,
  height = 1050,
  units = "px",
  dpi = 300
)


# assemble figure -------------------------------------------------------------

panel_ms <- ggdraw() + draw_image(readPNG("manuscript/pictures/tilt_microscope.png"))

panel_kymograph <- ggdraw() + draw_image(readPNG("manuscript/pictures/balacer_arrest_rebeat_kymograph.png"))

panel_graph <- ggdraw() + draw_image(readPNG("manuscript/pictures/arrest_rebeat_graph.png"))

panel_comparison <- ggdraw() + draw_image(readPNG("manuscript/pictures/map_comparison.png"))

layout <- "
ABB
###
CDD
"

Figure4 <- panel_ms + panel_kymograph + 
  panel_graph + panel_comparison +
  plot_layout(design = layout,
              heights = c(1, 0.1, 1),
              widths = c(1, 1, 1)) + 
  patchwork::plot_annotation(tag_levels = "A") &  
  ggplot2::theme(plot.tag = element_text(size = 12, 
                                         face='plain', color='black'))

ggsave("manuscript/figures/Figure4.png", limitsize = FALSE, 
       units = c("px"), Figure4, width = 2300, height = 1500, bg='white')  


ggsave("manuscript/figures/Figure4.pdf", limitsize = FALSE, 
       units = c("px"), Figure4, width = 2300, height = 1500) 


