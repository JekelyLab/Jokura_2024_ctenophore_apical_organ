# Code to generate Figure 4 of the Jokura et al 2024 Ctenophore apical organ connectome paper

# source packages and functions ------------------------------------------------
source("analysis/scripts/packages_and_functions.R")


# assemble figure -------------------------------------------------------------

panel_ms <- ggdraw() + draw_image(readPNG("manuscript/pictures/tilt_microscope.png"))

panel_map_s <- ggdraw() + draw_image(readPNG("manuscript/pictures/map_sagittal.png"))
panel_map_t <- ggdraw() + draw_image(readPNG("manuscript/pictures/map_tentacular.png"))
panel_kymograph <- ggdraw() + draw_image(readPNG("manuscript/pictures/balacer_arrest_rebeat_kymograph.png"))

panel_graph <- ggdraw() + draw_image(readPNG("manuscript/pictures/balacer_arrest_rebeat_graph.png"))

panel_comparison <- ggdraw() + draw_image(readPNG("manuscript/pictures/comparison.png"))

layout <- "
AAA#B#C
AAA#DDD
AAA#DDD
#######
EEFFFFF
"



Figure4 <- panel_ms + panel_map_s + panel_map_t +
  panel_kymograph + 
  panel_graph + panel_comparison +
  plot_layout(design = layout,
              heights = c(1, 1, 1, 0.1, 3),
              widths = c(1, 1, 1, 0.1, 1, 0.1, 1)) + 
  patchwork::plot_annotation(tag_levels = "A") &  
  ggplot2::theme(plot.tag = element_text(size = 12, 
                                         face='plain', color='black'))



ggsave("manuscript/figures/Figure4.png", limitsize = FALSE, 
       units = c("px"), Figure4, width = 2400, height = 1700, bg='white')  


ggsave("manuscript/figures/Figure4.pdf", limitsize = FALSE, 
       units = c("px"), Figure4, width = 2400, height = 1700) 


