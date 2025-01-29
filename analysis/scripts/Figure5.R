# Code to generate Figure 5 of the Jokura et al 2024 Ctenophore apical organ connectome paper

# source packages and functions ------------------------------------------------
source("analysis/scripts/packages_and_functions.R")


# assemble figure -------------------------------------------------------------

panel_ms <- ggdraw() + draw_image(readPNG("manuscript/pictures/tilt_microscope.png"))
panel_bal_s <- ggdraw() + draw_image(readPNG("manuscript/pictures/balancer_sagittal.png"))
panel_bal_t <- ggdraw() + draw_image(readPNG("manuscript/pictures/balancer_tentacular.png"))
panel_map_s <- ggdraw() + draw_image(readPNG("manuscript/pictures/map_sagittal.png"))
panel_map_t <- ggdraw() + draw_image(readPNG("manuscript/pictures/map_tentacular.png"))
panel_kymograph <- ggdraw() + draw_image(readPNG("manuscript/pictures/balacer_arrest_rebeat_kymograph.png"))
panel_graph <- ggdraw() + draw_image(readPNG("manuscript/pictures/balacer_arrest_rebeat_graph.png"))

layout <- "
AA#B#C
AA####
AA#D#E
######
FFFF#G
"



Figure5 <- panel_ms + panel_bal_s + panel_bal_t + panel_map_s + panel_map_t +
  panel_kymograph + panel_graph +
  plot_layout(design = layout,
              heights = c(1.5, 0.1, 1, 0.1, 3),
              widths = c(1, 1, 0.1, 1, 0.1, 1)) + 
  plot_annotation(tag_levels = "A") + 
  theme(plot.tag = element_text(size = 12, face='plain', color='black'))


ggsave("manuscript/figures/Figure5.png", limitsize = FALSE, 
       units = c("px"), Figure5, width = 2500, height = 1800, bg='white')  


ggsave("manuscript/figures/Figure5.pdf", limitsize = FALSE, 
       units = c("px"), Figure5, width = 2500, height = 1800) 

