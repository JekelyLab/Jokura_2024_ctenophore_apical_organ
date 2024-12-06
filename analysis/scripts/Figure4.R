# Code to generate Figure 4 of the Jokura et al 2024 Ctenophore apical organ connectome paper

# source packages and functions ------------------------------------------------
source("analysis/scripts/packages_and_functions.R")


# assemble figure -------------------------------------------------------------


panel_bri_3d <- ggdraw() + draw_image(readPNG("manuscript/pictures/bridge_text.png"))
panel_bri_mito <- ggdraw() + draw_image(readPNG("manuscript/pictures/mito_pos_bridge_text.png"))
panel_bri_syn <- ggdraw() + draw_image(readPNG("manuscript/pictures/bridge_syn_text.png"))
panel_grav_matrix <- ggdraw() + draw_image(readPNG("manuscript/pictures/mech_girdle_chaeMech_syn_matrix.png"))
panel_grav_graph_bridge <- ggdraw() + draw_image(readPNG("manuscript/pictures/gravity_neural_circuit_graph_bridge.png"))

layout <- "
AAAABBB
#######
CCCDDEE
"

Figure4 <- panel_bri_3d + panel_bri_mito + 
  panel_bri_syn + panel_grav_matrix + panel_grav_graph_bridge +
  plot_layout(design = layout,
              heights = c(1, 0.1, 1.2),
              widths = c(1, 1, 1, 1, 1, 1, 0.5)) + 
  plot_annotation(tag_levels = "A") + 
  theme(plot.tag = element_text(size = 12, face='plain', color='black'))


ggsave("manuscript/figures/Figure4.png", limitsize = FALSE, 
       units = c("px"), Figure4, width = 2400, height = 1000, bg='white')  


ggsave("manuscript/figures/Figure4.pdf", limitsize = FALSE, 
       units = c("px"), Figure4, width = 2400, height = 1000) 

