# Code to generate Figure 3 of the Jokura et al 2024 Ctenophore apical organ connectome paper

# source packages and functions ------------------------------------------------
source("analysis/scripts/packages_and_functions.R")


# assemble figure -------------------------------------------------------------


panel_AO_matrix <- ggdraw() + draw_image(readPNG("manuscript/pictures/AO_neurons_matrix.png"))
panel_bal_3d <- ggdraw() + draw_image(readPNG("manuscript/pictures/balancer.png"))
panel_bal_syn <- ggdraw() + draw_image(readPNG("manuscript/pictures/SSN_prepost_synapse_balancer.png"))
panel_grav_nc_graph <- ggdraw() + draw_image(readPNG("manuscript/pictures/gravity_neural_circuit_graph.png"))

layout <- "
AA#BBBBB
########
CCCCC#DD
"

Figure3 <- panel_AO_matrix + panel_bal_3d + 
  panel_bal_syn + panel_grav_nc_graph +
  plot_layout(design = layout,
              heights = c(1, 0.1, 1),
              widths = c(1, 1, 0.1, 1, 1, 0.1, 1, 1)) + 
  plot_annotation(tag_levels = "A") + 
  theme(plot.tag = element_text(size = 12, face='plain', color='black'))


ggsave("manuscript/figures/Figure3.png", limitsize = FALSE, 
       units = c("px"), Figure3, width = 2000, height = 1000, bg='white')  


ggsave("manuscript/figures/Figure3.pdf", limitsize = FALSE, 
       units = c("px"), Figure3, width = 2000, height = 1000) 

