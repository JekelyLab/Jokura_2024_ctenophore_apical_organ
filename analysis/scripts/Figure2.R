# Code to generate Figure 2 of the Jokura et al 2024 Ctenophore apical organ connectome paper

# source packages and functions ------------------------------------------------
source("analysis/scripts/packages_and_functions.R")


# assemble figure -------------------------------------------------------------


panel_synapse <- ggdraw() + draw_image(readPNG("manuscript/pictures/Figure_mito_syn_ves_syn.png"))
panel_SSN_Q1234 <- ggdraw() + draw_image(readPNG("manuscript/pictures/SSN_Q1234_panel.png"))
panel_SSN_Q12_Q34 <- ggdraw() + draw_image(readPNG("manuscript/pictures/SSN_Q1Q2_Q3Q4_panel.png"))
panel_mito <- ggdraw() + draw_image(readPNG("manuscript/pictures/mito_pos_SSN_anterior_trim.png"))

layout <- "
AA#BBB
######
DD#CCC
"

Figure2 <- panel_synapse + 
  panel_SSN_Q1234 + panel_SSN_Q12_Q34 +
  panel_mito +
  plot_layout(design = layout,
              heights = c(1, 0.1, 1),
              widths = c(1, 1, 0.1, 1, 1, 1)) + 
  plot_annotation(tag_levels = "A") + 
  theme(plot.tag = element_text(size = 12, face='plain', color='black'))


ggsave("manuscript/figures/Figure2.png", limitsize = FALSE, 
       units = c("px"), Figure2, width = 2400, height = 1000, bg='white')  


ggsave("manuscript/figures/Figure2.pdf", limitsize = FALSE, 
       units = c("px"), Figure2, width = 2400, height = 1000) 

