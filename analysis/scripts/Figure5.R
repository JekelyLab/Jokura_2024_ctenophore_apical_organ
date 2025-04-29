# Code to generate Figure 5 of the Jokura et al 2024 Ctenophore apical organ connectome paper

# source packages and functions ------------------------------------------------
source("analysis/scripts/packages_and_functions.R")

# assemble figure -------------------------------------------------------------

panel_comparison <- ggdraw() + draw_image(readPNG("manuscript/pictures/map_comparison.png"))

layout <- "
A
"

Figure5 <- panel_comparison +
  plot_layout(design = layout,
              heights = c(1),
              widths = c(1)) + 
  patchwork::plot_annotation(tag_levels = "A") &  
  ggplot2::theme(plot.tag = element_text(size = 12, 
                                         face='plain', color='black'))

ggsave("manuscript/figures/Figure5.png", limitsize = FALSE, 
       units = c("px"), Figure5, width = 2100, height = 1100, bg='white')  


ggsave("manuscript/figures/Figure5.pdf", limitsize = FALSE, 
       units = c("px"), Figure5, width = 2100, height = 1100) 


