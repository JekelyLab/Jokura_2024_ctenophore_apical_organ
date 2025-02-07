# Code to generate Figure 6 of the Jokura et al 2024 Ctenophore apical organ connectome paper

# source packages and functions ------------------------------------------------
source("analysis/scripts/packages_and_functions.R")


# assemble figure -------------------------------------------------------------

panel_sum <- ggdraw() + draw_image(readPNG("manuscript/pictures/discussion.png"))

layout <- "
A
"



Figure6 <- panel_sum +
  plot_layout(design = layout,
              heights = c(1),
              widths = c(1)) + 
  plot_annotation(tag_levels = "A") + 
  theme(plot.tag = element_text(size = 12, face='plain', color='black'))


ggsave("manuscript/figures/Figure6.png", limitsize = FALSE, 
       units = c("px"), Figure6, width = 2500, height = 1050, bg='white')  


ggsave("manuscript/figures/Figure6.pdf", limitsize = FALSE, 
       units = c("px"), Figure6, width = 2500, height = 1050) 

