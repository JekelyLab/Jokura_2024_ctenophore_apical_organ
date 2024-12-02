# Code to generate Figure 1 of the Jokura et al 2024 Ctenophore apical organ connectome paper

# source packages and functions -----------------
source("analysis/scripts/packages_and_functions.R")


# assemble figure -------------------------------------------------------------

library(jpeg)
library(magick)

panel_larva_pic <- ggdraw() + draw_image(image_read("manuscript/pictures/Mnemiopsis_larva_5dpf.png"))
panel_AO_pic <- ggdraw() + draw_image(image_read("manuscript/pictures/AO_mag_pics.png"))
panel_serial_sectioning <- ggdraw() + draw_image(readPNG("manuscript/pictures/serial sectioning.png"))
panel_catmaid_overview <- ggdraw() + draw_image(readPNG("manuscript/pictures/overview.png"))
panel_3d_all_cells <- ggdraw() + draw_image(readPNG("manuscript/pictures/all_cells_3_views.png"))

layout <- "
AABBBBB#CCC
###########
DDEEEEEEEEE
"

Figure1 <- panel_larva_pic + panel_AO_pic + 
  panel_serial_sectioning + panel_catmaid_overview + panel_3d_all_cells +
  plot_layout(design = layout,heights = c(1,0.1,1.4),widths = c(1,1.5,1,0.2,1,1,1,0.15,1,1,1)) + 
  plot_annotation(tag_levels = "A") + 
  theme(plot.tag = element_text(size = 12, face='plain', color='black'))


ggsave("manuscript/figures/Figure1.png", limitsize = FALSE, 
       units = c("px"), Figure1, width = 2600, height = 1000, bg='white')  


ggsave("manuscript/figures/Figure1.pdf", limitsize = FALSE, 
       units = c("px"), Figure1, width = 2600, height = 1000) 

