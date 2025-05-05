# Code to generate Figure 2 Supplement 1 of the Jokura et al 2024 Ctenophore apical organ connectome paper

# source packages and functions ------------------------------------------------
source("analysis/scripts/packages_and_functions.R")


# load Subepithelial nerve net -------------------------------------------------

SNN <- read_smooth_neuron("celltype:SNN")

# 3d plot Subepithelial nerve net -----------------------------------------------

close3d()
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  
#define the size of the rgl window, the view and zoom
par3d(windowRect = c(0, 0, 1200, 350))

#plot aboral view
plot_multinucleated_cell(SNN, lwd = 1.25, alpha = 0.8, col = Okabe_Ito[5])
plot3d(outline, add = TRUE, alpha = 0.025, col = Okabe_Ito[8])

#aboral view
aboral()
par3d(zoom = 0.7)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Sagittal plane
plot_multinucleated_cell(SNN, lwd = 1.25, alpha = 0.8, col = Okabe_Ito[5])
plot3d(outline, add = TRUE, alpha = 0.025, col = Okabe_Ito[8])

#lateral view of Sagittal plane
sagittal()
par3d(zoom = 0.7)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Tentacular plane
plot_multinucleated_cell(SNN, lwd = 1.25, alpha = 0.8, col = Okabe_Ito[5])
plot3d(outline, add = TRUE, alpha = 0.025, col = Okabe_Ito[8])

#lateral view of Tentacular plane
tentacular()
par3d(zoom = 0.7)

#move to next panel in rgl window
next3d(clear=F)


#make a snapshot to the working directory
rgl.snapshot("manuscript/pictures/3d_plot/plot_SNN.png")


close3d()


# assemble figure -------------------------------------------------------------

panel_SNN_3d <- ggdraw() + draw_image(readPNG("manuscript/pictures/3d_plot/plot_SNN.png"))
panel_SNN_EM <- ggdraw() + draw_image(readPNG("manuscript/pictures/SNN_EM.png"))

panel_CAT_pic <- ggdraw() + draw_image(readPNG("manuscript/pictures/figure_sup2_1.png"))

panel_mono_vs_poly <- ggdraw() + draw_image(readPNG("manuscript/pictures/figure_sup2_2.png"))

panel_multi_nucleus <- ggdraw() + draw_image(readPNG("manuscript/pictures/multi_nucleus.png"))

layout <- "
ABC
DE#
"

Fig2_Sup1 <- panel_CAT_pic + panel_mono_vs_poly + panel_multi_nucleus +
  panel_SNN_3d + panel_SNN_EM +
  plot_layout(design = layout,
              heights = c(),
              widths = c()) + 
  patchwork::plot_annotation(tag_levels = "A") &  
  ggplot2::theme(plot.tag = element_text(size = 12, 
                                         face='plain', color='black'))

ggsave("manuscript/figures/Figure2_Supplement1.png", limitsize = FALSE, 
       units = c("px"), Fig2_Sup1, width = 2700, height = 1000, bg='white')  

ggsave("manuscript/figures/Figure2_Supplement1.pdf", limitsize = FALSE, 
      units = c("px"), Fig2_Sup1, width = 2700, height = 1000) 












