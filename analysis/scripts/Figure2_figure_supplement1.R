# Code to generate Figure 2 Supplement 1 of the Jokura et al 2024 Ctenophore apical organ connectome paper

# source packages and functions ------------------------------------------------
source("analysis/scripts/packages_and_functions.R")


# load Subepithelial nerve net -------------------------------------------------

SNN <- read_smooth_neuron("celltype:SNN")
ANN <- read_smooth_neuron("celltype:SSN")

# 3d plot Subepithelial nerve net -----------------------------------------------

close3d()
# 3d plotting of cells
nopen3d()
mfrow3d(1, 3)
# define the size of the rgl window, the view and zoom
par3d(windowRect = c(0, 0, 2400, 700))

# plot aboral view
plot_multinucleated_cell(SNN, lwd = 1.5, alpha = 1, col = Okabe_Ito[2])
plot3d(outline, add = TRUE, alpha = 0.025, col = Okabe_Ito[8])

# aboral view
aboral()
par3d(zoom = 0.7)

# move to next panel in rgl window
next3d(clear = F)



# plot lateral view of Sagittal plane
plot_multinucleated_cell(SNN, lwd = 1.5, alpha = 1, col = Okabe_Ito[2])
plot3d(outline, add = TRUE, alpha = 0.025, col = Okabe_Ito[8])

# lateral view of Sagittal plane
sagittal()
par3d(zoom = 0.7)

# move to next panel in rgl window
next3d(clear = F)



# plot lateral view of Tentacular plane
plot_multinucleated_cell(SNN, lwd = 1.5, alpha = 1, col = Okabe_Ito[2])
plot3d(outline, add = TRUE, alpha = 0.025, col = Okabe_Ito[8])

# lateral view of Tentacular plane
tentacular()
par3d(zoom = 0.7)

# move to next panel in rgl window
next3d(clear = F)


# make a snapshot to the working directory
rgl.snapshot("manuscript/pictures/3d_plot/plot_SNN.png")

# add ANN cells to plot
plot_multinucleated_cell(
  ANN, lwd = 2, alpha = 1, 
  col = c(Okabe_Ito[6], Okabe_Ito[7], Okabe_Ito[5])
  )
next3d(clear = F)
plot_multinucleated_cell(
  ANN, lwd = 2, alpha = 1, 
  col = c(Okabe_Ito[6], Okabe_Ito[7], Okabe_Ito[5])
)
next3d(clear = F)
plot_multinucleated_cell(
  ANN, lwd = 2, alpha = 1, 
  col = c(Okabe_Ito[6], Okabe_Ito[7], Okabe_Ito[5])
)

rgl.snapshot("manuscript/pictures/3d_plot/plot_SNN_ANN.png")

close3d()


# assemble figure -------------------------------------------------------------

panel_SNN_3d <- ggdraw() + draw_image(readPNG("manuscript/pictures/3d_plot/plot_SNN.png")) +
  draw_label("subepithelial nerve net (SNN)", x = 0.5, y = 1.05, color = Okabe_Ito[5], size = 10, hjust = 0.5) +
  draw_label("aboral view", x = 0.18, y = 0.95, color = "black", size = 8, hjust = 0.5) +
  draw_label("lateral view", x = 0.67, y = 0.95, color = "black", size = 8, hjust = 0.5) +
  draw_label("Q1", x = 0.3, y = 0.9, color = "black", size = 9, hjust = 0.5, alpha = 0.5) +
  draw_label("Q2", x = 0.06, y = 0.79, color = "black", size = 9, hjust = 0.5, alpha = 0.5) +
  draw_label("Q3", x = 0.03, y = 0.08, color = "black", size = 9, hjust = 0.5, alpha = 0.5) +
  draw_label("Q4", x = 0.275, y = 0.11, color = "black", size = 9, hjust = 0.5, alpha = 0.5) +
  draw_label("sagittal plane", x = 0.5, y = 0.86, color = "black", size = 7.5, hjust = 0.5) +
  draw_label("tentacular plane", x = 0.83, y = 0.86, color = "black", size = 7.5, hjust = 0.5) +
  draw_label("A", x = 0.67, y = 0.85, size = 7, color = "black", hjust = 0.5) +
  draw_label("O", x = 0.67, y = 0.65, size = 7, color = "black", hjust = 0.5) +
  draw_line(x = c(0.67, 0.67), y = c(0.69, 0.81), color = "black", linewidth = 0.65) +
  draw_label(expression(paste("25 ", mu, "m")), x = 0.95, y = 0.1, color = "black", size = 7, hjust = 0.5) +
  draw_line(x = c(0.91, 0.99), y = c(0.05, 0.05), color = "black", linewidth = 0.7)

panel_ANN_SNN_EM <- ggdraw() + draw_image(readPNG("manuscript/pictures/comparison_ANN_SNN.png"))

panel_CAT_pic <- ggdraw() + draw_image(readPNG("manuscript/pictures/figure_sup2_1.png"))

panel_mono_vs_poly <- ggdraw() + draw_image(readPNG("manuscript/pictures/figure_sup2_2.png"))

panel_multi_nucleus <- ggdraw() + draw_image(readPNG("manuscript/pictures/multi_nucleus.png"))

layout <- "
A#BBB#C
#######
DDD#EEE
"

Fig2_Sup1 <- panel_CAT_pic + panel_mono_vs_poly + panel_multi_nucleus +
  panel_SNN_3d + panel_ANN_SNN_EM +
  plot_layout(
    design = layout,
    heights = c(1.1, 0.25, 1),
    widths = c(0.7, 0.1, 0.7, 0.1, 0.2, 0.1, 1)
  ) +
  patchwork::plot_annotation(tag_levels = "A") &
  ggplot2::theme(plot.tag = element_text(
    size = 12,
    face = "plain", color = "black"
  ))

ggsave("manuscript/figures/Figure2_Supplement1.png",
  limitsize = FALSE,
  units = c("px"), Fig2_Sup1, width = 2850, height = 1250, bg = "white"
)

ggsave("manuscript/figures/Figure2_Supplement1.pdf",
  limitsize = FALSE,
  units = c("px"), Fig2_Sup1, width = 2850, height = 1250
)
