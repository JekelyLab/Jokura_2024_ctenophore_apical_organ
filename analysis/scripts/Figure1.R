# Code to generate Figure 1 of the Jokura et al 2024 Ctenophore apical organ connectome paper

# source packages and functions ------------------------------------------------
source("analysis/scripts/packages_and_functions.R")

# load cells -------------------------------------------------------------------

balancer <- read_smooth_neuron("celltype:balancer")
cgroove <- read_smooth_neuron("celltype:Cgroove")
statolith <- read_smooth_neuron("statolith")

Q1 <- read_smooth_neuron("Q1")
Q2 <- read_smooth_neuron("Q2")
Q3 <- read_smooth_neuron("Q3")
Q4 <- read_smooth_neuron("Q4")

# plot quadrants ----------------------------------------------------------------

plot_views <- function(view_func) {
  plot3d(balancer, soma = FALSE, lwd = 2.5, add = TRUE, alpha = 0.75, col = "gray40", WithConnectors = FALSE)
  plot3d(statolith, soma = TRUE, lwd = 0, add = TRUE, alpha = 1, col = "gray", WithConnectors = FALSE)
  plot3d(Q1, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.4, col = Okabe_Ito[1])
  plot3d(Q2, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.4, col = Okabe_Ito[2])
  plot3d(Q3, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.4, col = Okabe_Ito[6])
  plot3d(Q4, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.4, col = Okabe_Ito[7])
  plot3d(outline, add = TRUE, alpha = 0.04, col = Okabe_Ito[8])
  view_func()
  par3d(zoom = 0.61)
  next3d(clear = FALSE)
}

close3d()
nopen3d()
mfrow3d(1, 3)
par3d(windowRect = c(0, 0, 2400, 700))

plot_views(aboral)
plot_views(sagittal)
plot_views(tentacular)

rgl.snapshot("manuscript/pictures/quadrants.png")
close3d()


# assemble figure -------------------------------------------------------------

panel_catmaid_overview <- ggdraw() + draw_image(readPNG("manuscript/pictures/overview.png")) +
  draw_label("serial EM volume", x = 0.5, y = 0.99, size = 10, fontface = "plain", hjust = 0.5) +
  draw_label("619 sections", x = 0.025, y = 0.81, color = "black", size = 9, hjust = 0) +
  draw_label("927 cells", x = 0.025, y = 0.75, color = "black", size = 9, hjust = 0) +
  draw_line(x = c(0.75, 0.9), y = c(0.05, 0.05), color = "black", size = 1) +
  draw_label(expression(paste("10 ", mu, "m")), x = 0.825, y = 0.09, color = "black", size = 8, hjust = 0.5)


panel_quadrants <- ggdraw() + draw_image(readPNG("manuscript/pictures/quadrants.png")) +
  draw_label("Q1", x = 0.3, y = 0.9, color = Okabe_Ito[1], size = 13, hjust = 0.5, alpha = 1) +
  draw_label("Q2", x = 0.06, y = 0.79, color = Okabe_Ito[2], size = 13, hjust = 0.5, alpha = 1) +
  draw_label("Q3", x = 0.03, y = 0.08, color = Okabe_Ito[6], size = 13, hjust = 0.5, alpha = 1) +
  draw_label("Q4", x = 0.275, y = 0.11, color = Okabe_Ito[7], size = 13, hjust = 0.5, alpha = 1) +
  draw_label("aboral view", x = 0.17, y = 0.99, color = "black", size = 10, hjust = 0.5) +
  draw_label("lateral view", x = 0.66, y = 0.99, color = "black", size = 10, hjust = 0.5) +
  draw_label("sagittal plane", x = 0.5, y = 0.95, color = "black", size = 10, hjust = 0.5) +
  draw_label("tentacular plane", x = 0.83, y = 0.95, color = "black", size = 10, hjust = 0.5)  +
  draw_label(expression(paste("25 ", mu, "m")), x = 0.95, y = 0.1, color = "black", size = 10, hjust = 0.5) +
  draw_line(x = c(0.91, 0.99), y = c(0.05, 0.05), color = "black", linewidth = 1) +
  draw_label("A", x = 0.66, y = 0.25, size = 10, color = "black", hjust = 0.5) +
  draw_label("O", x = 0.66, y = 0.05, size = 10, color = "black", hjust = 0.5) +
  draw_line(x = c(0.66, 0.66), y = c(0.09, 0.21), color = "black", linewidth = 0.75) +
  draw_label("li", x = 0.5, y = 0.88, color = "black", size = 10, hjust = 0) +
  draw_label("bal", x = 0.56, y = 0.88, color = "black", size = 10, hjust = 0) 

panel_larva_pic <- ggdraw() + draw_image(image_read("manuscript/pictures/Mnemiopsis_larva_24hpf.png"))

panel_AO_aov <- ggdraw() + draw_image(image_read("manuscript/pictures/AO_mag_pics_schem_aboral_view.png"))

panel_AO_pic <- ggdraw() + draw_image(image_read("manuscript/pictures/AO_mag_pics.png"))
panel_AO_pic_lv <- ggdraw() + draw_image(image_read("manuscript/pictures/AO_mag_pics_lateral_view.png")) +
  draw_label("lateral view", x = 0.5, y = 0.99, color = "black", size = 10, hjust = 0.5) +
  draw_label("sagittal plane", x = 0.25, y = 0.95, color = "black", size = 10, hjust = 0.5) +
  draw_label("tentacular plane", x = 0.75, y = 0.95, color = "black", size = 10, hjust = 0.5) 
  

panel_larva_schematic <- ggdraw() + draw_image(image_read("manuscript/pictures/larva_aboral_view_schematic.png"))

panel_AO_schem <- ggdraw() + draw_image(image_read("manuscript/pictures/AO_schematic.png"))
panel_AO_schem_lv <- ggdraw() + draw_image(image_read("manuscript/pictures/AO_schematic_lateral_view.png")) +
  draw_label("lateral view", x = 0.5, y = 0.99, color = "black", size = 10, hjust = 0.5) +
  draw_label("sagittal plane", x = 0.25, y = 0.95, color = "black", size = 10, hjust = 0.5) +
  draw_label("tentacular plane", x = 0.75, y = 0.95, color = "black", size = 10, hjust = 0.5) 


layout <- "
ABCC
####
DDEE
####
FGGG
"

Figure1 <- panel_larva_pic + panel_larva_schematic + panel_AO_aov +  
  panel_AO_pic_lv + panel_AO_schem_lv + 
  panel_catmaid_overview + panel_quadrants +
  plot_layout(
    design = layout,
    heights = c(1, 0.05, 1.2, 0.01, 1),
    widths = c(1, 1, 1, 1)
  ) + 
  patchwork::plot_annotation(tag_levels = "A") &  
  ggplot2::theme(plot.tag = element_text(size = 16, 
                                         face='plain', color='black'))


ggsave("manuscript/figures/Figure1.png",
  limitsize = FALSE,
  units = c("px"), Figure1, width = 3200, height = 2500, bg = "white"
)

ggsave("manuscript/figures/Figure1.pdf",
  limitsize = FALSE,
  units = c("px"), Figure1, width = 3200, height = 2500
)
