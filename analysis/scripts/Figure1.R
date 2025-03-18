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

# read pics
#panel_balancer <- ggdraw() + draw_image(readPNG("manuscript/pictures/balancer.png")) +
#  draw_label("balancer cells", x = 0.1, y = 0.98, size = 10, fontface = "plain") +
#  draw_label("aboral view", x = 0.1, y = 0.86, color = "black", size = 8) +
#  draw_label("lateral view of S plane", x = 0.45, y = 0.86, color = "black", size = 8) +
#  draw_label("lateral view of T plane", x = 0.75, y = 0.86, color = "black", size = 8) +
#  draw_line(x = c(0.85, 0.95), y = c(0.1, 0.1), color = "black", size = 0.5) +
#  draw_label(expression(paste("25 ", mu, " m")), x = 0.9, y = 0.14, color = "black", size = 7, hjust = 0.5) +
#  draw_label("S", x = 0.325, y = 0.16, size = 6, color = "black", hjust = 0.5) +
#  geom_segment(aes(x = 0.25, y = 0.16, xend = 0.31, yend = 0.16),
#    color = "black",
#    arrow = arrow(ends = "both", type = "closed", length = unit(0.1, "cm")),
#    lineend = "butt",
#    linejoin = "mitre",
#    arrow.fill = "black", linewidth = 0.175
#  ) +
#  draw_label("T", x = 0.28, y = 0.05, size = 6, color = "black", hjust = 0.5) +
#  geom_segment(aes(x = 0.28, y = 0.08, xend = 0.28, yend = 0.24),
#    color = "black",
#    arrow = arrow(ends = "both", type = "closed", length = unit(0.1, "cm")),
#    lineend = "butt",
#    linejoin = "mitre",
#    arrow.fill = "black", linewidth = 0.175
#  ) +
#  draw_label("A", x = 0.66, y = 0.25, size = 6, color = "black", hjust = 0.5) +
#  draw_label("O", x = 0.66, y = 0.05, size = 6, color = "black", hjust = 0.5) +
#  geom_segment(aes(x = 0.66, y = 0.09, xend = 0.66, yend = 0.21),
#    color = "black",
#    arrow = arrow(ends = "both", type = "closed", length = unit(0.1, "cm")),
#    lineend = "butt",
#    linejoin = "mitre",
#    arrow.fill = "black", linewidth = 0.175
#  )


#panel_larva_pic <- ggdraw() + draw_image(image_read("manuscript/pictures/Mnemiopsis_larva_24hpf.png")) +
#  draw_label("lateral view of TA plane", x = 0.5, y = 0.99, color = "black", size = 8) +
#  draw_line(x = c(0.75, 0.85), y = c(0.05, 0.05), color = "black", size = 0.5) +
#  draw_label(expression(paste("100 ", mu, " m")), x = 0.8, y = 0.09, color = "black", size = 7, hjust = 0.5)

#panel_AO_pic <- ggdraw() + draw_image(image_read("manuscript/pictures/AO_mag_pics.png")) +
#  draw_label("aboral view", x = 0.2, y = 0.99, color = "black", size = 8) +
#  draw_label("lateral view of SA plane", x = 0.5, y = 0.99, color = "black", size = 8) +
#  draw_label("lateral view of TA plane", x = 0.85, y = 0.99, color = "black", size = 8) +
#  draw_line(x = c(0.89, 0.99), y = c(0.05, 0.05), color = "black", size = 0.5) +
#  draw_label(expression(paste("20 ", mu, " m")), x = 0.94, y = 0.09, color = "black", size = 7, hjust = 0.5)

#panel_larva_schematic <- ggdraw() + draw_image(image_read("manuscript/pictures/larva_aboral_view_schematic.png")) +
#  draw_label("aboral view", x = 0.5, y = 0.97, color = "black", size = 8)

#panel_AO_schematic <- ggdraw() + draw_image(image_read("manuscript/pictures/AO_schematic.png")) +
#  draw_label("aboral view", x = 0.25, y = 0.97, color = "black", size = 8) +
#  draw_label("lateral view of SA plane", x = 0.5, y = 0.97, color = "black", size = 8) +
#  draw_label("lateral view of TA plane", x = 0.75, y = 0.97, color = "black", size = 8)

panel_catmaid_overview <- ggdraw() + draw_image(readPNG("manuscript/pictures/overview.png")) +
  draw_label("serial EM volume", x = 0.5, y = 0.925, size = 11, fontface = "plain", hjust = 0.5) +
  draw_label("619 sections", x = 0.025, y = 0.81, color = "black", size = 9, hjust = 0) +
  draw_label("927 cells", x = 0.025, y = 0.75, color = "black", size = 9, hjust = 0)

#panel_3d_all_cells <- ggdraw() + draw_image(readPNG("manuscript/pictures/all_cells_3_views.png"))

#panel_bal_ant <- ggdraw() + draw_image(readPNG("manuscript/pictures/balancer_Q1_4_aboral_view.png")) +
#  draw_label("four quadrants", x = 0.4, y = 0.98, size = 10, fontface = "plain", hjust = 0.5)

#panel_bal_side <- ggdraw() + draw_image(readPNG("manuscript/pictures/balancer_Q1_4_side_view.png"))

panel_quadrants <- ggdraw() + draw_image(readPNG("manuscript/pictures/quadrants.png")) +
#  draw_line(x = c(0, 1), y = c(0, 0), color = "black", linewidth = 0.25, alpha = 0.1) +
#  draw_line(x = c(0, 1), y = c(0.1, 0.1), color = "black", linewidth = 0.25, alpha = 0.1) +
#  draw_line(x = c(0, 1), y = c(0.2, 0.2), color = "black", linewidth = 0.25, alpha = 0.1) +
#  draw_line(x = c(0, 1), y = c(0.3, 0.3), color = "black", linewidth = 0.25, alpha = 0.1) +
#  draw_line(x = c(0, 1), y = c(0.4, 0.4), color = "black", linewidth = 0.25, alpha = 0.1) +
#  draw_line(x = c(0, 1), y = c(0.5, 0.5), color = "black", linewidth = 0.25, alpha = 0.1) +
#  draw_line(x = c(0, 1), y = c(0.6, 0.6), color = "black", linewidth = 0.25, alpha = 0.1) +
#  draw_line(x = c(0, 1), y = c(0.7, 0.7), color = "black", linewidth = 0.25, alpha = 0.1) +
#  draw_line(x = c(0, 1), y = c(0.8, 0.8), color = "black", linewidth = 0.25, alpha = 0.1) +
#  draw_line(x = c(0, 1), y = c(0.9, 0.9), color = "black", linewidth = 0.25, alpha = 0.1) +
#  draw_line(x = c(0, 1), y = c(1, 1), color = "black", linewidth = 0.25, alpha = 0.1) +
#  draw_line(x = c(0, 0), y = c(0, 1), color = "black", linewidth = 0.25, alpha = 0.1) +
#  draw_line(x = c(0.1, 0.1), y = c(0, 1), color = "black", linewidth = 0.25, alpha = 0.1) +
#  draw_line(x = c(0.2, 0.2), y = c(0, 1), color = "black", linewidth = 0.25, alpha = 0.1) +
#  draw_line(x = c(0.3, 0.3), y = c(0, 1), color = "black", linewidth = 0.25, alpha = 0.1) +
#  draw_line(x = c(0.4, 0.4), y = c(0, 1), color = "black", linewidth = 0.25, alpha = 0.1) +
#  draw_line(x = c(0.5, 0.5), y = c(0, 1), color = "black", linewidth = 0.25, alpha = 0.1) +
#  draw_line(x = c(0.6, 0.6), y = c(0, 1), color = "black", linewidth = 0.25, alpha = 0.1) +
#  draw_line(x = c(0.7, 0.7), y = c(0, 1), color = "black", linewidth = 0.25, alpha = 0.1) +
#  draw_line(x = c(0.8, 0.8), y = c(0, 1), color = "black", linewidth = 0.25, alpha = 0.1) +
#  draw_line(x = c(0.9, 0.9), y = c(0, 1), color = "black", linewidth = 0.25, alpha = 0.1) +
#  draw_line(x = c(1, 1), y = c(0, 1), color = "black", linewidth = 0.25, alpha = 0.1) +
  draw_label("Q1", x = 0.3, y = 0.9, color = Okabe_Ito[1], size = 13, hjust = 0.5, alpha = 1) +
  draw_label("Q2", x = 0.06, y = 0.79, color = Okabe_Ito[2], size = 13, hjust = 0.5, alpha = 1) +
  draw_label("Q3", x = 0.03, y = 0.08, color = Okabe_Ito[6], size = 13, hjust = 0.5, alpha = 1) +
  draw_label("Q4", x = 0.275, y = 0.11, color = Okabe_Ito[7], size = 13, hjust = 0.5, alpha = 1) +
  draw_label("aboral view", x = 0.17, y = 1, color = "black", size = 11, hjust = 0.5) +
  draw_label("sagittal plane", x = 0.5, y = 1, color = "black", size = 11, hjust = 0.5) +
  draw_label("tentacular plane", x = 0.83, y = 1, color = "black", size = 11, hjust = 0.5)  +
  draw_label(expression(paste("25 ", mu, "m")), x = 0.95, y = 0.1, color = "black", size = 10, hjust = 0.5) +
  draw_line(x = c(0.91, 0.99), y = c(0.05, 0.05), color = "black", linewidth = 1) +
#  draw_label("T", x = 0.33, y = 0.1, size = 10, color = "black", hjust = 0.5) +
#  draw_line(x = c(0.3, 0.32), y = c(0.1, 0.1), color = "black", linewidth = 0.75) +
#  draw_label("S", x = 0.31, y = 0.03, size = 10, color = "black", hjust = 0.5) +
#  draw_line(x = c(0.31, 0.31), y = c(0.07, 0.13), color = "black", linewidth = 0.75) +
  draw_label("A", x = 0.66, y = 0.25, size = 10, color = "black", hjust = 0.5) +
  draw_label("O", x = 0.66, y = 0.05, size = 10, color = "black", hjust = 0.5) +
  draw_line(x = c(0.66, 0.66), y = c(0.09, 0.21), color = "black", linewidth = 0.75) +
  draw_label("li", x = 0.5, y = 0.88, color = "black", size = 10, hjust = 0) +
  draw_label("bal", x = 0.56, y = 0.88, color = "black", size = 10, hjust = 0) 


# create panel A, B, create separately for precise alignment
#Fig1A_B <- plot_grid(panel_larva_pic, NULL, panel_AO_pic,
#                     ncol = 3,
#                     align = "h",
#                     labels = c("A", "", "B"),
#                     label_size = 12, label_y = 1.052, label_x = 0,
#                     label_fontfamily = "sans", label_fontface = "plain",
#                     # A negative rel_height shrinks space between elements
#                     rel_widths = c(1, 0.05, 3.747),
#                     rel_heights = c(1)
#)



#layout <- "
#AA
###
#BC
###
#DE
#"

#Figure1 <- Fig1A_B +
#  panel_catmaid_overview + panel_schematic +
#  panel_bal_ant + panel_balancer +
#  plot_layout(
#    design = layout,
#    heights = c(1, 0.05, 1, 0.05, 1),
#    widths = c(1,3)
#  ) +
#  plot_annotation(tag_levels = list(c("","C","D","E","F"))) &
#  theme(plot.tag = element_text(size = 12, face = "plain"))




panel_larva_pic <- ggdraw() + draw_image(image_read("manuscript/pictures/Mnemiopsis_larva_24hpf.png"))

panel_AO_pic <- ggdraw() + draw_image(image_read("manuscript/pictures/AO_mag_pics.png"))

panel_larva_schematic <- ggdraw() + draw_image(image_read("manuscript/pictures/larva_aboral_view_schematic.png"))

panel_AO_schematic <- ggdraw() + draw_image(image_read("manuscript/pictures/AO_schematic.png"))


layout <- "
ACCC
####
BDDD
####
EFFF
"

Figure1 <- panel_larva_pic + panel_larva_schematic +  
  panel_AO_pic + panel_AO_schematic + 
  panel_catmaid_overview + panel_quadrants +
  plot_layout(
    design = layout,
    heights = c(1, 0.1, 1, 0.1, 1),
    widths = c(1, 1, 1, 1)
  ) + 
  patchwork::plot_annotation(tag_levels = "A") &  
  ggplot2::theme(plot.tag = element_text(size = 16, 
                                         face='plain', color='black'))


ggsave("manuscript/figures/Figure1.png",
  limitsize = FALSE,
  units = c("px"), Figure1, width = 3200, height = 2200, bg = "white"
)

ggsave("manuscript/figures/Figure1.pdf",
  limitsize = FALSE,
  units = c("px"), Figure1, width = 3200, height = 2200
)
