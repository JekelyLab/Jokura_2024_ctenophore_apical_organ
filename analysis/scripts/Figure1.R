# Code to generate Figure 1 of the Jokura et al 2024 Ctenophore apical organ connectome paper

# source packages and functions -----------------
source("analysis/scripts/packages_and_functions.R")

# load cells ----------------------------------

balancer <- read_smooth_neuron("celltype:balancer")
Q1 <- read_smooth_neuron("Q1")
Q2 <- read_smooth_neuron("Q2")
Q3 <- read_smooth_neuron("Q3")
Q4 <- read_smooth_neuron("Q4")
with_soma <- read_smooth_neuron("with_soma")
lithocyte <- read_smooth_neuron("celltype:lithocyte")

# plot quadrants ----------------------

plot_background()
close3d()
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  #defines the two scenes
#define the size of the rgl window, the view and zoom
#par3d(windowRect = c(0, 0, 1200, 350))
par3d(windowRect = c(0, 0, 2400, 700))

#plot aboral view



#plot3d(
#  balancer,
#  soma = T, lwd = 1, add = T,
#  alpha = 0.9, col = Okabe_Ito[5],
#  WithConnectors = F
#)
plot3d(
  lithocyte,
  soma = T, lwd = 1, add = T,
  alpha = 0.5, col = Okabe_Ito[8],
  WithConnectors = F
)
plot3d(
  Q1,
  soma = T, lwd = 1, add = T,
  alpha = 0.15, col = Okabe_Ito[1]
)
plot3d(
  Q2,
  soma = T, lwd = 1, add = T,
  alpha = 0.15, col = Okabe_Ito[2]
)
plot3d(
  Q3,
  soma = T, lwd = 1, add = T,
  alpha = 0.15, col = Okabe_Ito[6]
)
plot3d(
  Q4,
  soma = T, lwd = 1, add = T,
  alpha = 0.15, col = Okabe_Ito[7]
)

#plot3d(
#  outline,
#  add = T,
#  alpha = 0.06, col = "grey50"
#)


nview3d("anterior", extramat = rotationMatrix(1.05, 250, -200, 1000))



text3d(15402, 56436, -200, "lithocyte", cex = 2.5)
#text3d(38002, 65436, 10000, "balancer", cex = 2.5)
text3d(23402, 29436, -200, "Q1", cex = 2.5)
text3d(11602, 54436, -200, "Q2", cex = 2.5)
text3d(23402, 65436, -200, "Q3", cex = 2.5)
text3d(36402, 35436, -200, "Q4", cex = 2.5)


#rgl.snapshot("manuscript/pictures/quadrants_aboral_view.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Sagittal plane


plot3d(
  lithocyte,
  soma = T, lwd = 1, add = T,
  alpha = 0.5, col = Okabe_Ito[8],
  WithConnectors = F
)
plot3d(
  Q1,
  soma = T, lwd = 1, add = T,
  alpha = 0.15, col = Okabe_Ito[1]
)
plot3d(
  Q2,
  soma = T, lwd = 1, add = T,
  alpha = 0.15, col = Okabe_Ito[2]
)
plot3d(
  Q3,
  soma = T, lwd = 1, add = T,
  alpha = 0.15, col = Okabe_Ito[6]
)
plot3d(
  Q4,
  soma = T, lwd = 1, add = T,
  alpha = 0.15, col = Okabe_Ito[7]
)

#plot3d(
#  outline,
#  add = T,
#  alpha = 0.06, col = "grey50"
#)


nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))

text3d(23402, 29436, -200, "Q1", cex = 2.5)
text3d(23402, 29436, 2000, "Q4", cex = 2.5)
text3d(33402, 70436, -7000, "Q2", cex = 2.5)
text3d(33402, 70436, -5000, "Q3", cex = 2.5)

#rgl.snapshot("manuscript/pictures/quadrants_viewing_sagittal_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Tentacular plane

plot3d(
  lithocyte,
  soma = T, lwd = 1, add = T,
  alpha = 0.5, col = Okabe_Ito[8],
  WithConnectors = F
)
plot3d(
  Q1,
  soma = T, lwd = 1, add = T,
  alpha = 0.15, col = Okabe_Ito[1]
)
plot3d(
  Q2,
  soma = T, lwd = 1, add = T,
  alpha = 0.15, col = Okabe_Ito[2]
)
plot3d(
  Q3,
  soma = T, lwd = 1, add = T,
  alpha = 0.15, col = Okabe_Ito[6]
)
plot3d(
  Q4,
  soma = T, lwd = 1, add = T,
  alpha = 0.15, col = Okabe_Ito[7]
)

#plot3d(
#  outline,
#  add = T,
#  alpha = 0.06, col = "grey50"
#)

#plot lateral view of Tentacular plane
nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))

text3d(23402, 25436, -200, "Q1", cex = 2.5)
text3d(23402, 25436, 2000, "Q2", cex = 2.5)
text3d(33402, 50436, -4000, "Q3", cex = 2.5)
text3d(33402, 50436, -2000, "Q4", cex = 2.5)

#rgl.snapshot("manuscript/pictures/quadrants_viewing_tentacular_plane.png")
par3d(zoom=0.61)

#make a snapshot to the working directory
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


panel_larva_pic <- ggdraw() + draw_image(image_read("manuscript/pictures/Mnemiopsis_larva_24hpf.png")) +
  draw_label("lateral view of TA plane", x = 0.5, y = 0.99, color = "black", size = 8) +
  draw_line(x = c(0.75, 0.85), y = c(0.05, 0.05), color = "black", size = 0.5) +
  draw_label(expression(paste("100 ", mu, " m")), x = 0.8, y = 0.09, color = "black", size = 7, hjust = 0.5)

panel_AO_pic <- ggdraw() + draw_image(image_read("manuscript/pictures/AO_mag_pics.png")) +
  draw_label("aboral view", x = 0.2, y = 0.99, color = "black", size = 8) +
  draw_label("lateral view of SA plane", x = 0.5, y = 0.99, color = "black", size = 8) +
  draw_label("lateral view of TA plane", x = 0.85, y = 0.99, color = "black", size = 8) +
  draw_line(x = c(0.89, 0.99), y = c(0.05, 0.05), color = "black", size = 0.5) +
  draw_label(expression(paste("20 ", mu, " m")), x = 0.94, y = 0.09, color = "black", size = 7, hjust = 0.5)

panel_larva_schematic <- ggdraw() + draw_image(image_read("manuscript/pictures/larva_aboral_view_schematic.png")) +
  draw_label("aboral view", x = 0.5, y = 0.97, color = "black", size = 8)

panel_AO_schematic <- ggdraw() + draw_image(image_read("manuscript/pictures/AO_schematic.png")) +
  draw_label("aboral view", x = 0.25, y = 0.97, color = "black", size = 8) +
  draw_label("lateral view of SA plane", x = 0.5, y = 0.97, color = "black", size = 8) +
  draw_label("lateral view of TA plane", x = 0.75, y = 0.97, color = "black", size = 8)


panel_catmaid_overview <- ggdraw() + draw_image(readPNG("manuscript/pictures/overview.png")) +
  draw_label("serial EM volume", x = 0.3, y = 0.98, size = 10, fontface = "plain") +
  draw_label("619 sections", x = 0.05, y = 0.9, color = "black", size = 8, hjust = 0) +
  draw_label("927 cells", x = 0.05, y = 0.8, color = "black", size = 8, hjust = 0)

#panel_3d_all_cells <- ggdraw() + draw_image(readPNG("manuscript/pictures/all_cells_3_views.png"))

#panel_bal_ant <- ggdraw() + draw_image(readPNG("manuscript/pictures/balancer_Q1_4_aboral_view.png")) +
#  draw_label("four quadrants", x = 0.4, y = 0.98, size = 10, fontface = "plain", hjust = 0.5)

#panel_bal_side <- ggdraw() + draw_image(readPNG("manuscript/pictures/balancer_Q1_4_side_view.png"))

panel_quadrants <- ggdraw() + draw_image(readPNG("manuscript/pictures/quadrants.png")) +
  draw_label("aboral view", x = 0.25, y = 0.97, color = "black", size = 8) +
  draw_label("lateral view of SA plane", x = 0.5, y = 0.97, color = "black", size = 8) +
  draw_label("lateral view of TA plane", x = 0.75, y = 0.97, color = "black", size = 8)  +
  draw_label("balancer cells", x = 0.5, y = 0.95, size = 8.5, fontface="bold", hjust = 0.5) +
  draw_label("aboral view", x = 0.01, y = 0.86, color="black", size = 6, fontface="plain", hjust = 0) +
  draw_label("lateral view of PA plane", x = 0.33, y = 0.86, color="black", size = 6, fontface="plain", hjust = 0) +
  draw_label("lateral view of TA plane", x = 0.69, y = 0.86, color="black", size = 6, fontface="plain", hjust = 0) +
  draw_line(x = c(0.85, 0.95), y = c(0.1, 0.1), color = "black", size = 0.5) +
  draw_label(expression(paste("25 ", mu, "m")), x = 0.9, y = 0.14, color = "black", size = 7, hjust = 0.5) +
  draw_label("TA", x = 0.325, y = 0.16, size = 6, color = "black", hjust = 0.5) +
  geom_segment(aes(x = 0.25, y = 0.16, xend = 0.31, yend = 0.16), color = "black", 
               arrow = arrow(ends = "both", type = "closed", length = unit(0.1,"cm")),
               lineend = "butt",
               linejoin = "mitre",
               arrow.fill = "black", linewidth = 0.175) +
  draw_label("PA", x = 0.28, y = 0.05, size = 6, color = "black", hjust = 0.5) +
  geom_segment(aes(x = 0.28, y = 0.08, xend = 0.28, yend = 0.24), color = "black", 
               arrow = arrow(ends = "both", type = "closed", length = unit(0.1,"cm")),
               lineend = "butt",
               linejoin = "mitre",
               arrow.fill = "black", linewidth = 0.175) +
  draw_label("A", x = 0.66, y = 0.25, size = 6, color = "black", hjust = 0.5) +
  draw_label("O", x = 0.66, y = 0.05, size = 6, color = "black", hjust = 0.5) +
  geom_segment(aes(x = 0.66, y = 0.09, xend = 0.66, yend = 0.21), color = "black", 
               arrow = arrow(ends = "both", type = "closed", length = unit(0.1,"cm")),
               lineend = "butt",
               linejoin = "mitre",
               arrow.fill = "black", linewidth = 0.175) 


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



layout <- "
A#BBB
#####
C#DDD
#####
E#FFF
"

Figure1 <- panel_larva_pic + panel_AO_pic + 
  panel_larva_schematic + panel_AO_schematic + 
  panel_catmaid_overview + panel_quadrants +
  plot_layout(
    design = layout,
    heights = c(1, 0.1, 1, 0.1, 1),
    widths = c(1, 0.1, 1, 1, 1)
  ) + 
  patchwork::plot_annotation(tag_levels = "A") &  
  ggplot2::theme(plot.tag = element_text(size = 12, 
                                         face='plain', color='black'))


ggsave("manuscript/figures/Figure1.png",
  limitsize = FALSE,
  units = c("px"), Figure1, width = 3200, height = 2250, bg = "white"
)

ggsave("manuscript/figures/Figure1.pdf",
  limitsize = FALSE,
  units = c("px"), Figure1, width = 3200, height = 2250
)
