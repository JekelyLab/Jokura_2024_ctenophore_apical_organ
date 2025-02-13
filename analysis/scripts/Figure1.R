# Code to generate Figure 1 of the Jokura et al 2024 Ctenophore apical organ connectome paper

# source packages and functions -----------------
source("analysis/scripts/packages_and_functions.R")

<<<<<<< Updated upstream

# assemble figure -------------------------------------------------------------

library(jpeg)
=======
# load cells ----------------------------------

balancer <- read_smooth_neuron("celltype:balancer")

Q1 <- nlapply(
  read.neurons.catmaid("Q1", pid = 35),
  function(x)
  smooth_neuron(x, sigma = 1000)
  )

Q2 <- nlapply(
  read.neurons.catmaid("Q2", pid = 35),
  function(x)
  smooth_neuron(x, sigma = 1000)
  )

Q3 <- nlapply(
  read.neurons.catmaid("Q3", pid = 35),
  function(x)
  smooth_neuron(x, sigma = 1000)
  )

Q4 <- nlapply(
  read.neurons.catmaid("Q4", pid = 35),
  function(x)
  smooth_neuron(x, sigma = 1000)
  )

with_soma <- nlapply(
  read.neurons.catmaid("with_soma", pid = 35),
  function(x)
  smooth_neuron(x, sigma = 1000)
  )


# plot balancer and quadrants ----------------------



# plot balancer -----------------------------------------------

plot_background()
close3d()
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  #defines the two scenes
#define the size of the rgl window, the view and zoom
#par3d(windowRect = c(0, 0, 1200, 350))
par3d(windowRect = c(0, 0, 2400, 700))

#plot aboral view
plot3d(balancer,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)

#for (object in objects) {
#  plot3d(object, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.05, col = Okabe_Ito[8])
#}


plot3d(
  with_soma, soma = T, lwd = 1, add = T, 
  alpha = 0.05, col = Okabe_Ito[8]
  )

#texts3d(58000,51000,5000, "INRGW", cex = 3, col = "#56B4E9")

#aboral view
nview3d("anterior", extramat = rotationMatrix(1.05, 250, -200, 1000))
#rgl.snapshot("manuscript/pictures/balancer_aboral_view.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Sagittal plane
plot3d(balancer,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)

plot3d(
  with_soma, soma = T, lwd = 1, add = T, 
  alpha = 0.05, col = Okabe_Ito[8]
  )

nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)

#plot lateral view of Tentacular plane
plot3d(balancer,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)

plot3d(
  with_soma, soma = T, lwd = 1, add = T, 
  alpha = 0.05, col = Okabe_Ito[8]
  )

nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)

#make a snapshot to the working directory
rgl.snapshot("manuscript/pictures/balancer.png")
close3d()



# assemble figure -------------------------------------------------------------

#read pics
panel_balancer <- ggdraw() + draw_image(readPNG("manuscript/pictures/balancer.png")) +
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


 
panel_larva_5dpf <- ggdraw() + draw_image(readJPEG("manuscript/pictures/Mnemiopsis_larva_5dpf.jpg"))

>>>>>>> Stashed changes
library(magick)

panel_larva_pic <- ggdraw() + draw_image(image_read("manuscript/pictures/Mnemiopsis_larva_5dpf.png"))
panel_AO_pic <- ggdraw() + draw_image(image_read("manuscript/pictures/AO_mag_pics.png"))
panel_serial_sectioning <- ggdraw() + draw_image(readPNG("manuscript/pictures/serial sectioning.png"))
panel_catmaid_overview <- ggdraw() + draw_image(readPNG("manuscript/pictures/overview.png"))
<<<<<<< Updated upstream
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
=======
panel_Q <- ggdraw() + draw_image(readPNG("manuscript/pictures/quadrants.png"))

layout <- "
AABBCCDDD
#########
EEEFFGGGG
"




Figure1 <- panel_larva_5dpf + panel_AO + panel_AO_schema + panel_AO_dimention_schema + 
  panel_catmaid_overview + panel_Q + panel_balancer +
  patchwork::plot_layout(design = layout, 
                         heights = c(1,0.1,1), 
                         widths = c(1,1,1,1,1,1,1)) + 
  patchwork::plot_annotation(tag_levels = "A") &  
  ggplot2::theme(plot.tag = element_text(size = 12, 
                                         face='plain', color='black'))


ggsave("manuscript/figures/Fig1.png", limitsize = FALSE, 
       units = c("px"), Figure1, width = 2400, height = 1000, bg='white')  
>>>>>>> Stashed changes


ggsave("manuscript/figures/Figure1.pdf", limitsize = FALSE, 
       units = c("px"), Figure1, width = 2600, height = 1000) 

