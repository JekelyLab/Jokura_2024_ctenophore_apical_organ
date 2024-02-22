# Code to generate Figure 1 of the Jokura et al 2024 Ctenophore apical organ connectome paper

# source packages and functions -----------------
source("analysis/scripts/packages_and_functions.R")

# load neurons ----------------------

Q1 <- nlapply(read.neurons.catmaid("Q1", pid = 35, conn = conn_http1),
              function(x) smooth_neuron(x, sigma = 1000)
              )
Q2 <- nlapply(read.neurons.catmaid("Q2", pid = 35, conn = conn_http1),
              function(x) smooth_neuron(x, sigma = 1000)
              )
Q3 <- nlapply(read.neurons.catmaid("Q3", pid = 35, conn = conn_http1),
              function(x) smooth_neuron(x, sigma = 1000)
              )
Q4 <- nlapply(read.neurons.catmaid("Q4", pid = 35, conn = conn_http1),
              function(x) smooth_neuron(x, sigma = 1000)
              )
balancer <- nlapply(read.neurons.catmaid("celltype:balancer", pid = 35, conn = conn_http1),
                     function(x) smooth_neuron(x, sigma = 1000)
                    )
LB <- nlapply(read.neurons.catmaid("celltype:lamellate", pid = 35, conn = conn_http1),
                    function(x) smooth_neuron(x, sigma = 1000)
                    )
syn_neuron <- nlapply(read.neurons.catmaid("celltype:neuron", pid = 35, conn = conn_http1),
                    function(x) smooth_neuron(x, sigma = 1000)
                    )
monocilia <- nlapply(read.neurons.catmaid("monociliated_cell", pid = 35, conn = conn_http1),
                    function(x) smooth_neuron(x, sigma = 1000)
                    )

balancer_Q1 <- nlapply(
  read.neurons.catmaid(skids_by_2annotations("celltype:balancer", "Q1"),
                       pid = 35, conn = conn_http1), 
  function(x) smooth_neuron(x, sigma = 1000)
  )
balancer_Q2 <- nlapply(
  read.neurons.catmaid(skids_by_2annotations("celltype:balancer", "Q2"),
                       pid = 35, conn = conn_http1), 
  function(x) smooth_neuron(x, sigma = 1000)
  )
balancer_Q3 <- nlapply(
  read.neurons.catmaid(skids_by_2annotations("celltype:balancer", "Q3"),
                       pid = 35, conn = conn_http1), 
  function(x) smooth_neuron(x, sigma = 1000)
  )
balancer_Q4 <- nlapply(
  read.neurons.catmaid(skids_by_2annotations("celltype:balancer", "Q4"),
                       pid = 35, conn = conn_http1), 
  function(x) smooth_neuron(x, sigma = 1000)
  )


# plot neurons -----------------------------------------------


plot_background()
plot3d(Q1, soma = T, lwd = c(2,6), add = T, alpha = c(0.5:1), col = Okabe_Ito[1])
plot3d(Q2, soma = T, lwd = c(2,6), add = T, alpha = c(0.5:1), col = Okabe_Ito[2])
plot3d(Q3, soma = T, lwd = c(2,6), add = T, alpha = c(0.5:1), col = Okabe_Ito[3])
plot3d(Q4, soma = T, lwd = c(2,6), add = T, alpha = c(0.5:1), col = Okabe_Ito[4])
rgl.snapshot("manuscript/pictures/Quadrants.png")
close3d()


plot_background()
plot3d(Q1, soma = T, lwd = c(2,6), add = T, alpha = 0.2, col = Okabe_Ito[1])
plot3d(Q2, soma = T, lwd = c(2,6), add = T, alpha = 0.2, col = Okabe_Ito[1])
plot3d(Q3, soma = T, lwd = c(2,6), add = T, alpha = 0.2, col = Okabe_Ito[1])
plot3d(Q4, soma = T, lwd = c(2,6), add = T, alpha = 0.2, col = Okabe_Ito[1])
plot3d(balancer_Q1, soma = T, lwd = c(2,6), add = T, alpha = 1, col = Okabe_Ito[2])
plot3d(balancer_Q2, soma = T, lwd = c(2,6), add = T, alpha = 1, col = Okabe_Ito[3])
plot3d(balancer_Q3, soma = T, lwd = c(2,6), add = T, alpha = 1, col = Okabe_Ito[5])
plot3d(balancer_Q4, soma = T, lwd = c(2,6), add = T, alpha = 1, col = Okabe_Ito[7])

rgl.snapshot("manuscript/pictures/balancer.png")

nview3d("ventral", extramat=rotationMatrix(0.9, 0.3, -0.1, 1))
rgl.snapshot("manuscript/pictures/balancer_ventral.png")

close3d()



# assemble figure -------------------------------------------------------------

#read pics
img_Q1_4 <- readPNG("manuscript/pictures/Quadrants.png")
img_balancer <- readPNG("manuscript/pictures/balancer.png")
img_balancer_v <- readPNG("manuscript/pictures/balancer_ventral.png")

#make panels
panel_Q1_4 <- ggdraw() + draw_image(img_Q1_4) +
  draw_label("Apical organ", x = 0.3, y = 0.98, size = 10, fontface = "italic") +
  draw_label("anterior", x = 0.6, y = 0.98, size = 10, fontface = "plain")

panel_balancer <- ggdraw() + draw_image(img_balancer) +
  draw_label("balancer", x = 0.3, y = 0.98, size = 10, fontface = "italic") +
  draw_label("anterior", x = 0.53, y = 0.98, size = 10, fontface = "plain")

panel_balancer_v <- ggdraw() + draw_image(img_balancer_v) +
  draw_label("balancer", x = 0.3, y = 0.98, size = 10, fontface = "italic") +
  draw_label("lateral", x = 0.53, y = 0.98, size = 10, fontface = "plain")

# define layout with textual representation for pathchwork assembly of figure
layout <- "
ABCD
EFGH
"

Figure1 <- panel_Q1_4 + panel_balancer + panel_balancer_v + panel_Q1_4 + 
  panel_Q1_4 + panel_Q1_4 + panel_Q1_4 + panel_Q1_4 +
  plot_layout(design = layout, heights = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 12, face = "plain"))

ggsave("manuscript/figures/Figure1.png", limitsize = FALSE, 
       units = c("px"), Figure1, width = 3300, height = 1600, bg = "white")


ggsave("manuscript/figures/Figure1.pdf", limitsize = FALSE, 
       units = c("px"), Figure1, width = 3300, height = 1600)



