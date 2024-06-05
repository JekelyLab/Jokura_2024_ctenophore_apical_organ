#network analysis and plotting

source("analysis/scripts/packages_and_functions.R")

# load cell type ----------------------------------

read_smooth_neuron <- function(annotation){
  nlapply(read.neurons.catmaid(annotation, pid = 35),
          function(x)
            smooth_neuron(x, sigma = 1000))
}

balancer <- read_smooth_neuron("celltype:balancer")
syn_neuron <- read_smooth_neuron("celltype:SSN")
Q1 <- read_smooth_neuron("Q1")
Q2 <- read_smooth_neuron("Q2")
Q3 <- read_smooth_neuron("Q3")
Q4 <- read_smooth_neuron("Q4")
all_celltypes <- read_smooth_neuron("celltype")

# retrieve connectors ----------------

conn_syn <- connectors(syn_neuron)
presyn_syn <- subset(conn_syn, prepost == 0)
postsyn_syn <- subset(conn_syn, prepost == 1)

# plot syn ----------------------
plot_background()

plot3d(
  syn_neuron, soma = FALSE, 
  color = c("#1f77b4", "#ff7f0e", "#2ca02c"), 
  alpha = 0.5, lwd = c(2,3,3)
)
plot3d(
  all_celltypes, soma = TRUE, color = "grey90", 
  alpha = 0.08, lwd = c(4,3,3)
)
plot3d(
  outline, WithNodes = F,
  add=T, alpha=0.07, col="#E2E2E2"
) 

# plot presynapses
plot3d(
  presyn_syn$x, 
  presyn_syn$y, 
  presyn_syn$z, 
  size = 7, alpha = 1, col = "red", 
  add = T
  )

# plot postsynapses
plot3d(
  postsyn_syn$x, 
  postsyn_syn$y, 
  postsyn_syn$z, 
  size = 10, alpha = 1, col = "magenta", 
  add = T
  )

anterior()
rgl.snapshot("manuscript/pictures/syn_prepost_ant.png")

left()
rgl.snapshot("manuscript/pictures/syn_prepost_left.png")

close3d()

# plot quadrants --------------------


# plot syn ----------------------
plot_background()

plot3d(
  Q1, soma = TRUE, 
  color = Okabe_Ito[1], 
  alpha = 0.7, lwd = c(1,2,3,4)
)
plot3d(
  Q2, soma = TRUE, 
  color = Okabe_Ito[2], 
  alpha = 0.7, lwd = c(1,2,3,4)
)
plot3d(
  Q3, soma = TRUE, 
  color = Okabe_Ito[3], 
  alpha = 0.7, lwd = c(1,2,3,4)
)
plot3d(
  Q4, soma = TRUE, 
  color = Okabe_Ito[4], 
  alpha = 0.7, lwd = c(1,2,3,4)
)



rgl.snapshot("manuscript/pictures/quadrants.png")

close3d()

# assemble figure ---------------------------------------------------------

panel_syn_oral <- ggdraw() + draw_image(readPNG("manuscript/pictures/syn_prepost_ant.png")) +
  draw_label("oral view", x = 0.75, y = 0.95, size = 10) +
  draw_label("nerve net", x = 0.1, y = 0.95, size = 10, hjust = 0)
panel_syn_side <- ggdraw() + draw_image(readPNG("manuscript/pictures/syn_prepost_left.png")) +
  draw_label("side view", x = 0.75, y = 0.95, size = 10)


layout <- "
AB"

Figure_nerve_net <-
  panel_syn_oral + panel_syn_side +
  plot_layout(design = layout) +
  patchwork::plot_layout(design = layout, 
                         heights = c(1,1,1,1,1), 
                         widths = c(1,1,1)) +
  patchwork::plot_annotation(tag_levels = c("A")) &
  ggplot2::theme(plot.tag = element_text(size = 12, face='plain', color='black'))


ggsave(
  "manuscript/figures/Figure_nerve_net.png", 
  limitsize = FALSE, units = c("px"), 
  Figure_nerve_net, width = 1600, height = 800, bg = "white"
  )

ggsave(
  "manuscript/figures/Figure_nerve_net.pdf", limitsize = FALSE, 
  units = c("px"), Figure_nerve_net, width = 1600, height = 800
  )


