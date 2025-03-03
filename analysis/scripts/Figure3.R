# Code to generate Figure 3 of the Jokura et al 2024 Ctenophore apical organ connectome paper

# source packages and functions ------------------------------------------------
source("analysis/scripts/packages_and_functions.R")


# load cell type ----------------------------------

balancer <- read_smooth_neuron("celltype:balancer")
bridge <- read_smooth_neuron("celltype:bridge")
bristle <- read_smooth_neuron("celltype:bristle")
Cgroove_sag <- read_smooth_neuron("celltype:Cgroove-sag")
Cgroove_tag <- read_smooth_neuron("celltype:Cgroove-tag")
dense_vesicle <- read_smooth_neuron("celltype:dense_vesicle")
dome <- read_smooth_neuron("celltype:dome")
epithelial_floor <- read_smooth_neuron("celltype:epithelial_floor")
intra_multi_ciliated <- read_smooth_neuron("celltype:intra-multi-ciliated")
lamellate <- read_smooth_neuron("celltype:lamellate")
lithocyte <- read_smooth_neuron("celltype:lithocyte")
plumose <- read_smooth_neuron("celltype:plumose")
SSN <- read_smooth_neuron("celltype:SSN")

monociliated <- read_smooth_neuron("celltype:monociliated")
biciliated <- read_smooth_neuron("celltype:biciliated")
multiciliated <- read_smooth_neuron("celltype:multiciliated")
nonciliated <- read_smooth_neuron("celltype:nonciliated")

Q1 <- read_smooth_neuron("Q1")
Q2 <- read_smooth_neuron("Q2")
Q3 <- read_smooth_neuron("Q3")
Q4 <- read_smooth_neuron("Q4")

all_celltypes <- list(balancer,
                      bridge,
                      bristle,
                      Cgroove_sag,
                      Cgroove_tag,
                      dense_vesicle,
                      dome,
                      epithelial_floor,
                      intra_multi_ciliated,
                      lamellate,
                      lithocyte,
                      plumose,
                      SSN,
                      monociliated,
                      biciliated,
                      multiciliated,
                      nonciliated)




balancer_Q1 <- read_smooth_neuron(get_skids_with_annot(pid = 35, c("celltype:balancer", "Q1")))
balancer_Q2 <- read_smooth_neuron(get_skids_with_annot(pid = 35, c("celltype:balancer", "Q2")))
balancer_Q3 <- read_smooth_neuron(get_skids_with_annot(pid = 35, c("celltype:balancer", "Q3")))
balancer_Q4 <- read_smooth_neuron(get_skids_with_annot(pid = 35, c("celltype:balancer", "Q4")))
















# plot balancer Gaspar-----------------------------------------------

# 3d plotting of cells
nopen3d()
mfrow3d(1, 3) # defines the two scenes
# define the size of the rgl window, the view and zoom
# par3d(windowRect = c(0, 0, 1200, 350))
par3d(windowRect = c(0, 0, 2400, 800))

# plot aboral view
plot3d(
  balancer,
  soma = T, lwd = 1, add = T,
  alpha = 0.5, col = Okabe_Ito[5],
  WithConnectors = F
)
plot3d(
  with_soma,
  soma = T, lwd = 1, add = T,
  alpha = 0.05, col = Okabe_Ito[8]
)

# aboral view
nview3d("anterior", extramat = rotationMatrix(1.05, 250, -200, 1000))
par3d(zoom = 0.7)
# y-axis clip
clipplanes3d(1, 0, 0, -11500)
# x-axis clip
clipplanes3d(0, 1, 0, -24000)

# move to next panel in rgl window
next3d(clear = F)
# plot lateral view of Sagittal plane
plot3d(
  balancer,
  soma = T, lwd = 1, add = T,
  alpha = 0.5, col = Okabe_Ito[5],
  WithConnectors = F
)
plot3d(
  with_soma,
  soma = T, lwd = 1, add = T,
  alpha = 0.05, col = Okabe_Ito[8]
)

nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
par3d(zoom = 0.7)
# y-axis clip
clipplanes3d(1, 0, 0, -11500)
# x-axis clip
clipplanes3d(0, 1, 0, -24000)

# move to next panel in rgl window
next3d(clear = F)

# plot lateral view of Tentacular plane
plot3d(
  balancer,
  soma = T, lwd = 1, add = T,
  alpha = 0.5, col = Okabe_Ito[5],
  WithConnectors = F
)
plot3d(
  with_soma,
  soma = T, lwd = 1, add = T,
  alpha = 0.05, col = Okabe_Ito[8]
)

nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))
par3d(zoom = 0.7)
# y-axis clip
clipplanes3d(1, 0, 0, -11500)
# x-axis clip
clipplanes3d(0, 1, 0, -24000)

# make a snapshot to the working directory
rgl.snapshot("manuscript/pictures/balancer.png")
close3d()




















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


for (i in 1:length(all_celltypes)) {
  print(i)
  plot3d(
    all_celltypes[[i]], soma = TRUE, lwd = 1, add = TRUE, 
    alpha = 0.05, col = Okabe_Ito[8]
  )
}

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

for (i in 1:length(all_celltypes)) {
  print(i)
  plot3d(
    all_celltypes[[i]], soma = TRUE, lwd = 1, add = TRUE, 
    alpha = 0.05, col = Okabe_Ito[8]
  )
}



nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
#rgl.snapshot("manuscript/pictures/balancer_sagittal_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Tentacular plane
plot3d(balancer,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)

for (i in 1:length(all_celltypes)) {
  print(i)
  plot3d(
    all_celltypes[[i]], soma = TRUE, lwd = 1, add = TRUE, 
    alpha = 0.05, col = Okabe_Ito[8]
  )
}


nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))
#rgl.snapshot("manuscript/pictures/balancer_tentacular_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)


#make a snapshot to the working directory
rgl.snapshot("manuscript/pictures/balancer.png")


close3d()






# plot bridge -----------------------------------------------

plot_background()
close3d()
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  #defines the two scenes
#define the size of the rgl window, the view and zoom
#par3d(windowRect = c(0, 0, 1200, 350))
par3d(windowRect = c(0, 0, 2400, 700))

#plot aboral view
plot3d(bridge,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)


for (i in 1:length(all_celltypes)) {
  print(i)
  plot3d(
    all_celltypes[[i]], soma = TRUE, lwd = 1, add = TRUE, 
    alpha = 0.05, col = Okabe_Ito[8]
  )
}

#texts3d(58000,51000,5000, "INRGW", cex = 3, col = "#56B4E9")

#aboral view
nview3d("anterior", extramat = rotationMatrix(1.05, 250, -200, 1000))
#rgl.snapshot("manuscript/pictures/balancer_aboral_view.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Sagittal plane
plot3d(bridge,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)

for (i in 1:length(all_celltypes)) {
  print(i)
  plot3d(
    all_celltypes[[i]], soma = TRUE, lwd = 1, add = TRUE, 
    alpha = 0.05, col = Okabe_Ito[8]
  )
}



nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
#rgl.snapshot("manuscript/pictures/balancer_sagittal_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Tentacular plane
plot3d(bridge,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)

for (i in 1:length(all_celltypes)) {
  print(i)
  plot3d(
    all_celltypes[[i]], soma = TRUE, lwd = 1, add = TRUE, 
    alpha = 0.05, col = Okabe_Ito[8]
  )
}


nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))
#rgl.snapshot("manuscript/pictures/balancer_tentacular_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)


#make a snapshot to the working directory
rgl.snapshot("manuscript/pictures/bridge.png")


close3d()












# assemble figure -------------------------------------------------------------


panel_AO_matrix <- ggdraw() + draw_image(readPNG("manuscript/pictures/AO_neurons_matrix.png"))
panel_bal_3d <- ggdraw() + draw_image(readPNG("manuscript/pictures/balancer.png"))
panel_bal_syn <- ggdraw() + draw_image(readPNG("manuscript/pictures/SSN_prepost_synapse_balancer.png"))
panel_grav_nc_graph <- ggdraw() + draw_image(readPNG("manuscript/pictures/gravity_neural_circuit_graph.png"))

layout <- "
AA#BBBBB
########
CCCCC#DD
"

Figure3 <- panel_AO_matrix + panel_bal_3d + 
  panel_bal_syn + panel_grav_nc_graph +
  plot_layout(design = layout,
              heights = c(1, 0.1, 1),
              widths = c(1, 1, 0.1, 1, 1, 0.1, 1, 1)) + 
  plot_annotation(tag_levels = "A") + 
  theme(plot.tag = element_text(size = 12, face='plain', color='black'))


ggsave("manuscript/figures/Figure3.png", limitsize = FALSE, 
       units = c("px"), Figure3, width = 2000, height = 1000, bg='white')  


ggsave("manuscript/figures/Figure3.pdf", limitsize = FALSE, 
       units = c("px"), Figure3, width = 2000, height = 1000) 









panel_bri_3d <- ggdraw() + draw_image(readPNG("manuscript/pictures/bridge.png"))
panel_bri_mito <- ggdraw() + draw_image(readPNG("manuscript/pictures/mito_pos_bridge_text.png"))
panel_bri_syn <- ggdraw() + draw_image(readPNG("manuscript/pictures/SSN_prepost_synapse_bridge.png"))
panel_grav_matrix <- ggdraw() + draw_image(readPNG("manuscript/pictures/mech_girdle_chaeMech_syn_matrix.png"))
panel_grav_graph_bridge <- ggdraw() + draw_image(readPNG("manuscript/pictures/gravity_neural_circuit_graph_bridge.png"))

layout <- "
AAAABBB
#######
CCCDDEE
"

Figure4 <- panel_bri_3d + panel_bri_mito + 
  panel_bri_syn + panel_grav_matrix + panel_grav_graph_bridge +
  plot_layout(design = layout,
              heights = c(1, 0.1, 1.2),
              widths = c(1, 1, 1, 1, 1, 1, 0.5)) + 
  plot_annotation(tag_levels = "A") + 
  theme(plot.tag = element_text(size = 12, face='plain', color='black'))


ggsave("manuscript/figures/Figure4.png", limitsize = FALSE, 
       units = c("px"), Figure4, width = 2400, height = 1000, bg='white')  


ggsave("manuscript/figures/Figure4.pdf", limitsize = FALSE, 
       units = c("px"), Figure4, width = 2400, height = 1000) 




