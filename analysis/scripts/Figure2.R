# Code to generate Figure 2 of the Jokura et al 2024 Ctenophore apical organ connectome paper

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

SSN_Q1Q2 <- read_smooth_neuron(get_skids_with_annot(pid = 35, c("celltype:SSN", "Q1Q2")))
SSN_Q3Q4 <- read_smooth_neuron(get_skids_with_annot(pid = 35, c("celltype:SSN", "Q3Q4")))
SSN_Q1Q2Q3Q4 <- read_smooth_neuron(get_skids_with_annot(pid = 35, c("celltype:SSN", "Q1Q2Q3Q4")))



with_soma <- read_smooth_neuron("with_soma")

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



# plot SSN Q1Q2 & Q3Q4 neuron -----------------------------------------------



close3d()
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  #defines the two scenes
#define the size of the rgl window, the view and zoom
par3d(windowRect = c(0, 0, 1200, 350))
#par3d(windowRect = c(0, 0, 2400, 700))

#plot aboral view
plot3d(SSN_Q1Q2,
       soma = T, lwd = 1, add = T, alpha = 1, col = Okabe_Ito[6],
       WithConnectors = F, WithNodes = F)

plot3d(SSN_Q3Q4,
       soma = T, lwd = 1, add = T, alpha = 1, col = Okabe_Ito[7],
       WithConnectors = F, WithNodes = F)

#for (i in 1:length(all_celltypes)) {
#  print(i)
#  plot3d(
#    all_celltypes[[i]], soma = TRUE, lwd = 0.5, add = TRUE, 
#    alpha = 0.05, col = Okabe_Ito[8]
#  )
#}

plot3d(with_soma,
       soma = T, lwd = 0.5, add = T, alpha = 0.025, col = Okabe_Ito[8],
       WithConnectors = F, WithNodes = F)

plot3d(outline,
       add = T, alpha = 0.05, col = "grey50"
)

#texts3d(58000,51000,5000, "INRGW", cex = 3, col = "#56B4E9")

#aboral view
nview3d("anterior", extramat = rotationMatrix(1.05, 250, -200, 1000))
#rgl.snapshot("manuscript/pictures/SSN_aboral_view.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Sagittal plane
plot3d(SSN_Q1Q2,
       soma = T, lwd = 1, add = T, alpha = 1, col = Okabe_Ito[6],
       WithConnectors = F, WithNodes = F)

plot3d(SSN_Q3Q4,
       soma = T, lwd = 1, add = T, alpha = 1, col = Okabe_Ito[7],
       WithConnectors = F, WithNodes = F)

#for (i in 1:length(all_celltypes)) {
#  print(i)
#  plot3d(
#    all_celltypes[[i]], soma = TRUE, lwd = 0.5, add = TRUE, 
#    alpha = 0.05, col = Okabe_Ito[8]
#  )
#}

plot3d(with_soma,
       soma = T, lwd = 0.5, add = T, alpha = 0.025, col = Okabe_Ito[8],
       WithConnectors = F, WithNodes = F)

plot3d(outline,
       add = T, alpha = 0.05, col = "grey50"
)

#texts3d(58000,51000,5000, "INRGW", cex = 3, col = "#56B4E9")

#lateral view of Sagittal plane
nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
#rgl.snapshot("manuscript/pictures/SSN_sagittal_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Tentacular plane
plot3d(SSN_Q1Q2,
       soma = T, lwd = 1, add = T, alpha = 1, col = Okabe_Ito[6],
       WithConnectors = F, WithNodes = F)

plot3d(SSN_Q3Q4,
       soma = T, lwd = 1, add = T, alpha = 1, col = Okabe_Ito[7],
       WithConnectors = F, WithNodes = F)

#for (i in 1:length(all_celltypes)) {
#  print(i)
#  plot3d(
#    all_celltypes[[i]], soma = TRUE, lwd = 0.5, add = TRUE, 
#    alpha = 0.05, col = Okabe_Ito[8]
#  )
#}

plot3d(with_soma,
       soma = T, lwd = 0.5, add = T, alpha = 0.025, col = Okabe_Ito[8],
       WithConnectors = F, WithNodes = F)

plot3d(outline,
       add = T, alpha = 0.05, col = "grey50"
)


#lateral view of Tentacular plane
nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))
#rgl.snapshot("manuscript/pictures/SSN_tentacular_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)


#make a snapshot to the working directory
rgl.snapshot("manuscript/pictures/SSN_Q12_Q34.png")


close3d()

# plot SSN Q1Q2Q3Q4 neuron -----------------------------------------------


close3d()
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  #defines the two scenes
#define the size of the rgl window, the view and zoom
par3d(windowRect = c(0, 0, 1200, 350))
#par3d(windowRect = c(0, 0, 2400, 700))

#plot aboral view
plot3d(SSN_Q1Q2Q3Q4,
       soma = T, lwd = 1, add = T, alpha = 1, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)

#for (i in 1:length(all_celltypes)) {
#  print(i)
#  plot3d(
#    all_celltypes[[i]], soma = TRUE, lwd = 0.5, add = TRUE, 
#    alpha = 0.05, col = Okabe_Ito[8]
#  )
#}

plot3d(with_soma,
       soma = T, lwd = 0.5, add = T, alpha = 0.025, col = Okabe_Ito[8],
       WithConnectors = F, WithNodes = F)

plot3d(outline,
       add = T, alpha = 0.05, col = "grey50"
)

#texts3d(58000,51000,5000, "INRGW", cex = 3, col = "#56B4E9")

#aboral view
nview3d("anterior", extramat = rotationMatrix(1.05, 250, -200, 1000))
#rgl.snapshot("manuscript/pictures/SSN_aboral_view.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Sagittal plane
plot3d(SSN_Q1Q2Q3Q4,
       soma = T, lwd = 1, add = T, alpha = 1, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)

#for (i in 1:length(all_celltypes)) {
#  print(i)
#  plot3d(
#    all_celltypes[[i]], soma = TRUE, lwd = 0.5, add = TRUE, 
#    alpha = 0.05, col = Okabe_Ito[8]
#  )
#}

plot3d(with_soma,
       soma = T, lwd = 0.5, add = T, alpha = 0.025, col = Okabe_Ito[8],
       WithConnectors = F, WithNodes = F)

plot3d(outline,
       add = T, alpha = 0.05, col = "grey50"
)

#texts3d(58000,51000,5000, "INRGW", cex = 3, col = "#56B4E9")

#lateral view of Sagittal plane
nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
#rgl.snapshot("manuscript/pictures/SSN_sagittal_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Tentacular plane
plot3d(SSN_Q1Q2Q3Q4,
       soma = T, lwd = 1, add = T, alpha = 1, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)

#for (i in 1:length(all_celltypes)) {
#  print(i)
#  plot3d(
#    all_celltypes[[i]], soma = TRUE, lwd = 0.5, add = TRUE, 
#    alpha = 0.05, col = Okabe_Ito[8]
#  )
#}

plot3d(with_soma,
       soma = T, lwd = 0.5, add = T, alpha = 0.025, col = Okabe_Ito[8],
       WithConnectors = F, WithNodes = F)

plot3d(outline,
       add = T, alpha = 0.05, col = "grey50"
)


#lateral view of Tentacular plane
nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))
#rgl.snapshot("manuscript/pictures/SSN_tentacular_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)


#make a snapshot to the working directory
rgl.snapshot("manuscript/pictures/SSN_Q1234.png")


close3d()




# retrieve connectors ----------------------------------------------------------

conn_syn <- connectors(SSN)
presyn_syn <- subset(conn_syn, prepost == 0)
postsyn_syn <- subset(conn_syn, prepost == 1)


# plot SSN synapses ---------------------------------------------------------------------


close3d()
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  #defines the two scenes
#define the size of the rgl window, the view and zoom
par3d(windowRect = c(0, 0, 1200, 350))
#par3d(windowRect = c(0, 0, 2400, 700))

#plot aboral view
plot3d(SSN_Q1Q2,
       soma = FALSE, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[6],
       WithConnectors = F, WithNodes = F)

plot3d(SSN_Q3Q4,
       soma = FALSE, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[7],
       WithConnectors = F, WithNodes = F)

plot3d(SSN_Q1Q2Q3Q4,
       soma = FALSE, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)

#for (i in 1:length(all_celltypes)) {
#  print(i)
#  plot3d(
#    all_celltypes[[i]], soma = TRUE, lwd = 0.5, add = TRUE, 
#    alpha = 0.05, col = Okabe_Ito[8]
#  )
#}

#plot3d(with_soma,
#       soma = T, lwd = 0.5, add = T, alpha = 0.025, col = Okabe_Ito[8],
#       WithConnectors = F, WithNodes = F)

plot3d(outline,
       add = T, alpha = 0.05, col = "grey50"
)

# plot presynapses
plot3d(
  presyn_syn$x, 
  presyn_syn$y, 
  presyn_syn$z, 
  size = 0.6, alpha = 0.5, col = "magenta2", 
  add = TRUE,
  point_antialias = TRUE,
  type = "s"
)

# plot postsynapses
plot3d(
  postsyn_syn$x, 
  postsyn_syn$y, 
  postsyn_syn$z, 
  size = 0.6, alpha = 0.5, col = "cyan2", 
  add = TRUE,
  point_antialias = TRUE,
  type = "s"
)


#texts3d(58000,51000,5000, "INRGW", cex = 3, col = "#56B4E9")

#aboral view
nview3d("anterior", extramat = rotationMatrix(1.05, 250, -200, 1000))
#rgl.snapshot("manuscript/pictures/balancer_aboral_view.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)




#plot lateral view of Sagittal plane
plot3d(SSN_Q1Q2,
       soma = FALSE, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[6],
       WithConnectors = F, WithNodes = F)

plot3d(SSN_Q3Q4,
       soma = FALSE, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[7],
       WithConnectors = F, WithNodes = F)

plot3d(SSN_Q1Q2Q3Q4,
       soma = FALSE, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)

#for (i in 1:length(all_celltypes)) {
#  print(i)
#  plot3d(
#    all_celltypes[[i]], soma = TRUE, lwd = 0.5, add = TRUE, 
#    alpha = 0.05, col = Okabe_Ito[8]
#  )
#}

#plot3d(with_soma,
#       soma = T, lwd = 0.5, add = T, alpha = 0.025, col = Okabe_Ito[8],
#       WithConnectors = F, WithNodes = F)

plot3d(outline,
       add = T, alpha = 0.05, col = "grey50"
)

# plot presynapses
plot3d(
  presyn_syn$x, 
  presyn_syn$y, 
  presyn_syn$z, 
  size = 0.6, alpha = 0.5, col = "magenta2", 
  add = TRUE,
  point_antialias = TRUE,
  type = "s"
)

# plot postsynapses
plot3d(
  postsyn_syn$x, 
  postsyn_syn$y, 
  postsyn_syn$z, 
  size = 0.6, alpha = 0.5, col = "cyan2", 
  add = TRUE,
  point_antialias = TRUE,
  type = "s"
)



nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
#rgl.snapshot("manuscript/pictures/balancer_sagittal_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Tentacular plane
plot3d(SSN_Q1Q2,
       soma = FALSE, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[6],
       WithConnectors = F, WithNodes = F)

plot3d(SSN_Q3Q4,
       soma = FALSE, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[7],
       WithConnectors = F, WithNodes = F)

plot3d(SSN_Q1Q2Q3Q4,
       soma = FALSE, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)

#for (i in 1:length(all_celltypes)) {
#  print(i)
#  plot3d(
#    all_celltypes[[i]], soma = TRUE, lwd = 0.5, add = TRUE, 
#    alpha = 0.05, col = Okabe_Ito[8]
#  )
#}

#plot3d(with_soma,
#       soma = T, lwd = 0.5, add = T, alpha = 0.025, col = Okabe_Ito[8],
#       WithConnectors = F, WithNodes = F)

plot3d(outline,
       add = T, alpha = 0.05, col = "grey50"
)

# plot presynapses
plot3d(
  presyn_syn$x, 
  presyn_syn$y, 
  presyn_syn$z, 
  size = 0.6, alpha = 0.5, col = "magenta2", 
  add = TRUE,
  point_antialias = TRUE,
  type = "s"
)

# plot postsynapses
plot3d(
  postsyn_syn$x, 
  postsyn_syn$y, 
  postsyn_syn$z, 
  size = 0.6, alpha = 0.5, col = "cyan2", 
  add = TRUE,
  point_antialias = TRUE,
  type = "s"
)


nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))
#rgl.snapshot("manuscript/pictures/balancer_tentacular_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)


#make a snapshot to the working directory
rgl.snapshot("manuscript/pictures/SSN_prepost_synapse.png")

close3d()




# assemble figure -------------------------------------------------------------

panel_mito_in_syn_assemble <- ggdraw() + draw_image(readPNG("manuscript/pictures/mito_in_syn_assemble.png"))
panel_mito_pos_SSN_anterior_trim <- ggdraw() + draw_image(readPNG("manuscript/pictures/mito_pos_SSN_anterior_trim.png"))
panel_mito_pos_SSN_left <- ggdraw() + draw_image(readPNG("manuscript/pictures/mito_pos_SSN_left.png"))


panel_SSN_Q1234 <- ggdraw() + draw_image(readPNG("manuscript/pictures/SSN_Q1234_from_CATMAID.png"))
panel_SSN_Q12_Q34 <- ggdraw() + draw_image(readPNG("manuscript/pictures/SSN_Q12_Q34_from_CATMAID.png"))

panel_synapse <- ggdraw() + draw_image(readPNG("manuscript/pictures/Figure_mito_syn_ves_syn.png"))

panel_SSN_prepost_synapse <- ggdraw() + draw_image(readPNG("manuscript/pictures/SSN_prepost_synapse.png"))

layout <- "
AABCDD
######
EEEFFF
######
GGGGHH
"

Figure2 <- panel_mito_in_syn_assemble + panel_mito_pos_SSN_anterior_trim + panel_mito_pos_SSN_anterior_trim + panel_mito_pos_SSN_left +
  panel_SSN_Q1234 + panel_SSN_Q12_Q34 +
  panel_SSN_prepost_synapse + panel_synapse + 
  plot_layout(design = layout,
              heights = c(1.25, 0.1, 1, 0.1, 1),
              widths = c(1, 1, 1, 1, 1, 1)) + 
  patchwork::plot_annotation(tag_levels = "A") &  
  ggplot2::theme(plot.tag = element_text(size = 12, 
                                         face='plain', color='black'))


ggsave("manuscript/figures/Figure2.png", limitsize = FALSE, 
       units = c("px"), Figure2, width = 2800, height = 2000, bg='white')  


ggsave("manuscript/figures/Figure2.pdf", limitsize = FALSE, 
       units = c("px"), Figure2, width = 2800, height = 2000) 

