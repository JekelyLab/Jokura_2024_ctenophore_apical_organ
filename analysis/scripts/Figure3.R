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


Q1 <- read_smooth_neuron("Q1")
Q2 <- read_smooth_neuron("Q2")
Q3 <- read_smooth_neuron("Q3")
Q4 <- read_smooth_neuron("Q4")


balancer_Q1 <- read_smooth_neuron(get_skids_with_annot(pid = 35, c("celltype:balancer", "Q1")))
balancer_Q2 <- read_smooth_neuron(get_skids_with_annot(pid = 35, c("celltype:balancer", "Q2")))
balancer_Q3 <- read_smooth_neuron(get_skids_with_annot(pid = 35, c("celltype:balancer", "Q3")))
balancer_Q4 <- read_smooth_neuron(get_skids_with_annot(pid = 35, c("celltype:balancer", "Q4")))

bridge_Q1Q2 <- read_smooth_neuron(get_skids_with_annot(pid = 35, c("celltype:bridge", "Q1Q2")))
bridge_Q3Q4 <- read_smooth_neuron(get_skids_with_annot(pid = 35, c("celltype:bridge", "Q3Q4")))

SSN_Q1Q2 <- read_smooth_neuron(get_skids_with_annot(pid = 35, c("celltype:SSN", "Q1Q2")))
SSN_Q3Q4 <- read_smooth_neuron(get_skids_with_annot(pid = 35, c("celltype:SSN", "Q3Q4")))
SSN_Q1Q2Q3Q4 <- read_smooth_neuron(get_skids_with_annot(pid = 35, c("celltype:SSN", "Q1Q2Q3Q4")))


# plot balancer -------------------------------------------------------------

close3d()

# 3d plotting of cells
nopen3d()
mfrow3d(1, 3) 
# define the size of the rgl window, the view and zoom
par3d(windowRect = c(0, 0, 1200, 350))
#par3d(windowRect = c(0, 0, 2400, 700))

# plot aboral view

plot3d(
  balancer_Q1,
  soma = T, lwd = 1.5, add = T,
  alpha = 0.5, col = Okabe_Ito[1]
)
plot3d(
  balancer_Q2,
  soma = T, lwd = 1.5, add = T,
  alpha = 0.5, col = Okabe_Ito[2]
)
plot3d(
  balancer_Q3,
  soma = T, lwd = 1.5, add = T,
  alpha = 0.5, col = Okabe_Ito[6]
)
plot3d(
  balancer_Q4,
  soma = T, lwd = 1.5, add = T,
  alpha = 0.5, col = Okabe_Ito[7]
)

plot3d(
  with_soma,
  soma = T, lwd = 1, add = T,
  alpha = 0.01, col = Okabe_Ito[8]
)

plot3d(
  outline,
  add = T, alpha = 0.05, col = "grey50"
)

#aboral view
aboral()
par3d(zoom = 0.61)

## y-axis clip
#clipplanes3d(1, 0, 0, -11500)
## x-axis clip
#clipplanes3d(0, 1, 0, -24000)


# move to next panel in rgl window
next3d(clear = F)

# plot lateral view of Sagittal plane

plot3d(
  balancer_Q1,
  soma = T, lwd = 1.5, add = T,
  alpha = 0.5, col = Okabe_Ito[1]
)
plot3d(
  balancer_Q2,
  soma = T, lwd = 1.5, add = T,
  alpha = 0.5, col = Okabe_Ito[2]
)
plot3d(
  balancer_Q3,
  soma = T, lwd = 1.5, add = T,
  alpha = 0.5, col = Okabe_Ito[6]
)
plot3d(
  balancer_Q4,
  soma = T, lwd = 1.5, add = T,
  alpha = 0.5, col = Okabe_Ito[7]
)

plot3d(
  with_soma,
  soma = T, lwd = 1, add = T,
  alpha = 0.01, col = Okabe_Ito[8]
)

plot3d(
  outline,
  add = T, alpha = 0.05, col = "grey50"
)

sagittal()
par3d(zoom = 0.61)

# y-axis clip
#clipplanes3d(1, 0, 0, -11500)
# x-axis clip
#clipplanes3d(0, 1, 0, -24000)

# move to next panel in rgl window
next3d(clear = F)

# plot lateral view of Tentacular plane

plot3d(
  balancer_Q1,
  soma = T, lwd = 1.5, add = T,
  alpha = 0.5, col = Okabe_Ito[1]
)
plot3d(
  balancer_Q2,
  soma = T, lwd = 1.5, add = T,
  alpha = 0.5, col = Okabe_Ito[2]
)
plot3d(
  balancer_Q3,
  soma = T, lwd = 1.5, add = T,
  alpha = 0.5, col = Okabe_Ito[6]
)
plot3d(
  balancer_Q4,
  soma = T, lwd = 1.5, add = T,
  alpha = 0.5, col = Okabe_Ito[7]
)

plot3d(
  with_soma,
  soma = T, lwd = 1, add = T,
  alpha = 0.01, col = Okabe_Ito[8]
)

plot3d(
  outline,
  add = T, alpha = 0.05, col = "grey50"
)

tentacular()
par3d(zoom = 0.61)

# y-axis clip
#clipplanes3d(1, 0, 0, -11500)
# x-axis clip
#clipplanes3d(0, 1, 0, -24000)

# make a snapshot to the working directory
rgl.snapshot("manuscript/pictures/balancer.png")
close3d()





# plot bridge -----------------------------------------------

close3d()

# 3d plotting of cells
nopen3d()
mfrow3d(1, 3) 
# define the size of the rgl window, the view and zoom
par3d(windowRect = c(0, 0, 1200, 350))
#par3d(windowRect = c(0, 0, 2400, 700))

# plot aboral view

plot3d(
  bridge_Q1Q2,
  soma = T, lwd = 1.5, add = T,
  alpha = 0.75, col = Okabe_Ito[2]
)
plot3d(
  bridge_Q3Q4,
  soma = T, lwd = 1.5, add = T,
  alpha = 0.75, col = Okabe_Ito[7]
)

plot3d(
  with_soma,
  soma = T, lwd = 1, add = T,
  alpha = 0.01, col = Okabe_Ito[8]
)

plot3d(
  outline,
  add = T, alpha = 0.05, col = "grey50"
)

#aboral view
aboral()
par3d(zoom = 0.61)

## y-axis clip
#clipplanes3d(1, 0, 0, -11500)
## x-axis clip
#clipplanes3d(0, 1, 0, -24000)


# move to next panel in rgl window
next3d(clear = F)

# plot lateral view of Sagittal plane

plot3d(
  bridge_Q1Q2,
  soma = T, lwd = 1.5, add = T,
  alpha = 0.75, col = Okabe_Ito[2]
)
plot3d(
  bridge_Q3Q4,
  soma = T, lwd = 1.5, add = T,
  alpha = 0.75, col = Okabe_Ito[7]
)

plot3d(
  with_soma,
  soma = T, lwd = 1, add = T,
  alpha = 0.01, col = Okabe_Ito[8]
)

plot3d(
  outline,
  add = T, alpha = 0.05, col = "grey50"
)

sagittal()
par3d(zoom = 0.61)

# y-axis clip
#clipplanes3d(1, 0, 0, -11500)
# x-axis clip
#clipplanes3d(0, 1, 0, -24000)

# move to next panel in rgl window
next3d(clear = F)

# plot lateral view of Tentacular plane

plot3d(
  bridge_Q1Q2,
  soma = T, lwd = 1.5, add = T,
  alpha = 0.75, col = Okabe_Ito[2]
)
plot3d(
  bridge_Q3Q4,
  soma = T, lwd = 1.5, add = T,
  alpha = 0.75, col = Okabe_Ito[7]
)

plot3d(
  with_soma,
  soma = T, lwd = 1, add = T,
  alpha = 0.01, col = Okabe_Ito[8]
)

plot3d(
  outline,
  add = T, alpha = 0.05, col = "grey50"
)

tentacular()
par3d(zoom = 0.61)

# y-axis clip
#clipplanes3d(1, 0, 0, -11500)
# x-axis clip
#clipplanes3d(0, 1, 0, -24000)


#make a snapshot to the working directory
rgl.snapshot("manuscript/pictures/bridge.png")


close3d()


# load of mitochondrial location information-------------------------------

mito_vesicle_info <- read.csv("analysis/data/mito_vesicle_info.csv")

# 3D plot mitochondria positions in bridge ----------------------------------------

close3d()
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  #defines the two scenes
#define the size of the rgl window, the view and zoom
par3d(windowRect = c(0, 0, 1200, 350))
#par3d(windowRect = c(0, 0, 2400, 700))

#plot aboral view

plot3d(bridge_Q1Q2,
       soma = FALSE, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[2],
       WithConnectors = F, WithNodes = F)

plot3d(bridge_Q3Q4,
       soma = FALSE, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[7],
       WithConnectors = F, WithNodes = F)

plot3d(outline,
       add = T, alpha = 0.05, col = "grey50"
)

# plot mitochondria

pos_ves_syn <- mito_vesicle_info |>
  filter(celltype == "bridge") |>
  filter(mito_type == "vesicles_syn") |>
  select(x, y, z)

plot3d(pos_ves_syn, 
       add = TRUE,
       col = "green3",
       size = 0.9, 
       alpha = 1,
       point_antialias = TRUE,
       type = "s"
)

pos_ves_no_syn <- mito_vesicle_info |>
  filter(celltype == "bridge") |>
  filter(mito_type != "vesicles_syn") |>
  select(x, y, z)

plot3d(pos_ves_no_syn, 
       add = TRUE,
       col = "black",
       size = 0.9, 
       alpha = 0.37,
       point_antialias = TRUE,
       type = "s"
)

#aboral view
aboral()

par3d(zoom=0.61)


#move to next panel in rgl window
next3d(clear=F)

#plot sagittal view

plot3d(bridge_Q1Q2,
       soma = FALSE, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[2],
       WithConnectors = F, WithNodes = F)

plot3d(bridge_Q3Q4,
       soma = FALSE, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[7],
       WithConnectors = F, WithNodes = F)

plot3d(outline,
       add = T, alpha = 0.05, col = "grey50"
)

# plot mitochondria

pos_ves_syn <- mito_vesicle_info |>
  filter(celltype == "bridge") |>
  filter(mito_type == "vesicles_syn") |>
  select(x, y, z)

plot3d(pos_ves_syn, 
       add = TRUE,
       col = "green3",
       size = 0.9, 
       alpha = 1,
       point_antialias = TRUE,
       type = "s"
)

pos_ves_no_syn <- mito_vesicle_info |>
  filter(celltype == "bridge") |>
  filter(mito_type != "vesicles_syn") |>
  select(x, y, z)

plot3d(pos_ves_no_syn, 
       add = TRUE,
       col = "black",
       size = 0.9, 
       alpha = 0.37,
       point_antialias = TRUE,
       type = "s"
)

#sagittal view
sagittal()

par3d(zoom=0.61)


#move to next panel in rgl window
next3d(clear=F)

#plot tentacular view

plot3d(bridge_Q1Q2,
       soma = FALSE, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[2],
       WithConnectors = F, WithNodes = F)

plot3d(bridge_Q3Q4,
       soma = FALSE, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[7],
       WithConnectors = F, WithNodes = F)

plot3d(outline,
       add = T, alpha = 0.05, col = "grey50"
)

# plot mitochondria

pos_ves_syn <- mito_vesicle_info |>
  filter(celltype == "bridge") |>
  filter(mito_type == "vesicles_syn") |>
  select(x, y, z)

plot3d(pos_ves_syn, 
       add = TRUE,
       col = "green3",
       size = 0.9, 
       alpha = 1,
       point_antialias = TRUE,
       type = "s"
)

pos_ves_no_syn <- mito_vesicle_info |>
  filter(celltype == "bridge") |>
  filter(mito_type != "vesicles_syn") |>
  select(x, y, z)

plot3d(pos_ves_no_syn, 
       add = TRUE,
       col = "black",
       size = 0.9, 
       alpha = 0.37,
       point_antialias = TRUE,
       type = "s"
)

#tentacular view
tentacular()

par3d(zoom=0.61)


#make a snapshot to the working directory
rgl.snapshot("manuscript/pictures/mito_pos_bridge.png")

close3d()




# retrieve connectors ----------------------------------------------------------

balancer_conn <- connectors(balancer)
presyn_balancer <- subset(balancer_conn, prepost == 0)
postsyn_balancer <- subset(balancer_conn, prepost == 1)

bridge_conn <- connectors(bridge)
presyn_bridge <- subset(bridge_conn, prepost == 0)
postsyn_bridge <- subset(bridge_conn, prepost == 1)


# 3D plot of synaptic input from SSN to balancer ---------------------------------------------------------------------

close3d()
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  
#define the size of the rgl window, the view and zoom
par3d(windowRect = c(0, 0, 1200, 350))
#par3d(windowRect = c(0, 0, 2400, 700))

#plot aboral view

plot3d(balancer_Q1,
       soma = FALSE, lwd = 1, add = T, alpha = 0.1, col = Okabe_Ito[8],
       WithConnectors = F, WithNodes = F)

plot3d(balancer_Q2,
       soma = FALSE, lwd = 1, add = T, alpha = 0.1, col = Okabe_Ito[8],
       WithConnectors = F, WithNodes = F)

plot3d(balancer_Q3,
       soma = FALSE, lwd = 1, add = T, alpha = 0.1, col = Okabe_Ito[8],
       WithConnectors = F, WithNodes = F)

plot3d(balancer_Q4,
       soma = FALSE, lwd = 1, add = T, alpha = 0.1, col = Okabe_Ito[8],
       WithConnectors = F, WithNodes = F)

plot3d(SSN_Q1Q2Q3Q4,
       soma = FALSE, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)

plot3d(SSN_Q1Q2,
       soma = FALSE, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[6],
       WithConnectors = F, WithNodes = F)

plot3d(SSN_Q3Q4,
       soma = FALSE, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[7],
       WithConnectors = F, WithNodes = F)



plot3d(outline,
       add = T, alpha = 0.05, col = "grey50"
)

# input from SSN to balancer
plot3d(
  postsyn_balancer$x, 
  postsyn_balancer$y, 
  postsyn_balancer$z, 
  size = 0.6, alpha = 1, col = "magenta2", 
  add = TRUE,
  point_antialias = TRUE,
  type = "s"
)

#aboral view
aboral()
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)


#plot lateral view of Sagittal plane

plot3d(balancer_Q1,
       soma = FALSE, lwd = 1, add = T, alpha = 0.1, col = Okabe_Ito[8],
       WithConnectors = F, WithNodes = F)

plot3d(balancer_Q2,
       soma = FALSE, lwd = 1, add = T, alpha = 0.1, col = Okabe_Ito[8],
       WithConnectors = F, WithNodes = F)

plot3d(balancer_Q3,
       soma = FALSE, lwd = 1, add = T, alpha = 0.1, col = Okabe_Ito[8],
       WithConnectors = F, WithNodes = F)

plot3d(balancer_Q4,
       soma = FALSE, lwd = 1, add = T, alpha = 0.1, col = Okabe_Ito[8],
       WithConnectors = F, WithNodes = F)

plot3d(SSN_Q1Q2Q3Q4,
       soma = FALSE, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)

plot3d(SSN_Q1Q2,
       soma = FALSE, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[6],
       WithConnectors = F, WithNodes = F)

plot3d(SSN_Q3Q4,
       soma = FALSE, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[7],
       WithConnectors = F, WithNodes = F)



plot3d(outline,
       add = T, alpha = 0.05, col = "grey50"
)

# input from SSN to balancer
plot3d(
  postsyn_balancer$x, 
  postsyn_balancer$y, 
  postsyn_balancer$z, 
  size = 0.6, alpha = 1, col = "magenta2", 
  add = TRUE,
  point_antialias = TRUE,
  type = "s"
)

sagittal()
par3d(zoom=0.61)
next3d(clear=F)



#plot lateral view of Tentacular plane
plot3d(balancer_Q1,
       soma = FALSE, lwd = 1, add = T, alpha = 0.1, col = Okabe_Ito[8],
       WithConnectors = F, WithNodes = F)

plot3d(balancer_Q2,
       soma = FALSE, lwd = 1, add = T, alpha = 0.1, col = Okabe_Ito[8],
       WithConnectors = F, WithNodes = F)

plot3d(balancer_Q3,
       soma = FALSE, lwd = 1, add = T, alpha = 0.1, col = Okabe_Ito[8],
       WithConnectors = F, WithNodes = F)

plot3d(balancer_Q4,
       soma = FALSE, lwd = 1, add = T, alpha = 0.1, col = Okabe_Ito[8],
       WithConnectors = F, WithNodes = F)

plot3d(SSN_Q1Q2Q3Q4,
       soma = FALSE, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)

plot3d(SSN_Q1Q2,
       soma = FALSE, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[6],
       WithConnectors = F, WithNodes = F)

plot3d(SSN_Q3Q4,
       soma = FALSE, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[7],
       WithConnectors = F, WithNodes = F)



plot3d(outline,
       add = T, alpha = 0.05, col = "grey50"
)

# input from SSN to balancer
plot3d(
  postsyn_balancer$x, 
  postsyn_balancer$y, 
  postsyn_balancer$z, 
  size = 0.6, alpha = 1, col = "magenta2", 
  add = TRUE,
  point_antialias = TRUE,
  type = "s"
)

tentacular()
par3d(zoom=0.61)

#make a snapshot to the working directory
rgl.snapshot("manuscript/pictures/SSN_prepost_synapse_balancer.png")

close3d()





# 3D plot of reciprocal synaptic inputs between SSN and bridge-----------------------------

close3d()
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  #defines the two scenes
#define the size of the rgl window, the view and zoom
par3d(windowRect = c(0, 0, 1200, 350))
#par3d(windowRect = c(0, 0, 2400, 700))

#plot aboral view
plot3d(bridge_Q1Q2,
       soma = FALSE, lwd = 1, add = T, alpha = 0.2, col = Okabe_Ito[8],
       WithConnectors = F, WithNodes = F)

plot3d(bridge_Q3Q4,
       soma = FALSE, lwd = 1, add = T, alpha = 0.2, col = Okabe_Ito[8],
       WithConnectors = F, WithNodes = F)

plot3d(SSN_Q1Q2Q3Q4,
       soma = FALSE, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)

plot3d(SSN_Q1Q2,
       soma = FALSE, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[6],
       WithConnectors = F, WithNodes = F)

plot3d(SSN_Q3Q4,
       soma = FALSE, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[7],
       WithConnectors = F, WithNodes = F)

plot3d(outline,
       add = T, alpha = 0.05, col = "grey50"
)

# inputs from SSN to bridge
plot3d(
  postsyn_bridge$x, 
  postsyn_bridge$y, 
  postsyn_bridge$z, 
  size = 0.8, alpha = 0.7, col = "magenta2", 
  add = TRUE,
  point_antialias = TRUE,
  type = "s"
)

# inputs from bridge to SNN
plot3d(
  presyn_bridge$x, 
  presyn_bridge$y, 
  presyn_bridge$z, 
  size = 0.8, alpha = 1, col = "cyan2", 
  add = TRUE,
  point_antialias = TRUE,
  type = "s"
)

aboral()
par3d(zoom=0.61)
next3d(clear=F)




#plot lateral view of Sagittal plane
plot3d(bridge_Q1Q2,
       soma = FALSE, lwd = 1, add = T, alpha = 0.2, col = Okabe_Ito[8],
       WithConnectors = F, WithNodes = F)

plot3d(bridge_Q3Q4,
       soma = FALSE, lwd = 1, add = T, alpha = 0.2, col = Okabe_Ito[8],
       WithConnectors = F, WithNodes = F)

plot3d(SSN_Q1Q2Q3Q4,
       soma = FALSE, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)

plot3d(SSN_Q1Q2,
       soma = FALSE, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[6],
       WithConnectors = F, WithNodes = F)

plot3d(SSN_Q3Q4,
       soma = FALSE, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[7],
       WithConnectors = F, WithNodes = F)

plot3d(outline,
       add = T, alpha = 0.05, col = "grey50"
)

# inputs from SSN to bridge
plot3d(
  postsyn_bridge$x, 
  postsyn_bridge$y, 
  postsyn_bridge$z, 
  size = 0.8, alpha = 0.7, col = "magenta2", 
  add = TRUE,
  point_antialias = TRUE,
  type = "s"
)

# inputs from bridge to SNN
plot3d(
  presyn_bridge$x, 
  presyn_bridge$y, 
  presyn_bridge$z, 
  size = 0.8, alpha = 1, col = "cyan2", 
  add = TRUE,
  point_antialias = TRUE,
  type = "s"
)

sagittal()
par3d(zoom=0.61)
next3d(clear=F)




#plot lateral view of Tentacular plane
plot3d(bridge_Q1Q2,
       soma = FALSE, lwd = 1, add = T, alpha = 0.2, col = Okabe_Ito[8],
       WithConnectors = F, WithNodes = F)

plot3d(bridge_Q3Q4,
       soma = FALSE, lwd = 1, add = T, alpha = 0.2, col = Okabe_Ito[8],
       WithConnectors = F, WithNodes = F)

plot3d(SSN_Q1Q2Q3Q4,
       soma = FALSE, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)

plot3d(SSN_Q1Q2,
       soma = FALSE, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[6],
       WithConnectors = F, WithNodes = F)

plot3d(SSN_Q3Q4,
       soma = FALSE, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[7],
       WithConnectors = F, WithNodes = F)

plot3d(outline,
       add = T, alpha = 0.05, col = "grey50"
)

# inputs from SSN to bridge
plot3d(
  postsyn_bridge$x, 
  postsyn_bridge$y, 
  postsyn_bridge$z, 
  size = 0.8, alpha = 0.7, col = "magenta2", 
  add = TRUE,
  point_antialias = TRUE,
  type = "s"
)

# inputs from bridge to SNN
plot3d(
  presyn_bridge$x, 
  presyn_bridge$y, 
  presyn_bridge$z, 
  size = 0.8, alpha = 1, col = "cyan2", 
  add = TRUE,
  point_antialias = TRUE,
  type = "s"
)

tentacular()
par3d(zoom=0.61)
next3d(clear=F)


#make a snapshot to the working directory
rgl.snapshot("manuscript/pictures/SSN_prepost_synapse_bridge.png")

close3d()


# assemble figure -------------------------------------------------------------

panel_bal <- ggdraw() + draw_image(readPNG("manuscript/pictures/balancer.png")) +
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
  draw_label("Q1", x = 0.25, y = 0.88, color = Okabe_Ito[1], size = 9, hjust = 0, alpha = 1) +
  draw_label("Q2", x = 0.1, y = 0.78, color = Okabe_Ito[2], size = 9, hjust = 1, alpha = 1) +
  draw_label("Q3", x = 0.075, y = 0.1, color = Okabe_Ito[6], size = 9, hjust = 1, alpha = 1) +
  draw_label("Q4", x = 0.24, y = 0.2, color = Okabe_Ito[7], size = 9, hjust = 0, alpha = 1) +
  draw_label("Balancer cells", x = 0.5, y = 1, color = "black", size = 10, hjust = 0.5)

panel_bri <- ggdraw() + draw_image(readPNG("manuscript/pictures/bridge.png")) +
  draw_label("Q1Q2", x = 0.2, y = 0.88, color = Okabe_Ito[2], size = 9, hjust = 0, alpha = 1) +
  draw_label("Q3Q4", x = 0.2, y = 0.18, color = Okabe_Ito[7], size = 9, hjust = 0, alpha = 1) +
  draw_label("Bridge cells", x = 0.5, y = 1, color = "black", size = 10, hjust = 0.5)

panel_bri_mito <- ggdraw() + draw_image(readPNG("manuscript/pictures/mito_pos_bridge.png")) +
  draw_label("mitochondria with synapses", x = 0.06, y = 0.99, color = "green4", size = 8.5, hjust = 0) +
  draw_label("mitochondria not forming synapses", x = 0.06, y = 0.91, color = "black", size = 8.5, hjust = 0, alpha = 0.7) +
  draw_label("bridge Q1Q2", x = 0.95, y = 1, color = Okabe_Ito[2], size = 7, hjust = 1) +
  draw_label("bridge Q3Q4", x = 0.95, y = 0.92, color = Okabe_Ito[7], size = 7, hjust = 1)

panel_SSN_bal <- ggdraw() + draw_image(readPNG("manuscript/pictures/SSN_prepost_synapse_balancer.png")) +
  draw_label("SSN to balancers", x = 0.06, y = 1, color = "#ff00ff", size = 8, hjust = 0)

panel_SSN_bri <- ggdraw() + draw_image(readPNG("manuscript/pictures/SSN_prepost_synapse_bridge.png")) +
  draw_label("SSN to bridges", x = 0.06, y = 1, color = "#ff00ff", size = 8, hjust = 0) +
  draw_label("bridges to SSN", x = 0.25, y = 1, color = "#00c9ff", size = 8, hjust = 0) +
  draw_label("SSN Q1-4", x = 0.85, y = 1.15, color = Okabe_Ito[5], size = 7, hjust = 0) +
  draw_label("SSN Q1Q2", x = 0.85, y = 1.07, color = Okabe_Ito[6], size = 7, hjust = 0) +
  draw_label("SSN Q3Q4", x = 0.85, y = 0.99, color = Okabe_Ito[7], size = 7, hjust = 0)

panel_matrix <- ggdraw() + draw_image(readPNG("manuscript/pictures/mech_girdle_chaeMech_syn_matrix.png"))
panel_map_graph <- ggdraw() + draw_image(readPNG("manuscript/pictures/map_balancer_bridge_SSN.png"))

layout <- "
AAABBBCCC
#########
DDDDFFFGG
####FFFGG
EEEEFFFGG
"

Figure3 <- panel_bal + panel_bri + panel_bri_mito +
  panel_SSN_bal + panel_SSN_bri + 
  panel_matrix + panel_map_graph +
  plot_layout(design = layout,
              heights = c(1, 0.15, 1, 0.1, 1),
              widths = c(1, 1, 1, 1, 1, 1, 1, 1, 1)) + 
  patchwork::plot_annotation(tag_levels = "A") &  
  ggplot2::theme(plot.tag = element_text(size = 12, 
                                         face='plain', color='black'))


ggsave("manuscript/figures/Figure3.png", limitsize = FALSE, 
       units = c("px"), Figure3, width = 3000, height = 1200, bg='white')  


ggsave("manuscript/figures/Figure3.pdf", limitsize = FALSE, 
       units = c("px"), Figure3, width = 3000, height = 1200) 


