
#Synaptic connection from SSN, bridge cells

# source packages and functions ------------------------------------------------

source("analysis/scripts/packages_and_functions.R")

# load cell type ---------------------------------------------------------------

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

all_celltypes_names <- list("balancer",
                      "bridge",
                      "bristle",
                      "Cgroove_sag",
                      "Cgroove_tag",
                      "dense_vesicle",
                      "dome",
                      "epithelial_floor",
                      "intra_multi_ciliated",
                      "lamellate",
                      "lithocyte",
                      "plumose",
                      "SSN",
                      "monociliated",
                      "biciliated",
                      "multiciliated",
                      "nonciliated")


# retrieve connectors ----------------------------------------------------------

conn_syn <- connectors(SSN)
presyn_syn <- subset(conn_syn, prepost == 0)
postsyn_syn <- subset(conn_syn, prepost == 1)

# plot SSN synapses ---------------------------------------------------------------------

plot_background()
close3d()
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  #defines the two scenes
#define the size of the rgl window, the view and zoom
#par3d(windowRect = c(0, 0, 1200, 350))
par3d(windowRect = c(0, 0, 2400, 700))

#plot aboral view
plot3d(SSN,
       soma = FALSE, lwd = 2, add = TRUE, alpha = 0.5, 
       col = Okabe_Ito[c(1, 5, 8)],
       WithConnectors = F, WithNodes = F)

# plot presynapses
plot3d(
  presyn_syn$x, 
  presyn_syn$y, 
  presyn_syn$z, 
  size = 5, alpha = 1, col = "magenta2", 
  add = TRUE,
  point_antialias = TRUE
)

# plot postsynapses
plot3d(
  postsyn_syn$x, 
  postsyn_syn$y, 
  postsyn_syn$z, 
  size = 5, alpha = 1, col = "cyan2", 
  add = TRUE,
  point_antialias = TRUE
)

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
plot3d(SSN,
       soma = FALSE, lwd = 2, add = TRUE, alpha = 0.5, 
       col = Okabe_Ito[c(1, 5, 8)],
       WithConnectors = F, WithNodes = F)

# plot presynapses
plot3d(
  presyn_syn$x, 
  presyn_syn$y, 
  presyn_syn$z, 
  size = 5, alpha = 1, col = "magenta2", 
  add = TRUE,
  point_antialias = TRUE
)

# plot postsynapses
plot3d(
  postsyn_syn$x, 
  postsyn_syn$y, 
  postsyn_syn$z, 
  size = 5, alpha = 1, col = "cyan2", 
  add = TRUE,
  point_antialias = TRUE
)

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
plot3d(SSN,
       soma = FALSE, lwd = 2, add = TRUE, alpha = 0.5, 
       col = Okabe_Ito[c(1, 5, 8)],
       WithConnectors = F, WithNodes = F)

# plot presynapses
plot3d(
  presyn_syn$x, 
  presyn_syn$y, 
  presyn_syn$z, 
  size = 5, alpha = 1, col = "magenta2", 
  add = TRUE,
  point_antialias = TRUE
)

# plot postsynapses
plot3d(
  postsyn_syn$x, 
  postsyn_syn$y, 
  postsyn_syn$z, 
  size = 5, alpha = 1, col = "cyan2", 
  add = TRUE,
  point_antialias = TRUE
)

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
rgl.snapshot("manuscript/pictures/SSN_prepost_synapse.png")

close3d()

# iterate through cell group neuron lists and get connectivity-------------------
# define empty synapse list with the right dimensions
synapse_list <- c()

for (i in 1:length(all_celltypes)) {
  for (j in 1:length(all_celltypes)) {
    # get connectors between two cell groups
    presyn_skids <- attr(all_celltypes[i][[1]], "df")$skid
    postsyn_skids <- attr(all_celltypes[j][[1]], "df")$skid
    connectivity <- catmaid_get_connectors_between(
      pre = presyn_skids,
      post = postsyn_skids, pid = 35
    )
    # check the number of synapses from group1 -> group2
    N_synapses <- dim(connectivity)[1]
    if(length(connectivity) == 0) {N_synapses = 0}
    # add value to synapse list
    synapse_list <- c(synapse_list, N_synapses)
  }
}
synapse_list
# convert synapse list into a matrix of appropriate dimensions
synapse_matrix <- matrix(
  unlist(synapse_list), byrow = TRUE, 
  nrow = length(all_celltypes)
  )

rownames(synapse_matrix) <- as.character(all_celltypes_names)
colnames(synapse_matrix) <- as.character(all_celltypes_names)

# with the make_graph function of igraph we turn it into a graph (input is the list of edge pairs)
AO_graph <- graph_from_adjacency_matrix(
  synapse_matrix,
  mode = c("directed"),
  weighted = TRUE, diag = TRUE
)

#AO_graph |> as_tbl_graph() 

# calculate node weighted degree -----------------------------------------------
degree=degree(
  AO_graph, v = V(AO_graph), mode = c("all"), 
  loops = TRUE, normalized = FALSE
  )

# retrieve connectors ----------------------------------------------------------

balancer_conn <- connectors(balancer)
presyn_balancer <- subset(balancer_conn, prepost == 0)
postsyn_balancer <- subset(balancer_conn, prepost == 1)

bridge_conn <- connectors(bridge)
presyn_bridge <- subset(bridge_conn, prepost == 0)
postsyn_bridge <- subset(bridge_conn, prepost == 1)


monoC_conn <- connectors(monociliated)
presyn_monoC <- subset(monoC_conn, prepost == 0)
postsyn_monoC <- subset(monoC_conn, prepost == 1)

biC_conn <- connectors(biciliated)
presyn_biC <- subset(biC_conn, prepost == 0)
postsyn_biC <- subset(biC_conn, prepost == 1)

nonC_conn <- connectors(nonciliated)
presyn_nonC <- subset(nonC_conn, prepost == 0)
postsyn_nonC <- subset(nonC_conn, prepost == 1)

# plot balancer ----------------------------------------------------------------

plot_background()
close3d()
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  #defines the two scenes
#define the size of the rgl window, the view and zoom
#par3d(windowRect = c(0, 0, 1200, 350))
par3d(windowRect = c(0, 0, 2400, 700))

#plot aboral view
plot3d(SSN,
       soma = FALSE, lwd = 2, add = TRUE, alpha = 0.5, 
       col = Okabe_Ito[c(1, 5, 8)],
       WithConnectors = F, WithNodes = F)

plot3d(
  balancer, soma = TRUE, color = Okabe_Ito[8], 
  alpha = 0.05, lwd = 2
)

# plot postsynapses
plot3d(
  postsyn_balancer$x, 
  postsyn_balancer$y, 
  postsyn_balancer$z, 
  size = 6, alpha = 1, col = "magenta2", 
  add = T
)

#texts3d(58000,51000,5000, "INRGW", cex = 3, col = "#56B4E9")

#aboral view
nview3d("anterior", extramat = rotationMatrix(1.05, 250, -200, 1000))
#rgl.snapshot("manuscript/pictures/balancer_aboral_view.png")
par3d(zoom=0.67)

#move to next panel in rgl window
next3d(clear=F)




#plot lateral view of Sagittal plane
plot3d(SSN,
       soma = FALSE, lwd = 2, add = TRUE, alpha = 0.5, 
       col = Okabe_Ito[c(1, 5, 8)],
       WithConnectors = F, WithNodes = F)

plot3d(
  balancer, soma = TRUE, color = Okabe_Ito[8], 
  alpha = 0.05, lwd = 2
)

# plot postsynapses
plot3d(
  postsyn_balancer$x, 
  postsyn_balancer$y, 
  postsyn_balancer$z, 
  size = 6, alpha = 1, col = "magenta2", 
  add = T
)

nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
#rgl.snapshot("manuscript/pictures/balancer_sagittal_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Tentacular plane
plot3d(SSN,
       soma = FALSE, lwd = 2, add = TRUE, alpha = 0.5, 
       col = Okabe_Ito[c(1, 5, 8)],
       WithConnectors = F, WithNodes = F)

plot3d(
  balancer, soma = TRUE, color = Okabe_Ito[8], 
  alpha = 0.05, lwd = 2
)

# plot postsynapses
plot3d(
  postsyn_balancer$x, 
  postsyn_balancer$y, 
  postsyn_balancer$z, 
  size = 6, alpha = 1, col = "magenta2", 
  add = T
)

nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))
#rgl.snapshot("manuscript/pictures/balancer_tentacular_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)


#make a snapshot to the working directory
rgl.snapshot("manuscript/pictures/SSN_prepost_synapse_balancer.png")

close3d()

# plot bridge ------------------------------------------------------------------

plot_background()
close3d()
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  #defines the two scenes
#define the size of the rgl window, the view and zoom
#par3d(windowRect = c(0, 0, 1200, 350))
par3d(windowRect = c(0, 0, 2400, 700))

#plot aboral view
plot3d(SSN,
       soma = FALSE, lwd = 2, add = TRUE, alpha = 0.5, 
       col = Okabe_Ito[c(1, 5, 8)],
       WithConnectors = F, WithNodes = F)

plot3d(
  bridge, soma = TRUE, color = Okabe_Ito[8], 
  alpha = 0.05, lwd = 2
)

# plot postsynapses
plot3d(
  postsyn_bridge$x, 
  postsyn_bridge$y, 
  postsyn_bridge$z, 
  size = 6, alpha = 1, col = "magenta2", 
  add = T
)

# plot presynapses
plot3d(
  presyn_bridge$x, 
  presyn_bridge$y, 
  presyn_bridge$z, 
  size = 6, alpha = 1, col = "cyan2", 
  add = T
)

#texts3d(58000,51000,5000, "INRGW", cex = 3, col = "#56B4E9")

#aboral view
nview3d("anterior", extramat = rotationMatrix(1.05, 250, -200, 1000))
#rgl.snapshot("manuscript/pictures/bridge_aboral_view.png")
par3d(zoom=0.67)

#move to next panel in rgl window
next3d(clear=F)




#plot lateral view of Sagittal plane
plot3d(SSN,
       soma = FALSE, lwd = 2, add = TRUE, alpha = 0.5, 
       col = Okabe_Ito[c(1, 5, 8)],
       WithConnectors = F, WithNodes = F)

plot3d(
  bridge, soma = TRUE, color = Okabe_Ito[8], 
  alpha = 0.05, lwd = 2
)

# plot postsynapses
plot3d(
  postsyn_bridge$x, 
  postsyn_bridge$y, 
  postsyn_bridge$z, 
  size = 6, alpha = 1, col = "magenta2", 
  add = T
)

# plot presynapses
plot3d(
  presyn_bridge$x, 
  presyn_bridge$y, 
  presyn_bridge$z, 
  size = 6, alpha = 1, col = "cyan2", 
  add = T
)

nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
#rgl.snapshot("manuscript/pictures/bridge_sagittal_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Tentacular plane
plot3d(SSN,
       soma = FALSE, lwd = 2, add = TRUE, alpha = 0.5, 
       col = Okabe_Ito[c(1, 5, 8)],
       WithConnectors = F, WithNodes = F)

plot3d(
  bridge, soma = TRUE, color = Okabe_Ito[8], 
  alpha = 0.05, lwd = 2
)

# plot postsynapses
plot3d(
  postsyn_bridge$x, 
  postsyn_bridge$y, 
  postsyn_bridge$z, 
  size = 6, alpha = 1, col = "magenta2", 
  add = T
)

# plot presynapses
plot3d(
  presyn_bridge$x, 
  presyn_bridge$y, 
  presyn_bridge$z, 
  size = 6, alpha = 1, col = "cyan2", 
  add = T
)

nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))
#rgl.snapshot("manuscript/pictures/bridge_tentacular_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)


#make a snapshot to the working directory
rgl.snapshot("manuscript/pictures/SSN_prepost_synapse_bridge.png")

close3d()

# plot monociliated syn --------------------------------------------------------

plot_background()
close3d()
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  #defines the two scenes
#define the size of the rgl window, the view and zoom
#par3d(windowRect = c(0, 0, 1200, 350))
par3d(windowRect = c(0, 0, 2400, 700))

#plot aboral view
plot3d(SSN,
       soma = FALSE, lwd = 2, add = TRUE, alpha = 0.5, 
       col = Okabe_Ito[c(1, 5, 8)],
       WithConnectors = F, WithNodes = F)

plot3d(
  monociliated, soma = TRUE, color = Okabe_Ito[8], 
  alpha = 0.05, lwd = 2
)

# plot postsynapses
plot3d(
  postsyn_monoC$x, 
  postsyn_monoC$y, 
  postsyn_monoC$z, 
  size = 6, alpha = 1, col = "magenta2", 
  add = T
)

#texts3d(58000,51000,5000, "INRGW", cex = 3, col = "#56B4E9")

#aboral view
nview3d("anterior", extramat = rotationMatrix(1.05, 250, -200, 1000))
#rgl.snapshot("manuscript/pictures/monoC_aboral_view.png")
par3d(zoom=0.67)

#move to next panel in rgl window
next3d(clear=F)




#plot lateral view of Sagittal plane
plot3d(SSN,
       soma = FALSE, lwd = 2, add = TRUE, alpha = 0.5, 
       col = Okabe_Ito[c(1, 5, 8)],
       WithConnectors = F, WithNodes = F)

plot3d(
  monociliated, soma = TRUE, color = Okabe_Ito[8], 
  alpha = 0.05, lwd = 2
)

# plot postsynapses
plot3d(
  postsyn_monoC$x, 
  postsyn_monoC$y, 
  postsyn_monoC$z, 
  size = 6, alpha = 1, col = "magenta2", 
  add = T
)

nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
#rgl.snapshot("manuscript/pictures/monoC_sagittal_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Tentacular plane
plot3d(SSN,
       soma = FALSE, lwd = 2, add = TRUE, alpha = 0.5, 
       col = Okabe_Ito[c(1, 5, 8)],
       WithConnectors = F, WithNodes = F)

plot3d(
  monociliated, soma = TRUE, color = Okabe_Ito[8], 
  alpha = 0.05, lwd = 2
)

# plot postsynapses
plot3d(
  postsyn_monoC$x, 
  postsyn_monoC$y, 
  postsyn_monoC$z, 
  size = 6, alpha = 1, col = "magenta2", 
  add = T
)

nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))
#rgl.snapshot("manuscript/pictures/monoC_tentacular_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)


#make a snapshot to the working directory
rgl.snapshot("manuscript/pictures/SSN_prepost_synapse_monoC.png")

close3d()

# plot biciliated syn ----------------------------------------------------------

plot_background()
close3d()
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  #defines the two scenes
#define the size of the rgl window, the view and zoom
#par3d(windowRect = c(0, 0, 1200, 350))
par3d(windowRect = c(0, 0, 2400, 700))

#plot aboral view
plot3d(SSN,
       soma = FALSE, lwd = 2, add = TRUE, alpha = 0.5, 
       col = Okabe_Ito[c(1, 5, 8)],
       WithConnectors = F, WithNodes = F)

plot3d(
  biciliated, soma = TRUE, color = Okabe_Ito[8], 
  alpha = 0.05, lwd = 2
)

# plot postsynapses
plot3d(
  postsyn_biC$x, 
  postsyn_biC$y, 
  postsyn_biC$z, 
  size = 6, alpha = 1, col = "magenta2", 
  add = T
)

#texts3d(58000,51000,5000, "INRGW", cex = 3, col = "#56B4E9")

#aboral view
nview3d("anterior", extramat = rotationMatrix(1.05, 250, -200, 1000))
#rgl.snapshot("manuscript/pictures/biC_aboral_view.png")
par3d(zoom=0.67)

#move to next panel in rgl window
next3d(clear=F)




#plot lateral view of Sagittal plane
plot3d(SSN,
       soma = FALSE, lwd = 2, add = TRUE, alpha = 0.5, 
       col = Okabe_Ito[c(1, 5, 8)],
       WithConnectors = F, WithNodes = F)

plot3d(
  biciliated, soma = TRUE, color = Okabe_Ito[8], 
  alpha = 0.05, lwd = 2
)

# plot postsynapses
plot3d(
  postsyn_biC$x, 
  postsyn_biC$y, 
  postsyn_biC$z, 
  size = 6, alpha = 1, col = "magenta2", 
  add = T
)

nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
#rgl.snapshot("manuscript/pictures/biC_sagittal_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Tentacular plane
plot3d(SSN,
       soma = FALSE, lwd = 2, add = TRUE, alpha = 0.5, 
       col = Okabe_Ito[c(1, 5, 8)],
       WithConnectors = F, WithNodes = F)

plot3d(
  biciliated, soma = TRUE, color = Okabe_Ito[8], 
  alpha = 0.05, lwd = 2
)

# plot postsynapses
plot3d(
  postsyn_biC$x, 
  postsyn_biC$y, 
  postsyn_biC$z, 
  size = 6, alpha = 1, col = "magenta2", 
  add = T
)

nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))
#rgl.snapshot("manuscript/pictures/biC_tentacular_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)


#make a snapshot to the working directory
rgl.snapshot("manuscript/pictures/SSN_prepost_synapse_biC.png")

close3d()


# plot nonciliated syn ---------------------------------------------------------

plot_background()
close3d()
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  #defines the two scenes
#define the size of the rgl window, the view and zoom
#par3d(windowRect = c(0, 0, 1200, 350))
par3d(windowRect = c(0, 0, 2400, 700))

#plot aboral view
plot3d(SSN,
       soma = FALSE, lwd = 2, add = TRUE, alpha = 0.5, 
       col = Okabe_Ito[c(1, 5, 8)],
       WithConnectors = F, WithNodes = F)

plot3d(
  nonciliated, soma = TRUE, color = Okabe_Ito[8], 
  alpha = 0.05, lwd = 2
)

# plot postsynapses
plot3d(
  postsyn_nonC$x, 
  postsyn_nonC$y, 
  postsyn_nonC$z, 
  size = 6, alpha = 1, col = "magenta2", 
  add = T
)

#texts3d(58000,51000,5000, "INRGW", cex = 3, col = "#56B4E9")

#aboral view
nview3d("anterior", extramat = rotationMatrix(1.05, 250, -200, 1000))
#rgl.snapshot("manuscript/pictures/nonC_aboral_view.png")
par3d(zoom=0.67)

#move to next panel in rgl window
next3d(clear=F)




#plot lateral view of Sagittal plane
plot3d(SSN,
       soma = FALSE, lwd = 2, add = TRUE, alpha = 0.5, 
       col = Okabe_Ito[c(1, 5, 8)],
       WithConnectors = F, WithNodes = F)

plot3d(
  nonciliated, soma = TRUE, color = Okabe_Ito[8], 
  alpha = 0.05, lwd = 2
)

# plot postsynapses
plot3d(
  postsyn_nonC$x, 
  postsyn_nonC$y, 
  postsyn_nonC$z, 
  size = 6, alpha = 1, col = "magenta2", 
  add = T
)

nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
#rgl.snapshot("manuscript/pictures/nonC_sagittal_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Tentacular plane
plot3d(SSN,
       soma = FALSE, lwd = 2, add = TRUE, alpha = 0.5, 
       col = Okabe_Ito[c(1, 5, 8)],
       WithConnectors = F, WithNodes = F)

plot3d(
  nonciliated, soma = TRUE, color = Okabe_Ito[8], 
  alpha = 0.05, lwd = 2
)

# plot postsynapses
plot3d(
  postsyn_nonC$x, 
  postsyn_nonC$y, 
  postsyn_nonC$z, 
  size = 6, alpha = 1, col = "magenta2", 
  add = T
)

nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))
#rgl.snapshot("manuscript/pictures/nonC_tentacular_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)


#make a snapshot to the working directory
rgl.snapshot("manuscript/pictures/SSN_prepost_synapse_nonC.png")

close3d()

# use visNetwork to plot the network -------------------------------------------

write_rds(AO_graph, "manuscript/source_data/celltype_synapse_matrix.rds")
AO_graph <- read_rds("manuscript/source_data/celltype_synapse_matrix.rds")

## convert to VisNetwork-list
AO_graph.visn <- toVisNetworkData(AO_graph)

#filter low-weight edges
AO_graph.visn$edges$weight
#AO_graph.visn$edges$weight[AO_graph.visn$edges$weight < 3]  <- 0
AO_graph.visn$edges$weight <- sqrt(AO_graph.visn$edges$weight)
AO_graph.visn$edges$weight

## copy column "weight" to new column "value" in list "edges"
AO_graph.visn$edges$value <- AO_graph.visn$edges$weight
AO_graph.visn$nodes$value <- degree

#define node color
AO_graph.visn$nodes$color <- Okabe_Ito[1:17]

#hierarchical layout - define level of nodes
AO_graph.visn$nodes$level <- c(2, 4, 1, 1, 3, 3, 3, 3, 3, 4)
#balancer, LB, syn_neuron, bridge, monociliated, biciliated, nonciliated, dense_vesicle

#hierarchical layout
visNet <- visNetwork(AO_graph.visn$nodes, AO_graph.visn$edges) %>%
  visIgraphLayout(
    layout = "layout_nicely", physics = FALSE
    ) %>%
  visHierarchicalLayout(
    levelSeparation=250, 
    nodeSpacing=200,
    direction='LR',
    sortMethod='hubsize',
    shakeTowards='roots'
    ) %>%
  visEdges(
    smooth = list(type = 'curvedCW', roundness=0.2),
    scaling=list(min=2, max=12),
    color = list(inherit=TRUE, opacity=0.7),
    arrows = list(
      to = list(enabled = TRUE, 
      scaleFactor = 1, type = 'arrow'))
    ) %>%
  visNodes(
    borderWidth=0.3, 
    color = list(background=AO_graph.visn$nodes$color, border='black'),
    opacity=0.9,
    shape='dot', 
    font=list(color='black', size=44),
    scaling = list(label=list(enabled=TRUE, min=22, max=80)),
    level= AO_graph.visn$nodes$level
    ) %>%
  visOptions(highlightNearest = TRUE, width = 1800, height = 1800)
visNet


saveNetwork(visNet, "manuscript/pictures/visNetwork_celltype_connectome.html")
webshot2::webshot(url = "manuscript/pictures/visNetwork_celltype_connectome.html",
                  file = "manuscript/pictures/visNetwork_celltype_connectome.png",
                  vwidth = 1800, vheight = 1800, #define the size of the browser window
                  cliprect = c(390, 250, 880, 1000), zoom = 1, delay = 1)


