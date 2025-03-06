# Code to generate all cell types in 3d plot and 
#save the views from three different directions as png files

# source packages and functions -----------------
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



# plot all cells ----------------------

colour_palettes <- c(bluepurple, oranges, Okabe_Ito)

nopen3d() 
mfrow3d(1, 3)  #defines the two scenes
#define the size of the rgl window, the view and zoom
#par3d(windowRect = c(0, 0, 1200, 350))
par3d(windowRect = c(0, 0, 2400, 700))


#plot aboral view
for (i in 1:length(all_celltypes)) {
  print(i)
  plot3d(
    all_celltypes[[i]], soma = TRUE, lwd = 1, add = TRUE, 
    alpha = 0.6, col = colour_palettes[i]
  )
}
nview3d("anterior", extramat = rotationMatrix(1.05, 250, -200, 1000))
par3d(zoom=0.61)
next3d(clear=F)

#plot lateral view of Sagittal plane
for (i in 1:length(all_celltypes)) {
  print(i)
  plot3d(
    all_celltypes[[i]], soma = TRUE, lwd = 1, add = TRUE, 
    alpha = 0.6, col = colour_palettes[i]
  )
}
nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
par3d(zoom=0.61)

next3d(clear=F)

# plot screen 3
for (i in 1:length(all_celltypes)) {
  print(i)
  plot3d(
    all_celltypes[[i]], soma = TRUE, lwd = 1, add = TRUE, 
    alpha = 0.6, col = colour_palettes[i]
  )
}
nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))
par3d(zoom=0.61)

#make a snapshot
rgl.snapshot("manuscript/pictures/all_cells_3_views.png")
close3d()

