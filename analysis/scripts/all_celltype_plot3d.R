# Code to generate all celltype 3d plot

# source packages and functions -----------------
source("analysis/scripts/packages_and_functions.R")

# load cell type ----------------------------------

read_smooth_neuron <- function(annotation){
  nlapply(read.neurons.catmaid(annotation, pid = 35),
          function(x)
            smooth_neuron(x, sigma = 1000))
}

balancer <- read_smooth_neuron("celltype:balancer")
bridge <- read_smooth_neuron("celltype:bridge")
bristle <- read_smooth_neuron("celltype:bristle")
dome <- read_smooth_neuron("celltype:dome")
Cgroove_tag <- read_smooth_neuron("celltype:Cgroove-tag")
Cgroove_sag <- read_smooth_neuron("celltype:Cgroove-sag")
intramulticilia <- read_smooth_neuron("celltype:intra-multi-ciliated")
lamellate <- read_smooth_neuron("celltype:lamellate")
lithocyte <- read_smooth_neuron("celltype:lithocyte")
neuron <- read_smooth_neuron("celltype:SSN")
plumose <- read_smooth_neuron("celltype:plumose")
dense_vesicle <- read_smooth_neuron("celltype:dense_vesicle")
monocilia <- read_smooth_neuron("celltype:monociliated")
bicilia <- read_smooth_neuron("celltype:biciliated")
non_cilia <- read_smooth_neuron("celltype:nonciliated")


Q1 <- nlapply(read.neurons.catmaid("Q1", pid = 35),
              function(x)
                smooth_neuron(x, sigma = 1000))
Q2 <- nlapply(read.neurons.catmaid("Q2", pid = 35),
              function(x)
                smooth_neuron(x, sigma = 1000))
Q3 <- nlapply(read.neurons.catmaid("Q3", pid = 35),
              function(x)
                smooth_neuron(x, sigma = 1000))
Q4 <- nlapply(read.neurons.catmaid("Q4", pid = 35),
              function(x)
                smooth_neuron(x, sigma = 1000))


balancer_Q1 <- nlapply(read.neurons.catmaid(skids_by_2annotations("celltype:balancer", "Q1"),
                                            pid = 35),
                       function(x)
                         smooth_neuron(x, sigma = 1000))
balancer_Q2 <- nlapply(read.neurons.catmaid(skids_by_2annotations("celltype:balancer", "Q2"),
                                            pid = 35),
                       function(x)
                         smooth_neuron(x, sigma = 1000))
balancer_Q3 <- nlapply(read.neurons.catmaid(skids_by_2annotations("celltype:balancer", "Q3"),
                                            pid = 35),
                       function(x)
                         smooth_neuron(x, sigma = 1000))
balancer_Q4 <- nlapply(read.neurons.catmaid(skids_by_2annotations("celltype:balancer", "Q4"),
                                            pid = 35),
                       function(x)
                         smooth_neuron(x, sigma = 1000))


all_celltypes <- list(balancer, bridge, bristle, dome, Cgroove_tag, Cgroove_sag,
                      intramulticilia, lamellate, lithocyte, neuron, 
                      plumose, dense_vesicle, monocilia, bicilia, non_cilia)


# plot all cells ----------------------

colour_palettes <- c(Okabe_Ito, bluepurple, oranges)

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

