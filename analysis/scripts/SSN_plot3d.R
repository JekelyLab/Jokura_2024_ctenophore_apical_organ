# Code to generate SSN 3d plot

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


balancer_Q1 <- nlapply(read.neurons.catmaid(get_skids_with_annot(pid = 35, c("celltype:balancer", "Q1")),
                                            pid = 35),
                       function(x)
                         smooth_neuron(x, sigma = 1000))
balancer_Q2 <- nlapply(read.neurons.catmaid(get_skids_with_annot(pid = 35, c("celltype:balancer", "Q2")),
                                            pid = 35),
                       function(x)
                         smooth_neuron(x, sigma = 1000))
balancer_Q3 <- nlapply(read.neurons.catmaid(get_skids_with_annot(pid = 35, c("celltype:balancer", "Q3")),
                                            pid = 35),
                       function(x)
                         smooth_neuron(x, sigma = 1000))
balancer_Q4 <- nlapply(read.neurons.catmaid(get_skids_with_annot(pid = 35, c("celltype:balancer", "Q4")),
                                            pid = 35),
                       function(x)
                         smooth_neuron(x, sigma = 1000))


all_celltypes <- list(balancer, bridge, bristle, dome, Cgroove_tag, Cgroove_sag,
                      intramulticilia, lamellate, lithocyte, neuron, 
                      plumose, dense_vesicle, monocilia, bicilia, non_cilia)


# plot SSN neuron -----------------------------------------------

plot_background()
close3d()
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  #defines the two scenes
#define the size of the rgl window, the view and zoom
#par3d(windowRect = c(0, 0, 1200, 350))
par3d(windowRect = c(0, 0, 2400, 700))

#plot aboral view
plot3d(neuron,
       soma = TRUE, lwd = 1, add = T, alpha = 0.5, color = Okabe_Ito[c(1,5,8)],
       WithConnectors = F, WithNodes = F)

#for (object in objects) {
#  plot3d(object, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.05, col = Okabe_Ito[8])
#}


plot3d(balancer, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(bridge, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(bristle, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(dome, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(Cgroove_tag, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(Cgroove_sag, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(intramulticilia, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(lamellate, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(lithocyte, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(plumose, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(dense_vesicle, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(monocilia, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(bicilia, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(non_cilia, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])

#texts3d(58000,51000,5000, "INRGW", cex = 3, col = "#56B4E9")

#aboral view
nview3d("anterior", extramat = rotationMatrix(1.05, 250, -200, 1000))
#rgl.snapshot("manuscript/pictures/balancer_aboral_view.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Sagittal plane
plot3d(neuron,
       soma = TRUE, lwd = 1, add = T, alpha = 0.5, color = Okabe_Ito[c(1,5,8)],
       WithConnectors = F, WithNodes = F)

#for (object in objects) {
#  plot3d(object, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.05, col = Okabe_Ito[8])
#}

plot3d(balancer, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(bridge, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(bristle, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(dome, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(Cgroove_tag, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(Cgroove_sag, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(intramulticilia, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(lamellate, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(lithocyte, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(plumose, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(dense_vesicle, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(monocilia, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(bicilia, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(non_cilia, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])




nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
#rgl.snapshot("manuscript/pictures/balancer_sagittal_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Tentacular plane
plot3d(neuron,
       soma = TRUE, lwd = 1, add = T, alpha = 0.5, color = Okabe_Ito[c(1,5,8)],
       WithConnectors = F, WithNodes = F)

#for (object in objects) {
#  plot3d(object, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.05, col = Okabe_Ito[8])
#}

plot3d(balancer, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(bridge, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(bristle, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(dome, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(Cgroove_tag, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(Cgroove_sag, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(intramulticilia, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(lamellate, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(lithocyte, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(plumose, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(dense_vesicle, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(monocilia, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(bicilia, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(non_cilia, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])



nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))
#rgl.snapshot("manuscript/pictures/balancer_tentacular_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)


#make a snapshot to the working directory
rgl.snapshot("manuscript/pictures/SSN.png")


close3d()

