

# source packages and functions -----------------
source("analysis/scripts/packages_and_functions.R")

# load cell type ----------------------------------

#thirteen celltypes
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

#four ciliated celltypes
monociliated <- read_smooth_neuron("celltype:monociliated")
biciliated <- read_smooth_neuron("celltype:biciliated")
multiciliated <- read_smooth_neuron("celltype:multiciliated")
nonciliated <- read_smooth_neuron("celltype:nonciliated")


# plot balancer -----------------------------------------------

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
rgl.snapshot("manuscript/pictures/plot_balancer.png")


close3d()


# plot bridge -----------------------------------------------



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
#rgl.snapshot("manuscript/pictures/bridge_aboral_view.png")
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
#rgl.snapshot("manuscript/pictures/bridge_sagittal_plane.png")
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
#rgl.snapshot("manuscript/pictures/bridge_tentacular_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)


#make a snapshot to the working directory
rgl.snapshot("manuscript/pictures/plot_bridge.png")


close3d()

# plot bristle -----------------------------------------------

close3d()
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  #defines the two scenes
#define the size of the rgl window, the view and zoom
#par3d(windowRect = c(0, 0, 1200, 350))
par3d(windowRect = c(0, 0, 2400, 700))

#plot aboral view
plot3d(bristle,
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
#rgl.snapshot("manuscript/pictures/bristle_aboral_view.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Sagittal plane
plot3d(bristle,
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
#rgl.snapshot("manuscript/pictures/bristle_sagittal_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Tentacular plane
plot3d(bristle,
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
#rgl.snapshot("manuscript/pictures/bristle_tentacular_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)


#make a snapshot to the working directory
rgl.snapshot("manuscript/pictures/plot_bristle.png")


close3d()

# plot ciliated groove -----------------------------------------------

close3d()
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  #defines the two scenes
#define the size of the rgl window, the view and zoom
#par3d(windowRect = c(0, 0, 1200, 350))
par3d(windowRect = c(0, 0, 2400, 700))

#plot aboral view
plot3d(Cgroove_sag,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)
plot3d(Cgroove_tag,
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
#rgl.snapshot("manuscript/pictures/Cgroove_aboral_view.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Sagittal plane
plot3d(Cgroove_sag,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)
plot3d(Cgroove_tag,
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
#rgl.snapshot("manuscript/pictures/Cgroove_sagittal_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Tentacular plane
plot3d(Cgroove_sag,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)
plot3d(Cgroove_tag,
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
#rgl.snapshot("manuscript/pictures/Cgroove_tentacular_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)


#make a snapshot to the working directory
rgl.snapshot("manuscript/pictures/plot_Cgroove.png")


close3d()

# plot dense vesicle -----------------------------------------------

close3d()
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  #defines the two scenes
#define the size of the rgl window, the view and zoom
#par3d(windowRect = c(0, 0, 1200, 350))
par3d(windowRect = c(0, 0, 2400, 700))

#plot aboral view
plot3d(dense_vesicle,
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
#rgl.snapshot("manuscript/pictures/dense_vesicle_aboral_view.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Sagittal plane
plot3d(dense_vesicle,
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
#rgl.snapshot("manuscript/pictures/dense_vesicle_sagittal_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Tentacular plane
plot3d(dense_vesicle,
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
#rgl.snapshot("manuscript/pictures/dense_vesicle_tentacular_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)


#make a snapshot to the working directory
rgl.snapshot("manuscript/pictures/plot_dense_vesicle.png")


close3d()


# plot dome -----------------------------------------------

close3d()
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  #defines the two scenes
#define the size of the rgl window, the view and zoom
#par3d(windowRect = c(0, 0, 1200, 350))
par3d(windowRect = c(0, 0, 2400, 700))

#plot aboral view
plot3d(dome,
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
#rgl.snapshot("manuscript/pictures/dome_aboral_view.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Sagittal plane
plot3d(dome,
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
#rgl.snapshot("manuscript/pictures/dome_sagittal_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Tentacular plane
plot3d(dome,
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
#rgl.snapshot("manuscript/pictures/dome_tentacular_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)


#make a snapshot to the working directory
rgl.snapshot("manuscript/pictures/plot_dome.png")


close3d()


# plot epithelial floor -----------------------------------------------

close3d()
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  #defines the two scenes
#define the size of the rgl window, the view and zoom
#par3d(windowRect = c(0, 0, 1200, 350))
par3d(windowRect = c(0, 0, 2400, 700))

#plot aboral view
plot3d(epithelial_floor,
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
#rgl.snapshot("manuscript/pictures/epithelial_floor_aboral_view.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Sagittal plane
plot3d(epithelial_floor,
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
#rgl.snapshot("manuscript/pictures/epithelial_floor_sagittal_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Tentacular plane
plot3d(epithelial_floor,
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
#rgl.snapshot("manuscript/pictures/epithelial_floor_tentacular_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)


#make a snapshot to the working directory
rgl.snapshot("manuscript/pictures/plot_epithelial_floor.png")


close3d()


# plot intra_multi_ciliated -----------------------------------------------

close3d()
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  #defines the two scenes
#define the size of the rgl window, the view and zoom
#par3d(windowRect = c(0, 0, 1200, 350))
par3d(windowRect = c(0, 0, 2400, 700))

#plot aboral view
plot3d(intra_multi_ciliated,
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
#rgl.snapshot("manuscript/pictures/intra_multi_ciliated_aboral_view.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Sagittal plane
plot3d(intra_multi_ciliated,
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
#rgl.snapshot("manuscript/pictures/intra_multi_ciliated_sagittal_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Tentacular plane
plot3d(intra_multi_ciliated,
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
#rgl.snapshot("manuscript/pictures/intra_multi_ciliated_tentacular_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)


#make a snapshot to the working directory
rgl.snapshot("manuscript/pictures/plot_intra_multi_ciliated.png")


close3d()


# plot lamellate -----------------------------------------------

close3d()
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  #defines the two scenes
#define the size of the rgl window, the view and zoom
#par3d(windowRect = c(0, 0, 1200, 350))
par3d(windowRect = c(0, 0, 2400, 700))

#plot aboral view
plot3d(lamellate,
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
#rgl.snapshot("manuscript/pictures/lamellate_aboral_view.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Sagittal plane
plot3d(lamellate,
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
#rgl.snapshot("manuscript/pictures/lamellate_sagittal_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Tentacular plane
plot3d(lamellate,
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
#rgl.snapshot("manuscript/pictures/lamellate_tentacular_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)


#make a snapshot to the working directory
rgl.snapshot("manuscript/pictures/plot_lamellate.png")


close3d()



# plot lithocyte -----------------------------------------------

close3d()
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  #defines the two scenes
#define the size of the rgl window, the view and zoom
#par3d(windowRect = c(0, 0, 1200, 350))
par3d(windowRect = c(0, 0, 2400, 700))

#plot aboral view
plot3d(lithocyte,
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
#rgl.snapshot("manuscript/pictures/lithocyte_aboral_view.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Sagittal plane
plot3d(lithocyte,
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
#rgl.snapshot("manuscript/pictures/lithocyte_sagittal_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Tentacular plane
plot3d(lithocyte,
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
#rgl.snapshot("manuscript/pictures/lithocyte_tentacular_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)


#make a snapshot to the working directory
rgl.snapshot("manuscript/pictures/plot_lithocyte.png")


close3d()


# plot plumose -----------------------------------------------

close3d()
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  #defines the two scenes
#define the size of the rgl window, the view and zoom
#par3d(windowRect = c(0, 0, 1200, 350))
par3d(windowRect = c(0, 0, 2400, 700))

#plot aboral view
plot3d(plumose,
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
#rgl.snapshot("manuscript/pictures/plumose_aboral_view.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Sagittal plane
plot3d(plumose,
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
#rgl.snapshot("manuscript/pictures/plumose_sagittal_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Tentacular plane
plot3d(plumose,
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
#rgl.snapshot("manuscript/pictures/plumose_tentacular_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)


#make a snapshot to the working directory
rgl.snapshot("manuscript/pictures/plot_plumose.png")


close3d()


# plot SSN -----------------------------------------------

close3d()
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  #defines the two scenes
#define the size of the rgl window, the view and zoom
#par3d(windowRect = c(0, 0, 1200, 350))
par3d(windowRect = c(0, 0, 2400, 700))

#plot aboral view
plot3d(SSN,
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
#rgl.snapshot("manuscript/pictures/SSN_aboral_view.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Sagittal plane
plot3d(SSN,
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
#rgl.snapshot("manuscript/pictures/SSN_sagittal_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Tentacular plane
plot3d(SSN,
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
#rgl.snapshot("manuscript/pictures/SSN_tentacular_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)


#make a snapshot to the working directory
rgl.snapshot("manuscript/pictures/plot_SSN.png")


close3d()


# plot monocilia -----------------------------------------------

close3d()
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  #defines the two scenes
#define the size of the rgl window, the view and zoom
#par3d(windowRect = c(0, 0, 1200, 350))
par3d(windowRect = c(0, 0, 2400, 700))

#plot aboral view
plot3d(monociliated,
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
#rgl.snapshot("manuscript/pictures/monociliated_aboral_view.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Sagittal plane
plot3d(monociliated,
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
#rgl.snapshot("manuscript/pictures/monociliated_sagittal_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Tentacular plane
plot3d(monociliated,
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
#rgl.snapshot("manuscript/pictures/monociliated_tentacular_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)


#make a snapshot to the working directory
rgl.snapshot("manuscript/pictures/plot_monociliated.png")


close3d()


# plot bicilia -----------------------------------------------

close3d()
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  #defines the two scenes
#define the size of the rgl window, the view and zoom
#par3d(windowRect = c(0, 0, 1200, 350))
par3d(windowRect = c(0, 0, 2400, 700))

#plot aboral view
plot3d(biciliated,
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
#rgl.snapshot("manuscript/pictures/biciliated_aboral_view.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Sagittal plane
plot3d(biciliated,
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
#rgl.snapshot("manuscript/pictures/biciliated_sagittal_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Tentacular plane
plot3d(biciliated,
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
#rgl.snapshot("manuscript/pictures/biciliated_tentacular_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)


#make a snapshot to the working directory
rgl.snapshot("manuscript/pictures/plot_biciliated.png")


close3d()


# plot multiciliated -----------------------------------------------

close3d()
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  #defines the two scenes
#define the size of the rgl window, the view and zoom
#par3d(windowRect = c(0, 0, 1200, 350))
par3d(windowRect = c(0, 0, 2400, 700))

#plot aboral view
plot3d(multiciliated,
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
#rgl.snapshot("manuscript/pictures/multiciliated_aboral_view.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Sagittal plane
plot3d(multiciliated,
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
#rgl.snapshot("manuscript/pictures/multiciliated_sagittal_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Tentacular plane
plot3d(multiciliated,
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
#rgl.snapshot("manuscript/pictures/multiciliated_tentacular_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)


#make a snapshot to the working directory
rgl.snapshot("manuscript/pictures/plot_multiciliated.png")


close3d()

# plot non_cilia -----------------------------------------------

close3d()
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  #defines the two scenes
#define the size of the rgl window, the view and zoom
#par3d(windowRect = c(0, 0, 1200, 350))
par3d(windowRect = c(0, 0, 2400, 700))

#plot aboral view
plot3d(nonciliated,
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
#rgl.snapshot("manuscript/pictures/nonciliated_aboral_view.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Sagittal plane
plot3d(nonciliated,
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
#rgl.snapshot("manuscript/pictures/nonciliated_sagittal_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Tentacular plane
plot3d(nonciliated,
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
#rgl.snapshot("manuscript/pictures/nonciliated_tentacular_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)


#make a snapshot to the working directory
rgl.snapshot("manuscript/pictures/plot_nonciliated.png")


close3d()


# assemble figure -------------------------------------------------------------

#read pics
panel_balancer <- ggdraw() + draw_image(readPNG("manuscript/pictures/plot_balancer.png")) +
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

panel_bridge <- ggdraw() + draw_image(readPNG("manuscript/pictures/plot_bridge.png")) +
  draw_label("bridge cells", x = 0.5, y = 0.95, size = 8, fontface="bold", hjust = 0.5)
panel_bristle <- ggdraw() + draw_image(readPNG("manuscript/pictures/plot_bristle.png")) +
  draw_label("bristle cells", x = 0.5, y = 0.95, size = 8, fontface="bold", hjust = 0.5)
panel_Cgroove <- ggdraw() + draw_image(readPNG("manuscript/pictures/plot_Cgroove.png")) +
  draw_label("ciliated groove cells", x = 0.5, y = 0.95, size = 8, fontface="bold", hjust = 0.5)
panel_dense_vesicle <- ggdraw() + draw_image(readPNG("manuscript/pictures/plot_dense_vesicle.png")) +
  draw_label("dense vesicle cells", x = 0.5, y = 0.95, size = 8, fontface="bold", hjust = 0.5)
panel_dome <- ggdraw() + draw_image(readPNG("manuscript/pictures/plot_dome.png")) +
  draw_label("dome cells", x = 0.5, y = 0.95, size = 8, fontface="bold", hjust = 0.5)
panel_intra_multi_ciliated <- ggdraw() + draw_image(readPNG("manuscript/pictures/plot_intra_multi_ciliated.png")) +
  draw_label("intra-multiciliated cells", x = 0.5, y = 0.95, size = 8, fontface="bold", hjust = 0.5)
panel_lamellate <- ggdraw() + draw_image(readPNG("manuscript/pictures/plot_lamellate.png")) +
  draw_label("lamellate bodies", x = 0.5, y = 0.95, size = 8, fontface="bold", hjust = 0.5)
panel_lithocyte <- ggdraw() + draw_image(readPNG("manuscript/pictures/plot_lithocyte.png")) +
  draw_label("lithocytes", x = 0.5, y = 0.95, size = 8, fontface="bold", hjust = 0.5)
panel_plumose <- ggdraw() + draw_image(readPNG("manuscript/pictures/plot_plumose.png")) +
  draw_label("plumose cells", x = 0.5, y = 0.95, size = 8, fontface="bold", hjust = 0.5)
panel_SSN <- ggdraw() + draw_image(readPNG("manuscript/pictures/plot_SSN.png")) +
  draw_label("AO neurons", x = 0.5, y = 0.95, size = 8, fontface="bold", hjust = 0.5)
panel_monociliated <- ggdraw() + draw_image(readPNG("manuscript/pictures/plot_monociliated.png")) +
  draw_label("monociliated cells", x = 0.5, y = 0.95, size = 8, fontface="bold", hjust = 0.5)
panel_biciliated <- ggdraw() + draw_image(readPNG("manuscript/pictures/plot_biciliated.png")) +
  draw_label("biciliated cells", x = 0.5, y = 0.95, size = 8, fontface="bold", hjust = 0.5)
panel_multiciliated <- ggdraw() + draw_image(readPNG("manuscript/pictures/plot_multiciliated.png")) +
  draw_label("multiciliated cells", x = 0.5, y = 0.95, size = 8, fontface="bold", hjust = 0.5)
panel_nonciliated <- ggdraw() + draw_image(readPNG("manuscript/pictures/plot_nonciliated.png")) +
  draw_label("nonciliated cells", x = 0.5, y = 0.95, size = 8, fontface="bold", hjust = 0.5)

layout <- "
ABC
DEF
GHI
JKL
MNO
"

all_assemble <- panel_balancer + panel_bridge + panel_bristle +
  panel_Cgroove + panel_dense_vesicle + panel_dome + 
  panel_intra_multi_ciliated + panel_lamellate + panel_lithocyte + 
  panel_plumose + panel_SSN + panel_monociliated + 
  panel_biciliated + panel_multiciliated + panel_nonciliated +
  patchwork::plot_layout(design = layout, 
                         heights = c(1,1,1,1,1), 
                         widths = c(1,1,1)) + 
  patchwork::plot_annotation(tag_levels = "A") &  
  ggplot2::theme(plot.tag = element_text(size = 12, 
                                         face='plain', color='black'))


ggsave("manuscript/pictures/all_assemble.png", limitsize = FALSE, 
       units = c("px"), all_assemble, width = 3800, height = 2400, bg='white')  

