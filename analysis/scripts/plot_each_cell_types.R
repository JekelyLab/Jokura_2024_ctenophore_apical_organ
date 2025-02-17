

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
par3d(windowRect = c(0, 0, 2400, 700))


#plot aboral view
plot3d(bridge,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)
for (object in objects) {
  plot3d(object, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.05, col = Okabe_Ito[8])
}

nview3d("anterior", extramat = rotationMatrix(1.05, 250, -200, 1000))
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)

#plot lateral view of Sagittal plane
plot3d(bridge,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)
for (object in objects) {
  plot3d(object, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.05, col = Okabe_Ito[8])
}

nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)

#plot lateral view of Tentacular plane
plot3d(bridge,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)
for (object in objects) {
  plot3d(object, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.05, col = Okabe_Ito[8])
}

nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))
par3d(zoom=0.61)

#make a snapshot
rgl.snapshot("manuscript/pictures/bridge.png")

# plot bristle -----------------------------------------------

close3d()
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  #defines the two scenes
#define the size of the rgl window, the view and zoom
par3d(windowRect = c(0, 0, 2400, 700))


#plot aboral view
plot3d(bristle,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)
for (object in objects) {
  plot3d(object, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.05, col = Okabe_Ito[8])
}

nview3d("anterior", extramat = rotationMatrix(1.05, 250, -200, 1000))
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)

#plot lateral view of Sagittal plane
plot3d(bristle,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)
for (object in objects) {
  plot3d(object, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.05, col = Okabe_Ito[8])
}

nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)

#plot lateral view of Tentacular plane
plot3d(bristle,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)
for (object in objects) {
  plot3d(object, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.05, col = Okabe_Ito[8])
}

nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))
par3d(zoom=0.61)


#make a snapshot
rgl.snapshot("manuscript/pictures/bristle.png")


# plot dome -----------------------------------------------

close3d()
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  #defines the two scenes
#define the size of the rgl window, the view and zoom
par3d(windowRect = c(0, 0, 2400, 700))


#plot aboral view
plot3d(dome,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)
for (object in objects) {
  plot3d(object, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.05, col = Okabe_Ito[8])
}

nview3d("anterior", extramat = rotationMatrix(1.05, 250, -200, 1000))
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)

#plot lateral view of Sagittal plane
plot3d(dome,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)
for (object in objects) {
  plot3d(object, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.05, col = Okabe_Ito[8])
}

nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)

#plot lateral view of Tentacular plane
plot3d(dome,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)
for (object in objects) {
  plot3d(object, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.05, col = Okabe_Ito[8])
}

nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))
par3d(zoom=0.61)


#make a snapshot
rgl.snapshot("manuscript/pictures/dome.png")

# plot groove -----------------------------------------------

close3d()
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  #defines the two scenes
#define the size of the rgl window, the view and zoom
par3d(windowRect = c(0, 0, 2400, 700))


#plot aboral view
plot3d(groove,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)
for (object in objects) {
  plot3d(object, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.05, col = Okabe_Ito[8])
}

nview3d("anterior", extramat = rotationMatrix(1.05, 250, -200, 1000))
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)

#plot lateral view of Sagittal plane
plot3d(groove,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)
for (object in objects) {
  plot3d(object, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.05, col = Okabe_Ito[8])
}

nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)

#plot lateral view of Tentacular plane
plot3d(groove,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)
for (object in objects) {
  plot3d(object, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.05, col = Okabe_Ito[8])
}

nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))
par3d(zoom=0.61)


#make a snapshot
rgl.snapshot("manuscript/pictures/groove.png")

# plot intramulticilia -----------------------------------------------


close3d()
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  #defines the two scenes
#define the size of the rgl window, the view and zoom
par3d(windowRect = c(0, 0, 2400, 700))


#plot aboral view
plot3d(intramulticilia,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)
for (object in objects) {
  plot3d(object, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.05, col = Okabe_Ito[8])
}

nview3d("anterior", extramat = rotationMatrix(1.05, 250, -200, 1000))
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)

#plot lateral view of Sagittal plane
plot3d(intramulticilia,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)
for (object in objects) {
  plot3d(object, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.05, col = Okabe_Ito[8])
}

nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)

#plot lateral view of Tentacular plane
plot3d(intramulticilia,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)
for (object in objects) {
  plot3d(object, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.05, col = Okabe_Ito[8])
}

nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))
par3d(zoom=0.61)


#make a snapshot
rgl.snapshot("manuscript/pictures/intramulticilia.png")

# plot lamellate -----------------------------------------------

close3d()
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  #defines the two scenes
#define the size of the rgl window, the view and zoom
par3d(windowRect = c(0, 0, 2400, 700))


#plot aboral view
plot3d(lamellate,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)
for (object in objects) {
  plot3d(object, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.05, col = Okabe_Ito[8])
}

nview3d("anterior", extramat = rotationMatrix(1.05, 250, -200, 1000))
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)

#plot lateral view of Sagittal plane
plot3d(lamellate,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)
for (object in objects) {
  plot3d(object, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.05, col = Okabe_Ito[8])
}

nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)

#plot lateral view of Tentacular plane
plot3d(lamellate,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)
for (object in objects) {
  plot3d(object, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.05, col = Okabe_Ito[8])
}

nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))
par3d(zoom=0.61)


#make a snapshot
rgl.snapshot("manuscript/pictures/lamellate.png")

# plot lithocyte -----------------------------------------------

close3d()
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  #defines the two scenes
#define the size of the rgl window, the view and zoom
par3d(windowRect = c(0, 0, 2400, 700))


#plot aboral view
plot3d(lithocyte,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)
for (object in objects) {
  plot3d(object, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.05, col = Okabe_Ito[8])
}

nview3d("anterior", extramat = rotationMatrix(1.05, 250, -200, 1000))
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)

#plot lateral view of Sagittal plane
plot3d(lithocyte,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)
for (object in objects) {
  plot3d(object, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.05, col = Okabe_Ito[8])
}

nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)

#plot lateral view of Tentacular plane
plot3d(lithocyte,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)
for (object in objects) {
  plot3d(object, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.05, col = Okabe_Ito[8])
}

nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))
par3d(zoom=0.61)


#make a snapshot
rgl.snapshot("manuscript/pictures/lithocyte.png")

# plot neuron -----------------------------------------------


close3d()
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  #defines the two scenes
#define the size of the rgl window, the view and zoom
par3d(windowRect = c(0, 0, 2400, 700))


#plot aboral view
plot3d(neuron,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)
for (object in objects) {
  plot3d(object, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.05, col = Okabe_Ito[8])
}

nview3d("anterior", extramat = rotationMatrix(1.05, 250, -200, 1000))
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)

#plot lateral view of Sagittal plane
plot3d(neuron,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)
for (object in objects) {
  plot3d(object, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.05, col = Okabe_Ito[8])
}

nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)

#plot lateral view of Tentacular plane
plot3d(neuron,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)
for (object in objects) {
  plot3d(object, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.05, col = Okabe_Ito[8])
}

nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))
par3d(zoom=0.61)


#make a snapshot
rgl.snapshot("manuscript/pictures/neuron.png")

# plot plumose -----------------------------------------------

close3d()
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  #defines the two scenes
#define the size of the rgl window, the view and zoom
par3d(windowRect = c(0, 0, 2400, 700))


#plot aboral view
plot3d(plumose,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)
for (object in objects) {
  plot3d(object, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.05, col = Okabe_Ito[8])
}

nview3d("anterior", extramat = rotationMatrix(1.05, 250, -200, 1000))
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)

#plot lateral view of Sagittal plane
plot3d(plumose,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)
for (object in objects) {
  plot3d(object, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.05, col = Okabe_Ito[8])
}

nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)

#plot lateral view of Tentacular plane
plot3d(plumose,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)
for (object in objects) {
  plot3d(object, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.05, col = Okabe_Ito[8])
}

nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))
par3d(zoom=0.61)


#make a snapshot
rgl.snapshot("manuscript/pictures/plumose.png")

# plot dense_vesicle -----------------------------------------------


close3d()
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  #defines the two scenes
#define the size of the rgl window, the view and zoom
par3d(windowRect = c(0, 0, 2400, 700))


#plot aboral view
plot3d(dense_vesicle,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)
for (object in objects) {
  plot3d(object, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.05, col = Okabe_Ito[8])
}

nview3d("anterior", extramat = rotationMatrix(1.05, 250, -200, 1000))
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)

#plot lateral view of Sagittal plane
plot3d(dense_vesicle,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)
for (object in objects) {
  plot3d(object, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.05, col = Okabe_Ito[8])
}

nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)

#plot lateral view of Tentacular plane
plot3d(dense_vesicle,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)
for (object in objects) {
  plot3d(object, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.05, col = Okabe_Ito[8])
}

nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))
par3d(zoom=0.61)


#make a snapshot
rgl.snapshot("manuscript/pictures/dense_vesicle.png")

# plot monocilia -----------------------------------------------

close3d()
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  #defines the two scenes
#define the size of the rgl window, the view and zoom
par3d(windowRect = c(0, 0, 2400, 700))


#plot aboral view
plot3d(monocilia,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[3],
       WithConnectors = F, WithNodes = F)
for (object in objects) {
  plot3d(object, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.05, col = Okabe_Ito[8])
}

nview3d("anterior", extramat = rotationMatrix(1.05, 250, -200, 1000))
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)

#plot lateral view of Sagittal plane
plot3d(monocilia,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[3],
       WithConnectors = F, WithNodes = F)
for (object in objects) {
  plot3d(object, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.05, col = Okabe_Ito[8])
}

nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)

#plot lateral view of Tentacular plane
plot3d(monocilia,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[3],
       WithConnectors = F, WithNodes = F)
for (object in objects) {
  plot3d(object, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.05, col = Okabe_Ito[8])
}

nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))
par3d(zoom=0.61)


#make a snapshot
rgl.snapshot("manuscript/pictures/monocilia.png")

# plot bicilia -----------------------------------------------

close3d()
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  #defines the two scenes
#define the size of the rgl window, the view and zoom
par3d(windowRect = c(0, 0, 2400, 700))


#plot aboral view
plot3d(bicilia,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[7],
       WithConnectors = F, WithNodes = F)
for (object in objects) {
  plot3d(object, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.05, col = Okabe_Ito[8])
}

nview3d("anterior", extramat = rotationMatrix(1.05, 250, -200, 1000))
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)

#plot lateral view of Sagittal plane
plot3d(bicilia,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[7],
       WithConnectors = F, WithNodes = F)
for (object in objects) {
  plot3d(object, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.05, col = Okabe_Ito[8])
}

nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)

#plot lateral view of Tentacular plane
plot3d(bicilia,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[7],
       WithConnectors = F, WithNodes = F)
for (object in objects) {
  plot3d(object, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.05, col = Okabe_Ito[8])
}

nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))
par3d(zoom=0.61)


#make a snapshot
rgl.snapshot("manuscript/pictures/bicilia.png")


# plot non_cilia -----------------------------------------------

close3d()
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  #defines the two scenes
#define the size of the rgl window, the view and zoom
par3d(windowRect = c(0, 0, 2400, 700))


#plot aboral view
plot3d(non_cilia,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[4],
       WithConnectors = F, WithNodes = F)
for (object in objects) {
  plot3d(object, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.05, col = Okabe_Ito[8])
}

nview3d("anterior", extramat = rotationMatrix(1.05, 250, -200, 1000))
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)

#plot lateral view of Sagittal plane
plot3d(non_cilia,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[4],
       WithConnectors = F, WithNodes = F)
for (object in objects) {
  plot3d(object, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.05, col = Okabe_Ito[8])
}

nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)

#plot lateral view of Tentacular plane
plot3d(non_cilia,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[4],
       WithConnectors = F, WithNodes = F)
for (object in objects) {
  plot3d(object, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.05, col = Okabe_Ito[8])
}

nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))
par3d(zoom=0.61)


#make a snapshot
rgl.snapshot("manuscript/pictures/non_cilia.png")

close3d()

# assemble figure -------------------------------------------------------------

#read pics
panel_balancer <- ggdraw() + draw_image(readPNG("manuscript/pictures/balancer.png")) +
  draw_label("balancer cells", x = 0.5, y = 0.95, size = 8.5, fontface="bold", hjust = 0.5) +
  draw_label("aboral view", x = 0.01, y = 0.86, color="black", size = 6, fontface="plain", hjust = 0) +
  draw_label("lateral view of PA plane", x = 0.33, y = 0.86, color="black", size = 6, fontface="plain", hjust = 0) +
  draw_label("lateral view of TA plane", x = 0.69, y = 0.86, color="black", size = 6, fontface="plain", hjust = 0) +
  draw_line(x = c(0.85, 0.95), y = c(0.1, 0.1), color = "black", size = 0.5) +
  draw_label(expression(paste("25 ", mu, "m")), x = 0.9, y = 0.14, color = "black", size = 7, hjust = 0.5) +
  draw_label("TA", x = 0.325, y = 0.16, size = 6, color = "black", hjust = 0.5) +
  geom_segment(aes(x = 0.25, y = 0.16, xend = 0.31, yend = 0.16), data = arrow_TA, color = "black", 
               arrow = arrow(ends = "both", type = "closed", length = unit(0.1,"cm")),
               lineend = "butt",
               linejoin = "mitre",
               arrow.fill = "black", size = 0.175) +
  draw_label("PA", x = 0.28, y = 0.05, size = 6, color = "black", hjust = 0.5) +
  geom_segment(aes(x = 0.28, y = 0.08, xend = 0.28, yend = 0.24), data = arrow_PA, color = "black", 
               arrow = arrow(ends = "both", type = "closed", length = unit(0.1,"cm")),
               lineend = "butt",
               linejoin = "mitre",
               arrow.fill = "black", size = 0.175) +
  draw_label("A", x = 0.66, y = 0.25, size = 6, color = "black", hjust = 0.5) +
  draw_label("O", x = 0.66, y = 0.05, size = 6, color = "black", hjust = 0.5) +
  geom_segment(aes(x = 0.66, y = 0.09, xend = 0.66, yend = 0.21), data = arrow_AO, color = "black", 
               arrow = arrow(ends = "both", type = "closed", length = unit(0.1,"cm")),
               lineend = "butt",
               linejoin = "mitre",
               arrow.fill = "black", size = 0.175)
 

panel_bridge <- ggdraw() + draw_image(readPNG("manuscript/pictures/bridge.png")) +
  draw_label("bridge cells", x = 0.5, y = 0.95, size = 8, fontface="bold", hjust = 0.5)
panel_bristle <- ggdraw() + draw_image(readPNG("manuscript/pictures/bristle.png")) +
  draw_label("bristle cells", x = 0.5, y = 0.95, size = 8, fontface="bold", hjust = 0.5)
panel_dome <- ggdraw() + draw_image(readPNG("manuscript/pictures/dome.png")) +
  draw_label("dome cells", x = 0.5, y = 0.95, size = 8, fontface="bold", hjust = 0.5)
panel_groove <- ggdraw() + draw_image(readPNG("manuscript/pictures/groove.png")) +
  draw_label("ciliated groove cells", x = 0.5, y = 0.95, size = 8, fontface="bold", hjust = 0.5)
panel_intramulticilia <- ggdraw() + draw_image(readPNG("manuscript/pictures/intramulticilia.png")) +
  draw_label("intra-multiciliated cells", x = 0.5, y = 0.95, size = 8, fontface="bold", hjust = 0.5)
panel_lamellate <- ggdraw() + draw_image(readPNG("manuscript/pictures/lamellate.png")) +
  draw_label("lamellate bodis", x = 0.5, y = 0.95, size = 8, fontface="bold", hjust = 0.5)
panel_lithocyte <- ggdraw() + draw_image(readPNG("manuscript/pictures/lithocyte.png")) +
  draw_label("lithocytes", x = 0.5, y = 0.95, size = 8, fontface="bold", hjust = 0.5)
panel_neuron <- ggdraw() + draw_image(readPNG("manuscript/pictures/neuron.png")) +
  draw_label("synaptic neurons", x = 0.5, y = 0.95, size = 8, fontface="bold", hjust = 0.5)
panel_plumose <- ggdraw() + draw_image(readPNG("manuscript/pictures/plumose.png")) +
  draw_label("plumose cells", x = 0.5, y = 0.95, size = 8, fontface="bold", hjust = 0.5)
panel_dense_vesicle <- ggdraw() + draw_image(readPNG("manuscript/pictures/dense_vesicle.png")) +
  draw_label("dense vesicle cells", x = 0.5, y = 0.95, size = 8, fontface="bold", hjust = 0.5)
panel_monocilia <- ggdraw() + draw_image(readPNG("manuscript/pictures/monocilia.png")) +
  draw_label("monociliated cells", x = 0.5, y = 0.95, size = 8, fontface="bold", hjust = 0.5)
panel_bicilia <- ggdraw() + draw_image(readPNG("manuscript/pictures/bicilia.png")) +
  draw_label("biciliated cells", x = 0.5, y = 0.95, size = 8, fontface="bold", hjust = 0.5)
panel_non_cilia <- ggdraw() + draw_image(readPNG("manuscript/pictures/non_cilia.png")) +
  draw_label("non-ciliated cells", x = 0.5, y = 0.95, size = 8, fontface="bold", hjust = 0.5)

layout <- "
ABC
DEF
GHI
JK#
LMN
"

Figure1 <- panel_balancer + panel_bridge + panel_bristle +
  panel_dome + panel_groove + panel_intramulticilia + panel_lamellate +
  panel_lithocyte + panel_neuron + panel_plumose +
  panel_dense_vesicle + panel_monocilia + panel_bicilia + panel_non_cilia +
  patchwork::plot_layout(design = layout, 
                         heights = c(1,1,1,1,1), 
                         widths = c(1,1,1)) + 
  patchwork::plot_annotation(tag_levels = "A") &  
  ggplot2::theme(plot.tag = element_text(size = 12, 
                                         face='plain', color='black'))


ggsave("manuscript/figures/Fig2_sup1.png", limitsize = FALSE, 
       units = c("px"), Figure1, width = 3800, height = 2400, bg='white')  


ggsave("manuscript/figures/Fig2_sup1.pdf", limitsize = FALSE, 
       units = c("px"), Figure1, width = 3800, height = 2400) 






















img_bridge_aboral <- readPNG("manuscript/pictures/bridge_aboral_view.png")
img_bridge_tentacular <- readPNG("manuscript/pictures/bridge_tentacular_plane.png")
img_bridge_sagittal <- readPNG("manuscript/pictures/bridge_sagittal_plane.png")
img_bristle_aboral <- readPNG("manuscript/pictures/bristle_aboral_view.png")
img_bristle_tentacular <- readPNG("manuscript/pictures/bristle_tentacular_plane.png")
img_bristle_sagittal <- readPNG("manuscript/pictures/bristle_sagittal_plane.png")
img_dome_aboral <- readPNG("manuscript/pictures/dome_aboral_view.png")
img_dome_tentacular <- readPNG("manuscript/pictures/dome_tentacular_plane.png")
img_dome_sagittal <- readPNG("manuscript/pictures/dome_sagittal_plane.png")
img_groove_aboral <- readPNG("manuscript/pictures/groove_aboral_view.png")
img_groove_tentacular <- readPNG("manuscript/pictures/groove_tentacular_plane.png")
img_groove_sagittal <- readPNG("manuscript/pictures/groove_sagittal_plane.png")
img_intramulticilia_aboral <- readPNG("manuscript/pictures/intramulticilia_aboral_view.png")
img_intramulticilia_tentacular <- readPNG("manuscript/pictures/intramulticilia_tentacular_plane.png")
img_intramulticilia_sagittal <- readPNG("manuscript/pictures/intramulticilia_sagittal_plane.png")
img_lamellate_aboral <- readPNG("manuscript/pictures/lamellate_aboral_view.png")
img_lamellate_tentacular <- readPNG("manuscript/pictures/lamellate_tentacular_plane.png")
img_lamellate_sagittal <- readPNG("manuscript/pictures/lamellate_sagittal_plane.png")
img_lithocyte_aboral <- readPNG("manuscript/pictures/lithocyte_aboral_view.png")
img_lithocyte_tentacular <- readPNG("manuscript/pictures/lithocyte_tentacular_plane.png")
img_lithocyte_sagittal <- readPNG("manuscript/pictures/lithocyte_sagittal_plane.png")
img_neuron_aboral <- readPNG("manuscript/pictures/neuron_aboral_view.png")
img_neuron_tentacular <- readPNG("manuscript/pictures/neuron_tentacular_plane.png")
img_neuron_sagittal <- readPNG("manuscript/pictures/neuron_sagittal_plane.png")
img_plumose_aboral <- readPNG("manuscript/pictures/plumose_aboral_view.png")
img_plumose_tentacular <- readPNG("manuscript/pictures/plumose_tentacular_plane.png")
img_plumose_sagittal <- readPNG("manuscript/pictures/plumose_sagittal_plane.png")
img_dense_vesicle_aboral <- readPNG("manuscript/pictures/dense_vesicle_aboral_view.png")
img_dense_vesicle_tentacular <- readPNG("manuscript/pictures/dense_vesicle_tentacular_plane.png")
img_dense_vesicle_sagittal <- readPNG("manuscript/pictures/dense_vesicle_sagittal_plane.png")
img_monocilia_aboral <- readPNG("manuscript/pictures/monocilia_aboral_view.png")
img_monocilia_tentacular <- readPNG("manuscript/pictures/monocilia_tentacular_plane.png")
img_monocilia_sagittal <- readPNG("manuscript/pictures/monocilia_sagittal_plane.png")
img_bicilia_aboral <- readPNG("manuscript/pictures/bicilia_aboral_view.png")
img_bicilia_tentacular <- readPNG("manuscript/pictures/bicilia_tentacular_plane.png")
img_bicilia_sagittal <- readPNG("manuscript/pictures/bicilia_sagittal_plane.png")
img_non_cilia_aboral <- readPNG("manuscript/pictures/non_cilia_aboral_view.png")
img_non_cilia_tentacular <- readPNG("manuscript/pictures/non_cilia_tentacular_plane.png")
img_non_cilia_sagittal <- readPNG("manuscript/pictures/non_cilia_sagittal_plane.png")


#make panels

panel_balancer <- ggdraw() + draw_image(img_balancer_aboral) +
  draw_label("aboral view", x = 0.3, y = 0.62, size = 6, fontface = "plain", hjust = 0.5) +
  ggdraw() + draw_image(img_balancer_sagittal) +
  draw_label("saggital plane", x = 0.3, y = 0.62, size = 6, fontface = "plain", hjust = 0.5) + 
  ggdraw() + draw_image(img_balancer_tentacular) +
  draw_label("tentacular plane", x = 0.3, y = 0.62, size = 6, fontface = "plain", hjust = 0.5) +
  draw_line(x = c(0.58, 0.82), y = c(0.395, 0.395), color = "black", size = 0.5) +
  draw_label(expression(paste("20 ", mu, "m")), x = 0.7, y = 0.415, color = "black", size = 6, hjust = 0.5)

panel_bridge <- ggdraw() + draw_image(img_bridge_aboral) +
  ggdraw() + draw_image(img_bridge_sagittal) + 
  ggdraw() + draw_image(img_bridge_tentacular)

panel_bristle <- ggdraw() + draw_image(img_bristle_aboral) +
  ggdraw() + draw_image(img_bristle_sagittal) + 
  ggdraw() + draw_image(img_bristle_tentacular)

panel_dome <- ggdraw() + draw_image(img_dome_aboral) +
  ggdraw() + draw_image(img_dome_sagittal) + 
  ggdraw() + draw_image(img_dome_tentacular)

panel_groove <- ggdraw() + draw_image(img_groove_aboral) +
  ggdraw() + draw_image(img_groove_sagittal) + 
  ggdraw() + draw_image(img_groove_tentacular)

panel_intramulticilia <- ggdraw() + draw_image(img_intramulticilia_aboral) +
  ggdraw() + draw_image(img_intramulticilia_sagittal) + 
  ggdraw() + draw_image(img_intramulticilia_tentacular)

panel_lamellate <- ggdraw() + draw_image(img_lamellate_aboral) +
  ggdraw() + draw_image(img_lamellate_sagittal) + 
  ggdraw() + draw_image(img_lamellate_tentacular)

panel_lithocyte <- ggdraw() + draw_image(img_lithocyte_aboral) +
  ggdraw() + draw_image(img_lithocyte_sagittal) + 
  ggdraw() + draw_image(img_lithocyte_tentacular)

panel_neuron <- ggdraw() + draw_image(img_neuron_aboral) +
  ggdraw() + draw_image(img_neuron_sagittal) + 
  ggdraw() + draw_image(img_neuron_tentacular)

panel_plumose <- ggdraw() + draw_image(img_plumose_aboral) +
  ggdraw() + draw_image(img_plumose_sagittal) + 
  ggdraw() + draw_image(img_plumose_tentacular)

panel_dense_vesicle <- ggdraw() + draw_image(img_dense_vesicle_aboral) +
  ggdraw() + draw_image(img_dense_vesicle_sagittal) + 
  ggdraw() + draw_image(img_dense_vesicle_tentacular)

panel_monocilia <- ggdraw() + draw_image(img_monocilia_aboral) +
  ggdraw() + draw_image(img_monocilia_sagittal) + 
  ggdraw() + draw_image(img_monocilia_tentacular)

panel_bicilia <- ggdraw() + draw_image(img_bicilia_aboral) +
  ggdraw() + draw_image(img_bicilia_sagittal) + 
  ggdraw() + draw_image(img_bicilia_tentacular)

panel_non_cilia <- ggdraw() + draw_image(img_non_cilia_aboral) +
  ggdraw() + draw_image(img_non_cilia_sagittal) + 
  ggdraw() + draw_image(img_non_cilia_tentacular)




# define layout with textual representation for pathchwork assembly of figure
layout <- "
ABC
DEF
GHI
JKL
M##
"

Figure1 <-
  panel_balancer + panel_bridge + panel_bristle + panel_dome + 
  panel_groove + panel_intramulticilia + panel_lamellate + 
  panel_lithocyte + panel_neuron + panel_plumose + 
  panel_dense_vesicle + panel_monocilia + panel_bicilia + panel_non_cilia +
  plot_layout(design = layout) +
  patchwork::plot_layout(design = layout, 
                         heights = c(1,1,1,1,1), 
                         widths = c(1,1,1)) +
  patchwork::plot_annotation(tag_levels = "A") &
  ggplot2::theme(plot.tag = element_text(size = 12, face='plain', color='black'))


ggsave("manuscript/figures/Fig2_sup2.png", limitsize = FALSE, 
       units = c("px"), Figure1, width = 3200, height = 3200, bg = "white")


ggsave(
  "manuscript/figures/Fig2_sup2.pdf",
  limitsize = FALSE,
  units = c("px"),
  Figure1,
  width = 3300,
  height = 1600
)










# assemble figure -------------------------------------------------------------

#read pics
panel_balancer <- ggdraw() + draw_image(readPNG("manuscript/pictures/balancer.png")) +
  draw_label("balancer cells", x = 0.5, y = 0.95, size = 8.5, fontface="bold", hjust = 0.5) +
  draw_label("aboral view", x = 0.01, y = 0.86, color="black", size = 6, fontface="plain", hjust = 0) +
  draw_label("lateral view of PA plane", x = 0.33, y = 0.86, color="black", size = 6, fontface="plain", hjust = 0) +
  draw_label("lateral view of TA plane", x = 0.69, y = 0.86, color="black", size = 6, fontface="plain", hjust = 0) +
  draw_line(x = c(0.85, 0.95), y = c(0.1, 0.1), color = "black", size = 0.5) +
  draw_label(expression(paste("25 ", mu, "m")), x = 0.9, y = 0.14, color = "black", size = 7, hjust = 0.5) +
  draw_label("TA", x = 0.325, y = 0.16, size = 6, color = "black", hjust = 0.5) +
  geom_segment(aes(x = 0.25, y = 0.16, xend = 0.31, yend = 0.16), data = arrow_TA, color = "black", 
               arrow = arrow(ends = "both", type = "closed", length = unit(0.1,"cm")),
               lineend = "butt",
               linejoin = "mitre",
               arrow.fill = "black", size = 0.175) +
  draw_label("PA", x = 0.28, y = 0.05, size = 6, color = "black", hjust = 0.5) +
  geom_segment(aes(x = 0.28, y = 0.08, xend = 0.28, yend = 0.24), data = arrow_PA, color = "black", 
               arrow = arrow(ends = "both", type = "closed", length = unit(0.1,"cm")),
               lineend = "butt",
               linejoin = "mitre",
               arrow.fill = "black", size = 0.175) +
  draw_label("A", x = 0.66, y = 0.25, size = 6, color = "black", hjust = 0.5) +
  draw_label("O", x = 0.66, y = 0.05, size = 6, color = "black", hjust = 0.5) +
  geom_segment(aes(x = 0.66, y = 0.09, xend = 0.66, yend = 0.21), data = arrow_AO, color = "black", 
               arrow = arrow(ends = "both", type = "closed", length = unit(0.1,"cm")),
               lineend = "butt",
               linejoin = "mitre",
               arrow.fill = "black", size = 0.175)


img_bridge_aboral <- readPNG("manuscript/pictures/bridge_aboral_view.png")
img_bridge_tentacular <- readPNG("manuscript/pictures/bridge_tentacular_plane.png")
img_bridge_sagittal <- readPNG("manuscript/pictures/bridge_sagittal_plane.png")
img_bristle_aboral <- readPNG("manuscript/pictures/bristle_aboral_view.png")
img_bristle_tentacular <- readPNG("manuscript/pictures/bristle_tentacular_plane.png")
img_bristle_sagittal <- readPNG("manuscript/pictures/bristle_sagittal_plane.png")
img_dome_aboral <- readPNG("manuscript/pictures/dome_aboral_view.png")
img_dome_tentacular <- readPNG("manuscript/pictures/dome_tentacular_plane.png")
img_dome_sagittal <- readPNG("manuscript/pictures/dome_sagittal_plane.png")
img_groove_aboral <- readPNG("manuscript/pictures/groove_aboral_view.png")
img_groove_tentacular <- readPNG("manuscript/pictures/groove_tentacular_plane.png")
img_groove_sagittal <- readPNG("manuscript/pictures/groove_sagittal_plane.png")
img_intramulticilia_aboral <- readPNG("manuscript/pictures/intramulticilia_aboral_view.png")
img_intramulticilia_tentacular <- readPNG("manuscript/pictures/intramulticilia_tentacular_plane.png")
img_intramulticilia_sagittal <- readPNG("manuscript/pictures/intramulticilia_sagittal_plane.png")
img_lamellate_aboral <- readPNG("manuscript/pictures/lamellate_aboral_view.png")
img_lamellate_tentacular <- readPNG("manuscript/pictures/lamellate_tentacular_plane.png")
img_lamellate_sagittal <- readPNG("manuscript/pictures/lamellate_sagittal_plane.png")
img_lithocyte_aboral <- readPNG("manuscript/pictures/lithocyte_aboral_view.png")
img_lithocyte_tentacular <- readPNG("manuscript/pictures/lithocyte_tentacular_plane.png")
img_lithocyte_sagittal <- readPNG("manuscript/pictures/lithocyte_sagittal_plane.png")
img_neuron_aboral <- readPNG("manuscript/pictures/neuron_aboral_view.png")
img_neuron_tentacular <- readPNG("manuscript/pictures/neuron_tentacular_plane.png")
img_neuron_sagittal <- readPNG("manuscript/pictures/neuron_sagittal_plane.png")
img_plumose_aboral <- readPNG("manuscript/pictures/plumose_aboral_view.png")
img_plumose_tentacular <- readPNG("manuscript/pictures/plumose_tentacular_plane.png")
img_plumose_sagittal <- readPNG("manuscript/pictures/plumose_sagittal_plane.png")
img_dense_vesicle_aboral <- readPNG("manuscript/pictures/dense_vesicle_aboral_view.png")
img_dense_vesicle_tentacular <- readPNG("manuscript/pictures/dense_vesicle_tentacular_plane.png")
img_dense_vesicle_sagittal <- readPNG("manuscript/pictures/dense_vesicle_sagittal_plane.png")
img_monocilia_aboral <- readPNG("manuscript/pictures/monocilia_aboral_view.png")
img_monocilia_tentacular <- readPNG("manuscript/pictures/monocilia_tentacular_plane.png")
img_monocilia_sagittal <- readPNG("manuscript/pictures/monocilia_sagittal_plane.png")
img_bicilia_aboral <- readPNG("manuscript/pictures/bicilia_aboral_view.png")
img_bicilia_tentacular <- readPNG("manuscript/pictures/bicilia_tentacular_plane.png")
img_bicilia_sagittal <- readPNG("manuscript/pictures/bicilia_sagittal_plane.png")
img_non_cilia_aboral <- readPNG("manuscript/pictures/non_cilia_aboral_view.png")
img_non_cilia_tentacular <- readPNG("manuscript/pictures/non_cilia_tentacular_plane.png")
img_non_cilia_sagittal <- readPNG("manuscript/pictures/non_cilia_sagittal_plane.png")


#make panels

panel_balancer <- ggdraw() + draw_image(img_balancer_aboral) +
  draw_label("aboral view", x = 0.3, y = 0.62, size = 6, fontface = "plain", hjust = 0.5) +
  ggdraw() + draw_image(img_balancer_sagittal) +
  draw_label("saggital plane", x = 0.3, y = 0.62, size = 6, fontface = "plain", hjust = 0.5) + 
  ggdraw() + draw_image(img_balancer_tentacular) +
  draw_label("tentacular plane", x = 0.3, y = 0.62, size = 6, fontface = "plain", hjust = 0.5) +
  draw_line(x = c(0.58, 0.82), y = c(0.395, 0.395), color = "black", size = 0.5) +
  draw_label(expression(paste("20 ", mu, "m")), x = 0.7, y = 0.415, color = "black", size = 6, hjust = 0.5)

panel_bridge <- ggdraw() + draw_image(img_bridge_aboral) +
  ggdraw() + draw_image(img_bridge_sagittal) + 
  ggdraw() + draw_image(img_bridge_tentacular)

panel_bristle <- ggdraw() + draw_image(img_bristle_aboral) +
  ggdraw() + draw_image(img_bristle_sagittal) + 
  ggdraw() + draw_image(img_bristle_tentacular)

panel_dome <- ggdraw() + draw_image(img_dome_aboral) +
  ggdraw() + draw_image(img_dome_sagittal) + 
  ggdraw() + draw_image(img_dome_tentacular)

panel_groove <- ggdraw() + draw_image(img_groove_aboral) +
  ggdraw() + draw_image(img_groove_sagittal) + 
  ggdraw() + draw_image(img_groove_tentacular)

panel_intramulticilia <- ggdraw() + draw_image(img_intramulticilia_aboral) +
  ggdraw() + draw_image(img_intramulticilia_sagittal) + 
  ggdraw() + draw_image(img_intramulticilia_tentacular)

panel_lamellate <- ggdraw() + draw_image(img_lamellate_aboral) +
  ggdraw() + draw_image(img_lamellate_sagittal) + 
  ggdraw() + draw_image(img_lamellate_tentacular)

panel_lithocyte <- ggdraw() + draw_image(img_lithocyte_aboral) +
  ggdraw() + draw_image(img_lithocyte_sagittal) + 
  ggdraw() + draw_image(img_lithocyte_tentacular)

panel_neuron <- ggdraw() + draw_image(img_neuron_aboral) +
  ggdraw() + draw_image(img_neuron_sagittal) + 
  ggdraw() + draw_image(img_neuron_tentacular)

panel_plumose <- ggdraw() + draw_image(img_plumose_aboral) +
  ggdraw() + draw_image(img_plumose_sagittal) + 
  ggdraw() + draw_image(img_plumose_tentacular)

panel_dense_vesicle <- ggdraw() + draw_image(img_dense_vesicle_aboral) +
  ggdraw() + draw_image(img_dense_vesicle_sagittal) + 
  ggdraw() + draw_image(img_dense_vesicle_tentacular)

panel_monocilia <- ggdraw() + draw_image(img_monocilia_aboral) +
  ggdraw() + draw_image(img_monocilia_sagittal) + 
  ggdraw() + draw_image(img_monocilia_tentacular)

panel_bicilia <- ggdraw() + draw_image(img_bicilia_aboral) +
  ggdraw() + draw_image(img_bicilia_sagittal) + 
  ggdraw() + draw_image(img_bicilia_tentacular)

panel_non_cilia <- ggdraw() + draw_image(img_non_cilia_aboral) +
  ggdraw() + draw_image(img_non_cilia_sagittal) + 
  ggdraw() + draw_image(img_non_cilia_tentacular)


