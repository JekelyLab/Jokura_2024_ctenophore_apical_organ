# Code to generate Figure 1 of the Jokura et al 2024 Ctenophore apical organ connectome paper

# source packages and functions -----------------
source("analysis/scripts/packages_and_functions.R")

# load cell type ----------------------------------

balancer <-
  nlapply(read.neurons.catmaid("celltype:balancer", pid = 35),
          function(x)
            smooth_neuron(x, sigma = 1000))

bridge <-
  nlapply(read.neurons.catmaid("celltype:bridge", pid = 35),
          function(x)
            smooth_neuron(x, sigma = 1000))

bristle <-
  nlapply(read.neurons.catmaid("celltype:bristle", pid = 35),
          function(x)
            smooth_neuron(x, sigma = 1000))

dome <-
  nlapply(read.neurons.catmaid("celltype:dome", pid = 35),
          function(x)
            smooth_neuron(x, sigma = 1000))

groove <-
  nlapply(read.neurons.catmaid("celltype:groove", pid = 35),
          function(x)
            smooth_neuron(x, sigma = 1000))

intramulticilia <-
  nlapply(read.neurons.catmaid("celltype:intra-multi-ciliated", pid = 35),
          function(x)
            smooth_neuron(x, sigma = 1000))

lamellate <-
  nlapply(read.neurons.catmaid("celltype:lamellate", pid = 35),
          function(x)
            smooth_neuron(x, sigma = 1000))

lithocyte <-
  nlapply(read.neurons.catmaid("celltype:lithocyte", pid = 35),
          function(x)
            smooth_neuron(x, sigma = 1000))

neuron <-
  nlapply(read.neurons.catmaid("celltype:neuron", pid = 35),
          function(x)
            smooth_neuron(x, sigma = 1000))

plumose <-
  nlapply(read.neurons.catmaid("celltype:plumose", pid = 35),
          function(x)
            smooth_neuron(x, sigma = 1000))

dense_vesicle <-
  nlapply(read.neurons.catmaid("celltype:dense_vesicle", pid = 35),
          function(x)
            smooth_neuron(x, sigma = 1000))

monocilia <-
  nlapply(read.neurons.catmaid("celltype:monociliated", pid = 35),
          function(x)
            smooth_neuron(x, sigma = 1000))

bicilia <-
  nlapply(read.neurons.catmaid("celltype:biciliated", pid = 35),
          function(x)
            smooth_neuron(x, sigma = 1000))

non_cilia <-
  nlapply(read.neurons.catmaid("celltype:nonciliated", pid = 35),
          function(x)
            smooth_neuron(x, sigma = 1000))


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


# plot balancer -----------------------------------------------

plot_background()

plot3d(balancer,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)

plot3d(bridge, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(bristle, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(dome, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(groove, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(intramulticilia, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(lamellate, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(lithocyte, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(neuron, 
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
rgl.snapshot("manuscript/pictures/balancer_aboral_view.png")


#lateral view of Tentacular plane
nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))
rgl.snapshot("manuscript/pictures/balancer_tentacular_plane.png")

#lateral view of Sagittal plane        
nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
rgl.snapshot("manuscript/pictures/balancer_sagittal_plane.png")

close3d()



# plot bridge -----------------------------------------------

plot_background()

plot3d(bridge,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)

plot3d(balancer, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(bristle, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(dome, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(groove, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(intramulticilia, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(lamellate, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(lithocyte, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(neuron, 
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
rgl.snapshot("manuscript/pictures/bridge_aboral_view.png")


#lateral view of Tentacular plane
nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))
rgl.snapshot("manuscript/pictures/bridge_tentacular_plane.png")

#lateral view of Sagittal plane        
nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
rgl.snapshot("manuscript/pictures/bridge_sagittal_plane.png")

close3d()

# plot bristle -----------------------------------------------

plot_background()

plot3d(bristle,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)

plot3d(balancer, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(bridge, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(dome, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(groove, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(intramulticilia, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(lamellate, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(lithocyte, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(neuron, 
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
rgl.snapshot("manuscript/pictures/bristle_aboral_view.png")


#lateral view of Tentacular plane
nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))
rgl.snapshot("manuscript/pictures/bristle_tentacular_plane.png")

#lateral view of Sagittal plane        
nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
rgl.snapshot("manuscript/pictures/bristle_sagittal_plane.png")

close3d()


# plot dome -----------------------------------------------

plot_background()

plot3d(dome,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)

plot3d(balancer, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(bridge, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(bristle, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(groove, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(intramulticilia, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(lamellate, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(lithocyte, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(neuron, 
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
rgl.snapshot("manuscript/pictures/dome_aboral_view.png")


#lateral view of Tentacular plane
nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))
rgl.snapshot("manuscript/pictures/dome_tentacular_plane.png")

#lateral view of Sagittal plane        
nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
rgl.snapshot("manuscript/pictures/dome_sagittal_plane.png")

close3d()

# plot groove -----------------------------------------------

plot_background()

plot3d(groove,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)

plot3d(balancer, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(bridge, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(bristle, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(dome, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(intramulticilia, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(lamellate, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(lithocyte, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(neuron, 
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
rgl.snapshot("manuscript/pictures/groove_aboral_view.png")


#lateral view of Tentacular plane
nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))
rgl.snapshot("manuscript/pictures/groove_tentacular_plane.png")

#lateral view of Sagittal plane        
nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
rgl.snapshot("manuscript/pictures/groove_sagittal_plane.png")

close3d()


# plot intramulticilia -----------------------------------------------

plot_background()

plot3d(intramulticilia,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)

plot3d(balancer, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(bridge, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(bristle, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(dome, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(groove, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(lamellate, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(lithocyte, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(neuron, 
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
rgl.snapshot("manuscript/pictures/intramulticilia_aboral_view.png")


#lateral view of Tentacular plane
nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))
rgl.snapshot("manuscript/pictures/intramulticilia_tentacular_plane.png")

#lateral view of Sagittal plane        
nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
rgl.snapshot("manuscript/pictures/intramulticilia_sagittal_plane.png")

close3d()

# plot lamellate -----------------------------------------------

plot_background()

plot3d(lamellate,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)

plot3d(balancer, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(bridge, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(bristle, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(dome, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(groove, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(intramulticilia, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(lithocyte, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(neuron, 
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
rgl.snapshot("manuscript/pictures/lamellate_aboral_view.png")


#lateral view of Tentacular plane
nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))
rgl.snapshot("manuscript/pictures/lamellate_tentacular_plane.png")

#lateral view of Sagittal plane        
nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
rgl.snapshot("manuscript/pictures/lamellate_sagittal_plane.png")

close3d()


# plot lithocyte -----------------------------------------------

plot_background()

plot3d(lithocyte,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)

plot3d(balancer, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(bridge, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(bristle, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(dome, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(groove, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(intramulticilia, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(lamellate, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(neuron, 
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
rgl.snapshot("manuscript/pictures/lithocyte_aboral_view.png")


#lateral view of Tentacular plane
nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))
rgl.snapshot("manuscript/pictures/lithocyte_tentacular_plane.png")

#lateral view of Sagittal plane        
nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
rgl.snapshot("manuscript/pictures/lithocyte_sagittal_plane.png")

close3d()


# plot neuron -----------------------------------------------

plot_background()

plot3d(neuron,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)

plot3d(balancer, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(bridge, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(bristle, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(dome, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(groove, 
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
rgl.snapshot("manuscript/pictures/neuron_aboral_view.png")


#lateral view of Tentacular plane
nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))
rgl.snapshot("manuscript/pictures/neuron_tentacular_plane.png")

#lateral view of Sagittal plane        
nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
rgl.snapshot("manuscript/pictures/neuron_sagittal_plane.png")

close3d()


# plot plumose -----------------------------------------------

plot_background()

plot3d(plumose,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)

plot3d(balancer, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(bridge, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(bristle, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(dome, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(groove, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(intramulticilia, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(lamellate, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(lithocyte, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(neuron, 
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
rgl.snapshot("manuscript/pictures/plumose_aboral_view.png")

#lateral view of Tentacular plane
nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))
rgl.snapshot("manuscript/pictures/plumose_tentacular_plane.png")

#lateral view of Sagittal plane        
nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
rgl.snapshot("manuscript/pictures/plumose_sagittal_plane.png")

close3d()


# plot dense_vesicle -----------------------------------------------

plot_background()

plot3d(dense_vesicle,
       soma = T, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)

plot3d(balancer, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(bridge, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(bristle, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(dome, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(groove, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(intramulticilia, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(lamellate, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(lithocyte, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(neuron, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(plumose, 
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
rgl.snapshot("manuscript/pictures/dense_vesicle_aboral_view.png")

#lateral view of Tentacular plane
nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))
rgl.snapshot("manuscript/pictures/dense_vesicle_tentacular_plane.png")

#lateral view of Sagittal plane        
nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
rgl.snapshot("manuscript/pictures/dense_vesicle_sagittal_plane.png")

close3d()


# plot monocilia -----------------------------------------------

plot_background()

plot3d(monocilia,
       soma = T, lwd = 1, add = T, alpha = 0.3, col = Okabe_Ito[3],
       WithConnectors = F, WithNodes = F)

plot3d(balancer, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(bridge, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(bristle, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(dome, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(groove, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(intramulticilia, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(lamellate, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(lithocyte, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(neuron, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(plumose, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(dense_vesicle, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])

plot3d(bicilia, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(non_cilia, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])

#texts3d(58000,51000,5000, "INRGW", cex = 3, col = "#56B4E9")

#aboral view
nview3d("anterior", extramat = rotationMatrix(1.05, 250, -200, 1000))
rgl.snapshot("manuscript/pictures/monocilia_aboral_view.png")

#lateral view of Tentacular plane
nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))
rgl.snapshot("manuscript/pictures/monocilia_tentacular_plane.png")

#lateral view of Sagittal plane        
nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
rgl.snapshot("manuscript/pictures/monocilia_sagittal_plane.png")

close3d()

# plot bicilia -----------------------------------------------

plot_background()

plot3d(bicilia,
       soma = T, lwd = 1, add = T, alpha = 0.3, col = Okabe_Ito[7],
       WithConnectors = F, WithNodes = F)

plot3d(balancer, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(bridge, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(bristle, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(dome, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(groove, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(intramulticilia, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(lamellate, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(lithocyte, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(neuron, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(plumose, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(dense_vesicle, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])

plot3d(monocilia, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(non_cilia, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])

#texts3d(58000,51000,5000, "INRGW", cex = 3, col = "#56B4E9")

#aboral view
nview3d("anterior", extramat = rotationMatrix(1.05, 250, -200, 1000))
rgl.snapshot("manuscript/pictures/bicilia_aboral_view.png")

#lateral view of Tentacular plane
nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))
rgl.snapshot("manuscript/pictures/bicilia_tentacular_plane.png")

#lateral view of Sagittal plane        
nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
rgl.snapshot("manuscript/pictures/bicilia_sagittal_plane.png")

close3d()

# plot non_cilia -----------------------------------------------

plot_background()

plot3d(non_cilia,
       soma = T, lwd = 1, add = T, alpha = 0.3, col = Okabe_Ito[1],
       WithConnectors = F, WithNodes = F)

plot3d(balancer, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(bridge, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(bristle, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(dome, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(groove, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(intramulticilia, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(lamellate, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(lithocyte, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(neuron, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(plumose, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(dense_vesicle, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])

plot3d(monocilia, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])
plot3d(non_cilia, 
       soma = T, lwd = 1, add = T, alpha = 0.05, col = Okabe_Ito[8])

#texts3d(58000,51000,5000, "INRGW", cex = 3, col = "#56B4E9")

#aboral view
nview3d("anterior", extramat = rotationMatrix(1.05, 250, -200, 1000))
rgl.snapshot("manuscript/pictures/non_cilia_aboral_view.png")

#lateral view of Tentacular plane
nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))
rgl.snapshot("manuscript/pictures/non_cilia_tentacular_plane.png")

#lateral view of Sagittal plane        
nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
rgl.snapshot("manuscript/pictures/non_cilia_sagittal_plane.png")

close3d()

Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
               "#0072B2", "#D55E00", "#CC79A7", "#000000")



# assemble figure -------------------------------------------------------------


#read pics

img_balancer_aboral <- readPNG("manuscript/pictures/balancer_aboral_view.png")
img_balancer_tentacular <- readPNG("manuscript/pictures/balancer_tentacular_plane.png")
img_balancer_sagittal <- readPNG("manuscript/pictures/balancer_sagittal_plane.png")
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


ggsave("manuscript/figures/Figure1.png", limitsize = FALSE, 
       units = c("px"), Figure1, width = 3200, height = 3200, bg = "white")


ggsave(
  "manuscript/figures/Figure1.pdf",
  limitsize = FALSE,
  units = c("px"),
  Figure1,
  width = 3300,
  height = 1600
)
