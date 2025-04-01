source("analysis/scripts/packages_and_functions.R")

stats_synapse <- read.csv("analysis/data/stats_synapse.csv")

SSN_Q1Q2 <- read_smooth_neuron("SSN_Q1Q2")[[1]]
skid_Q1Q2 <- SSN_Q1Q2$skid
SSN_Q3Q4 <- read_smooth_neuron("SSN_Q3Q4")[[1]]
skid_Q3Q4 <- SSN_Q3Q4$skid
SSN_Q1Q2Q3Q4 <- read_smooth_neuron("SSN_Q1Q2Q3Q4")[[1]]
skid_Q1Q2Q3Q4 <- SSN_Q1Q2Q3Q4$skid

small_outgoing_syn <- stats_synapse %>% filter(skid == skid_Q1Q2 | skid == skid_Q3Q4) %>% filter(prepost == 0) %>% select(x, y, z)
big_outgoing_syn <- stats_synapse %>% filter(skid == skid_Q1Q2Q3Q4) %>% filter(prepost == 0) %>% select(x, y, z)

# 3d plot SSN Q1Q2 & Q3Q4 neuron -----------------------------------------------
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  
#define the size of the rgl window, the view and zoom
par3d(windowRect = c(0, 0, 1200, 350))
#par3d(windowRect = c(0, 0, 2400, 700))

#plot aboral view
plot_multinucleated_cell(SSN_Q1Q2,
                         lwd = 1, alpha = 1, col = Okabe_Ito[6])

plot_multinucleated_cell(SSN_Q3Q4,
                         lwd = 1, alpha = 1, col = Okabe_Ito[7])

plot3d(outline,
       add = T, alpha = 0.05, col = "grey50"
)
plot3d(small_outgoing_syn, add = T, size = 5)

#aboral view
aboral()
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)


#plot lateral view of Sagittal plane
plot_multinucleated_cell(SSN_Q1Q2,
                         lwd = 1, alpha = 1, col = Okabe_Ito[6])

plot_multinucleated_cell(SSN_Q3Q4,
                         lwd = 1, alpha = 1, col = Okabe_Ito[7])

plot3d(outline,
       add = T, alpha = 0.05, col = "grey50")
plot3d(small_outgoing_syn, add = T, size = 5)


#lateral view of Sagittal plane
sagittal()
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)


#plot lateral view of Tentacular plane
plot_multinucleated_cell(SSN_Q1Q2,
                         lwd = 1, alpha = 1, col = Okabe_Ito[6])

plot_multinucleated_cell(SSN_Q3Q4,
                         lwd = 1, alpha = 1, col = Okabe_Ito[7])

plot3d(outline,
       add = T, alpha = 0.05, col = "grey50")
plot3d(small_outgoing_syn, add = T, size = 5)


#lateral view of Tentacular plane
tentacular()
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)

#make a snapshot to the working directory
rgl.snapshot("manuscript/pictures/synapses_ANN_small.png")


close3d()

# 3d plot SSN Q1Q2Q3Q4 neuron -----------------------------------------------

# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  
#define the size of the rgl window, the view and zoom
par3d(windowRect = c(0, 0, 1200, 350))

#plot aboral view
plot3d(bridge, lwd = 1, alpha = 1, col = Okabe_Ito[5])

plot3d(outline,
       add = T, alpha = 0.05, col = "grey50")
plot3d(big_outgoing_syn, add = T, size = 5)


#aboral view
aboral()
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)


#plot lateral view of Sagittal plane
plot_multinucleated_cell(SSN_Q1Q2Q3Q4,
                         lwd = 1, alpha = 1, col = Okabe_Ito[5])

plot3d(outline,
       add = T, alpha = 0.05, col = "grey50")

plot3d(big_outgoing_syn, add = T, size = 5)
#lateral view of Sagittal plane
sagittal()
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)

#plot lateral view of Tentacular plane
plot_multinucleated_cell(SSN_Q1Q2Q3Q4,
                         lwd = 1, alpha = 1, col = Okabe_Ito[5])

plot3d(outline,
       add = T, alpha = 0.05, col = "grey50"
)
plot3d(big_outgoing_syn, add = T, size = 5)

#lateral view of Tentacular plane
tentacular()
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)

#make a snapshot to the working directory
rgl.snapshot("mnuscript/pictures/synapses_ANN_big.png")
close3d()


##### plot bridge cells

bridge <- read_smooth_neuron("celltype:bridge")
skids_bridge <- names(bridge)

bridge_outgoing_syn <- stats_synapse %>% filter(skid %in% skids_bridge) %>% filter(prepost == 0) %>% select(x, y, z)


# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  
#define the size of the rgl window, the view and zoom
par3d(windowRect = c(0, 0, 1200, 350))

#plot aboral view
plot3d(bridge, soma = T, lwd = 1, alpha = 1, col = Okabe_Ito[4])
plot3d(outline,
       add = T, alpha = 0.05, col = "grey50")
plot3d(bridge_outgoing_syn, add = T, size = 5)

#aboral view
aboral()
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)


#plot lateral view of Sagittal plane
plot3d(bridge, soma = T, lwd = 1, alpha = 1, col = Okabe_Ito[4])
plot3d(outline,
       add = T, alpha = 0.05, col = "grey50")
plot3d(bridge_outgoing_syn, add = T, size = 5)
#lateral view of Sagittal plane
sagittal()
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)

#plot lateral view of Tentacular plane
plot3d(bridge, soma = T, lwd = 1, alpha = 1, col = Okabe_Ito[4])
plot3d(outline,
       add = T, alpha = 0.05, col = "grey50")
plot3d(bridge_outgoing_syn, add = T, size = 5)

#lateral view of Tentacular plane
tentacular()
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)

#make a snapshot to the working directory
rgl.snapshot("pictures/synapses_bridge.png")
close3d()
