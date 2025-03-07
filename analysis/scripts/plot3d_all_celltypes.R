# Code to generate all cell types in 3d plot and 
#save the views from three different directions as png files

# source packages and functions -----------------
source("analysis/scripts/packages_and_functions.R")

# get celltypes
celltypes <- get_celltypes(35)
celltypes_list <- list()
for (celltype in celltypes) {
  celltypes_list[[celltype]] <- read.neurons.catmaid(paste("celltype:", celltype, sep = ""), pid = 35)
}

colour_palettes <- c(Okabe_Ito, bluepurple, oranges)

nopen3d() 
mfrow3d(1, 3)  #defines the two scenes
par3d(windowRect = c(0, 0, 2400, 700))


start.time <- Sys.time()
#plot aboral view
for (i in seq_along(celltypes_list)) {
  # using i because counter is needed for color palette
  # only use plot_multinucleated_cell for SSN, because it's much slower than normal plotting
  if (names(celltypes_list)[[i]] == "SSN") {
    plot_multinucleated_cell(celltypes_list[[i]], alpha = 0.6, color = colour_palettes[i])
  } else {
    plot3d(celltypes_list[[i]], alpha = 0.6, color = colour_palettes[i], soma = TRUE)
  }
  
  i = i + 1
}
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken


anterior()
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
sagittal()
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
tentacular()
par3d(zoom=0.61)

#make a snapshot
rgl.snapshot("manuscript/pictures/all_cells_3_views_alt.png")
close3d()

