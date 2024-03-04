# Code to generate Figure balancer of the Jokura et al 2024 Ctenophore apical organ connectome paper

# source packages and functions ------------------------------------------------
source("analysis/scripts/packages_and_functions.R")

# load cell types from catmaid -------------------------------------------------

read_smooth_cell <- function(annotation){
  nlapply(read.neurons.catmaid(annotation, pid = 35),
          function(x)
            smooth_neuron(x, sigma = 500))
}

#balancer <- read_smooth_cell("celltype:balancer")
#bridge <- read_smooth_cell("celltype:bridge")
#bristle <- read_smooth_cell("celltype:bristle")
#dome <- read_smooth_cell("celltype:dome")
#groove <- read_smooth_cell("celltype:groove")
#intramulticilia <- read_smooth_cell("celltype:intra-multi-ciliated")
#lamellate <- read_smooth_cell("celltype:lamellate")
#lithocyte <- read_smooth_cell("celltype:lithocyte")
#cell <- read_smooth_cell("celltype:neuron")
#plumose <- read_smooth_cell("celltype:plumose")
#dense_vesicle <- read_smooth_cell("celltype:dense_vesicle")
#monocilia <- read_smooth_cell("celltype:monociliated")
#bicilia <- read_smooth_cell("celltype:biciliated")
#non_cilia <- read_smooth_cell("celltype:nonciliated")


balancer_Q1 <- read_smooth_cell(skids_by_2annotations("celltype:balancer", "Q1"))
balancer_Q2 <- read_smooth_cell(skids_by_2annotations("celltype:balancer", "Q2"))
balancer_Q3 <- read_smooth_cell(skids_by_2annotations("celltype:balancer", "Q3"))
balancer_Q4 <- read_smooth_cell(skids_by_2annotations("celltype:balancer", "Q4"))
intramulticilia <- read_smooth_cell("celltype:intra-multi-ciliated")
lamellate <- read_smooth_cell("celltype:lamellate")
bridge <- read_smooth_cell("celltype:bridge")
statolith <- read_smooth_cell("statolith")
Cgroove_tag <- read_smooth_cell("celltype:Cgroove-tag")
Cgroove_sag <- read_smooth_cell("celltype:Cgroove-sag")
lamellate_Q3 <- read_smooth_cell(skids_by_2annotations("celltype:lamellate", "Q3"))
intramulticilia_Q3 <- read_smooth_cell(skids_by_2annotations("celltype:intra-multi-ciliated", "Q3"))
Cgroove_tag_Q3 <- read_smooth_cell(skids_by_2annotations("celltype:Cgroove-tag", "Q3"))
Cgroove_sag_Q3 <- read_smooth_cell(skids_by_2annotations("celltype:Cgroove-sag", "Q3"))
monociliated <- nlapply(read.neurons.catmaid("monociliated_cell", pid = 35), function(x) smooth_neuron(x, sigma = 1000))
bounding_dots <- read.neurons.catmaid("bounding dot", pid = 35)

# 3D plotting -------------------------------------------------------------

plot_bal_cells <- function(){
  plot3d(bounding_dots, alpha = 0, lwd = 0)
  plot3d(balancer_Q1, color = Okabe_Ito[1], lwd = 2, soma = TRUE, alpha = 0.6)
  plot3d(balancer_Q2, color = Okabe_Ito[1], lwd = 2, soma = TRUE, alpha = 0.6)
  plot3d(balancer_Q3, color = Okabe_Ito[1], lwd = 2, soma = TRUE, alpha = 0.6)
  plot3d(balancer_Q4, color = Okabe_Ito[1], lwd = 2, soma = TRUE, alpha = 0.6)
  plot3d(intramulticilia, color = Okabe_Ito[2], lwd = 3, soma = TRUE, alpha = 0.6)
  plot3d(lamellate, color = Okabe_Ito[3], lwd = 3, soma = TRUE, alpha = 0.6)
#  plot3d(bridge, color = Okabe_Ito[4], lwd = 5, soma = TRUE, alpha = 1)
  plot3d(Cgroove_tag, color = Okabe_Ito[5], lwd = 3, soma = TRUE, alpha = 0.6)
  plot3d(Cgroove_sag, color = Okabe_Ito[6], lwd = 3, soma = TRUE, alpha = 1)
  plot3d(statolith, color = Okabe_Ito[8], lwd = 3, soma = TRUE, alpha = 1)
  plot3d(monociliated, soma = TRUE, color =  "grey50", alpha = 0.05, lwd = 3
)
}

{
nopen3d() 
mfrow3d(1, 3)  #defines the two scenes
#define the size of the rgl window, the view and zoom
par3d(windowRect = c(0, 0, 2400, 800))
        
plot_bal_cells()
nview3d("anterior", 
        extramat = rotationMatrix(2.54, 0.1, 0, 1)
        )
plot3d(bounding_dots, alpha = 0, lwd = 0)
par3d(zoom=0.6)

#move to next panel in rgl window
next3d(clear=F)
#plot lateral view of Sagittal plane
plot_bal_cells()
nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
plot3d(bounding_dots, alpha = 0, lwd = 0)
par3d(zoom=0.6)

#move to next panel in rgl window
next3d(clear=F)
#plot lateral view of Tentacular plane
plot_bal_cells()
nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))
plot3d(bounding_dots, alpha = 0, lwd = 0)
par3d(zoom=0.6)
}

#make a snapshot
rgl.snapshot("manuscript/pictures/balancer_statolith.png")
close3d()

# close-up view ----------------

nopen3d() 
#define the size of the rgl window
par3d(windowRect = c(0, 0, 800, 800))
nview3d("anterior", extramat = rotationMatrix(2, 1, 1, 1))

plot3d(bounding_dots, alpha = 0, lwd = 0)
par3d(zoom=0.6)       
plot3d(balancer_Q3, color = Okabe_Ito[1], lwd = 2, soma = TRUE, alpha = 0.6)
plot3d(intramulticilia_Q3, color = Okabe_Ito[2], lwd = 3, soma = TRUE, alpha = 0.6)
plot3d(lamellate_Q3, color = Okabe_Ito[3], lwd = 3, soma = TRUE, alpha = 0.6)
plot3d(Cgroove_tag_Q3, color = Okabe_Ito[5], lwd = 3, soma = TRUE, alpha = 0.6)
plot3d(Cgroove_sag_Q3, color = Okabe_Ito[6], lwd = 3, soma = TRUE, alpha = 1)
plot3d(statolith, color = Okabe_Ito[8], lwd = 3, soma = TRUE, alpha = 1)
plot3d(monociliated, soma = TRUE, color =  "grey50", alpha = 0.05, lwd = 3)
nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))

#make a snapshot
rgl.snapshot("manuscript/pictures/balancer_statolith_Q3.png")
close3d()


# crop EM pictures ---------------

crop_catmaid("nice LB", 2000, 2000, 1, 0, "manuscript/pictures/", 35, 28)

#list all files starting with crop
system("ls ./pictures/crop_*")

dcv1 <- ggdraw() + draw_image(magick::image_read("pictures/crop_6l53xv.tiff"))
