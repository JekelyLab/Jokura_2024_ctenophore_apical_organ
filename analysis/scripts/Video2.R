#Video2 of the Jokura et al ctenophore AO paper
#Gaspar Jekely

# source packages and functions ------------------------------------------------
source("analysis/scripts/packages_and_functions.R")

# load cell types to plot individually -------------------------------------

balancer <- read_smooth_neuron("celltype:balancer")
lithocyte <- read_smooth_neuron("statolith")

balancer_Q1 <- read_smooth_neuron(get_skids_with_annot(pid = 35, c("celltype:balancer", "Q1")))
balancer_Q2 <- read_smooth_neuron(get_skids_with_annot(pid = 35, c("celltype:balancer", "Q2")))
balancer_Q3 <- read_smooth_neuron(get_skids_with_annot(pid = 35, c("celltype:balancer", "Q3")))
balancer_Q4 <- read_smooth_neuron(get_skids_with_annot(pid = 35, c("celltype:balancer", "Q4")))

# read basal bodies -------------------------------------------------------

#get skids
balancer_skids <- catmaid_skids("celltype:balancer", pid=35)

#read all nodes
nodes <- lapply(
  balancer_skids, 
  function(x) catmaid_get_treenode_table(x, pid = 35)
)

#get tags
tags_with_id <- lapply(nodes, function(x) 
  catmaid_get_labels(treenodes = x$id, pid = 35))

#select basal body tags
bb_nodes_id <- lapply(tags_with_id, function(x) x %>%
                        as_tibble() %>%
                        filter("basal body" == label) %>%
                        select(id) %>%
                        pull()
)

#ids of nodes with bb tag
bb_nodes_id <- unlist(bb_nodes_id)
bb_nodes_id

#coordinates
bb_xyz <- lapply(nodes, function(k) k %>%
                   filter(id %in% bb_nodes_id) %>%
                   select(x,y,z)
)
bb_xyz <- bind_rows(bb_xyz)

# open image with three panels -------------------------------------------------------------

nopen3d()
mfrow3d(1, 3) 
# define the size of the rgl window, the view and zoom
#par3d(windowRect = c(0, 0, 1200, 350))
par3d(windowRect = c(0, 0, 2400, 700))

# function to plot balancer cells  ---------------

plot_cells <- function(){
  plot3d(
    balancer_Q1,
    soma = T, lwd = 1.5, add = T,
    alpha = 0.5, col = Okabe_Ito[1]
  )
  plot3d(
    balancer_Q2,
    soma = T, lwd = 1.5, add = T,
    alpha = 0.5, col = Okabe_Ito[2]
  )
  plot3d(
    balancer_Q3,
    soma = T, lwd = 1.5, add = T,
    alpha = 0.5, col = Okabe_Ito[6]
  )
  plot3d(
    balancer_Q4,
    soma = T, lwd = 1.5, add = T,
    alpha = 0.5, col = Okabe_Ito[7]
  )
  
  plot3d(
    outline,
    add = T, alpha = 0.05, col = "grey80"
  )
  plot3d(
    lithocyte,
    soma = T, lwd = 0, add = T,
    alpha = 0.6, col = "grey50"
  )
  # input from SSN to balancer
  plot3d(
    bb_xyz$x, 
    bb_xyz$y, 
    bb_xyz$z, 
    size = 0.5, alpha = 1, col = "red", 
    add = TRUE,
    point_antialias = TRUE,
    type = "s"
  )
}

# plot cells in three different views ------------

plot_cells()
#aboral view
aboral()
par3d(zoom = 0.61)
# move to next panel in rgl window
next3d(clear = F)

# plot lateral view of Sagittal plane
plot_cells()
sagittal()
par3d(zoom = 0.61)
# move to next panel in rgl window
next3d(clear = F)

# plot lateral view of Tentacular plane
plot_cells()
tentacular()
par3d(zoom = 0.61)

# make a snapshot to the working directory
rgl.snapshot("manuscript/pictures/balancer_with_bb.png")

# plot only one panel -----------------

nopen3d()
# define the size of the rgl window, the view and zoom
par3d(windowRect = c(0, 0, 800, 700))
# plot lateral view of Tentacular plane
plot_cells()
tentacular()
par3d(zoom = 0.61)

#create temp dir to store video frames -----------

mainDir = getwd()
dir.create(file.path(mainDir, "videoframes"), showWarnings = FALSE)


# rotation -----------------


# get rotation matrix
um1 <- par3d()$userMatrix
rotation=300

for (l in 1:90){
  #rotate in a loop (with l e.g. 1:90 for a 180 turn)
  nview3d(userMatrix = um1 %*%rotationMatrix(pi*l/90, 0, 0, 1)
          %*%rotationMatrix(-0.3*l/90, 0, 1, 0)
          %*%rotationMatrix(0.2*l/90, 1, 0, 0)
  )
  print (l)
  #save a snapshot
  filename <- paste("./videoframes/Video2_",
                    rotation,
                    formatC(l, digits = 2, flag = "0"),
                    ".png", sep = "")
  rgl.snapshot(filename)
  rotation = rotation + 1 
}


close3d()


#read png files and write video ----------------


av::av_encode_video(
  paste('videoframes/', list.files("videoframes/", '*.png'), 
        sep = ""),
  framerate = 10,
  output = 'manuscript/videos/Video2.mp4'
)

# delete temp folder -----------------
unlink("videoframes", recursive = T)



