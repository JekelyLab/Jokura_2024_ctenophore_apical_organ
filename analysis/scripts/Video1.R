#Video1 of the Jokura et al ctenophore AO paper

source("analysis/scripts/packages_and_functions.R")

#create temp dir to store video frames
mainDir = getwd()
dir.create(file.path(mainDir, "videoframes"), showWarnings = FALSE)


# load cell types from catmaid

balancer <- nlapply(
  read.neurons.catmaid(
    "celltype:balancer", pid = 35
    ),
  function(x) smooth_neuron(x, sigma = 1000)
)

intramulticilia_Q1 <- read_smooth_cell(
  get_skids_with_annot(pid = 35, c(
    "celltype:intra-multi-ciliated", "Q1"))
  )

LB <- nlapply(
  read.neurons.catmaid(
    "celltype:lamellate", pid = 35
    ),
  function(x) smooth_neuron(x, sigma = 1000)
)

syn_neuron <- nlapply(
  read.neurons.catmaid(
    "celltype:neuron", pid = 35
    ),
  function(x) smooth_neuron(x, sigma = 1000)
)

monociliated <- nlapply(
  read.neurons.catmaid(
    "monociliated_cell", pid = 35
  ),
  function(x) smooth_neuron(x, sigma = 1000)
)

bridge <- nlapply(
  read.neurons.catmaid(
    "celltype:bridge", pid = 35
  ),
  function(x) smooth_neuron(x, sigma = 1000)
)

# cell types
AO_celltypes <- list(
  balancer, LB, syn_neuron, bridge, monociliated
  )
AO_celltype_names <- list(
  "balancer", "LB", "syn_neuron", "bridge", "monociliated"
  )

bounding_dots <- read.neurons.catmaid("bounding dot", pid = 35)

# plot cells ----------------

nopen3d() # opens a pannable 3d window
mfrow3d(1, 2)  #defines the two scenes
par3d(windowRect = c(0, 0, 1600, 800)) #to define the size of the rgl window
par3d(zoom=0.9)
#nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
#nview3d("ventral", extramat=rotationMatrix(1.2, 0, 0, 1))
nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))


plot3d(
  balancer, soma = TRUE, color = Okabe_Ito[1], 
  alpha = 0.6, lwd = 3
  )
plot3d(bounding_dots, alpha = 0, lwd = 0)
par3d(zoom=0.6)
texts3d(
  14000, 32000, 1000, text = "balancer", col='black', cex = 2
  )

#go to next scene
next3d(clear=FALSE)
nview3d("anterior", 
        extramat = rotationMatrix(2.54, 0.1, 0, 1)
        )
plot3d(
  balancer, soma = TRUE, color = Okabe_Ito[1], 
  alpha = 0.6, lwd = 3
)
plot3d(bounding_dots, alpha = 0, lwd = 0)
par3d(zoom=0.6)

for(i in 101:120){
  rgl.snapshot(paste("videoframes/Video1_", i, ".png", sep = ""))
  }

next3d(clear=FALSE)
#remove text
rgl.pop(id = as_tibble(ids3d()) |> filter(type =="text") |> pull(id))

plot3d(
  intramulticilia, soma = TRUE, color = Okabe_Ito[2], 
  alpha = 0.6, lwd = 2
)
texts3d(
  15000, 32000, 1000, text = "intra-multiciliated", col='black', cex = 2
)
next3d(clear=FALSE)
plot3d(
  intramulticilia, soma = TRUE, color = Okabe_Ito[2], 
  alpha = 0.6, lwd = 2
)

for(i in 121:140){
  rgl.snapshot(paste("videoframes/Video1_", i, ".png", sep = ""))
}

next3d(clear=FALSE)
#remove text
rgl.pop(id = as_tibble(ids3d()) |> filter(type =="text") |> pull(id))

plot3d(
  LB, soma = TRUE, color = Okabe_Ito[3], 
  alpha = 0.6, lwd = 2
)
texts3d(
  15000, 32000, 1000, text = "lamellate body", col='black', cex = 2
)

next3d(clear=FALSE)
plot3d(
  LB, soma = TRUE, color = Okabe_Ito[3], 
  alpha = 0.6, lwd = 2
)

for(i in 141:160){
  rgl.snapshot(paste("videoframes/Video1_", i, ".png", sep = ""))
}

next3d(clear=FALSE)
rgl.pop(id = as_tibble(ids3d()) |> filter(type =="text") |> pull(id))

plot3d(
  syn_neuron, soma = TRUE, color = Okabe_Ito[4], 
  alpha = 0.6, lwd = 4
)
texts3d(
  15000, 32000, 1000, text = "nerve net", col='black', cex = 2
)

next3d(clear=FALSE)
plot3d(
  syn_neuron, soma = TRUE, color = Okabe_Ito[4], 
  alpha = 0.6, lwd = 4
)


for(i in 161:180){
  rgl.snapshot(paste("videoframes/Video1_", i, ".png", sep = ""))
}

next3d(clear=FALSE)
rgl.pop(id = as_tibble(ids3d()) |> filter(type =="text") |> pull(id))


plot3d(
  bridge, soma = TRUE, color = Okabe_Ito[8], 
  alpha = 0.6, lwd = 3
)
texts3d(
  15000, 32000, 1000, text = "bridge", col='black', cex = 2
)

next3d(clear=FALSE)
plot3d(
  bridge, soma = TRUE, color = Okabe_Ito[8], 
  alpha = 0.6, lwd = 3
)

for(i in 181:200){
  rgl.snapshot(paste("videoframes/Video1_", i, ".png", sep = ""))
}
next3d(clear=FALSE)
rgl.pop(id = as_tibble(ids3d()) |> filter(type =="text") |> pull(id))

plot3d(
  monociliated, soma = TRUE, color = "grey50", 
  alpha = 0.1, lwd = 1
)

next3d(clear=FALSE)
plot3d(
  monociliated, soma = TRUE, color = "grey50", 
  alpha = 0.1, lwd = 1
)

for(i in 201:220){
  rgl.snapshot(paste("videoframes/Video1_", i, ".png", sep = ""))
}


# get rotation matrix
next3d(clear=F)

um1 <- par3d()$userMatrix
next3d(clear=F)
#nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
um2 <- par3d()$userMatrix
next3d(clear=F)

rotation=300

for (l in 1:90){
  #rotate in a loop (with l e.g. 1:90 for a 180 turn)
  nview3d(userMatrix = um1 %*%rotationMatrix(pi*l/90, 0, 0, 1)
          %*%rotationMatrix(0, 0, 1, 0)
          %*%rotationMatrix(0.2, 1, 0, 0)
          )
  next3d(clear=F)
  nview3d(userMatrix = um2 %*%rotationMatrix(pi*l/90, 0, 0, 1)
          %*%rotationMatrix(0, 0, 1, 0)
          %*%rotationMatrix(0, 1, 0, 0)
  )
  next3d(clear=F)
  print (l)
  #save a snapshot
  filename <- paste("./videoframes/Video1_",
                    rotation,
                    formatC(l, digits = 2, flag = "0"),
                    ".png", sep = "")
  rgl.snapshot(filename)
  rotation = rotation + 1 
}


close3d()


#read png files and write video
av::av_encode_video(
  paste('videoframes/', list.files("videoframes/", '*.png'), 
        sep = ""),
  framerate = 10,
  output = 'manuscript/videos/Video1.mp4'
  )

unlink("videoframes", recursive = T)


