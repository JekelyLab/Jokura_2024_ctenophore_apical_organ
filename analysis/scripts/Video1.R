#Video1 of the Jokura et al ctenophore AO paper
#Sanja Jasek, Kei Jokura, Gaspar Jekely

#load packages and functions ------------
source("analysis/scripts/packages_and_functions.R")

#create temp dir to store video frames -----------
mainDir = getwd()
dir.create(file.path(mainDir, "videoframes"), showWarnings = FALSE)

# load cell types from CATMAID ---------

#with_soma <- read_smooth_neuron("with_soma")
bounding_dots <- read.neurons.catmaid("bounding dot", pid = 35)
balancer <- read_smooth_neuron("celltype:balancer")
lithocyte <- read_smooth_neuron("celltype:lithocyte")

balancer_Q1 <- read_smooth_neuron(get_skids_with_annot(pid = 35, c("celltype:balancer", "Q1")))
balancer_Q2 <- read_smooth_neuron(get_skids_with_annot(pid = 35, c("celltype:balancer", "Q2")))
balancer_Q3 <- read_smooth_neuron(get_skids_with_annot(pid = 35, c("celltype:balancer", "Q3")))
balancer_Q4 <- read_smooth_neuron(get_skids_with_annot(pid = 35, c("celltype:balancer", "Q4")))

bridge_Q1Q2 <- read_smooth_neuron(get_skids_with_annot(pid = 35, c("celltype:bridge", "Q1Q2")))
bridge_Q3Q4 <- read_smooth_neuron(get_skids_with_annot(pid = 35, c("celltype:bridge", "Q3Q4")))

SSN_Q1Q2 <- read_smooth_neuron("SSN_Q1Q2")[[1]]
SSN_Q3Q4 <- read_smooth_neuron("SSN_Q3Q4")[[1]]
SSN_Q1Q2Q3Q4 <- read_smooth_neuron("SSN_Q1Q2Q3Q4")[[1]]


# retrieve connectors ----------------------------------------------------------

balancer_conn <- connectors(balancer)
presyn_balancer <- subset(balancer_conn, prepost == 0)
postsyn_balancer <- subset(balancer_conn, prepost == 1)

# plot cells ----------------

nopen3d() # opens a pannable 3d window
mfrow3d(1, 2)  #defines the two scenes
par3d(windowRect = c(0, 0, 1600, 800)) #to define the size of the rgl window
par3d(zoom=0.9)
#nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
#nview3d("ventral", extramat=rotationMatrix(1.2, 0, 0, 1))
nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))



plot3d(
  balancer_Q1,
  soma = T, lwd = 1.5, add = T,
  alpha = 0.6, col = Okabe_Ito[1]
)
plot3d(
  balancer_Q2,
  soma = T, lwd = 1.5, add = T,
  alpha = 0.6, col = Okabe_Ito[2]
)
plot3d(
  balancer_Q3,
  soma = T, lwd = 1.5, add = T,
  alpha = 0.6, col = Okabe_Ito[6]
)
plot3d(
  balancer_Q4,
  soma = T, lwd = 1.5, add = T,
  alpha = 0.6, col = Okabe_Ito[7]
)
plot3d(
  lithocyte,
  soma = T, lwd = 1.5, add = T,
  alpha = 0.9, col = Okabe_Ito[8]
)
#plot3d(
#  with_soma,
#  soma = T, lwd = 1, add = T,
#  alpha = 0.01, col = Okabe_Ito[8]
#)

plot3d(
  outline,
  add = T, alpha = 0.05, col = "grey50"
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
  balancer_Q1,
  soma = T, lwd = 1.5, add = T,
  alpha = 0.6, col = Okabe_Ito[1]
)
plot3d(
  balancer_Q2,
  soma = T, lwd = 1.5, add = T,
  alpha = 0.6, col = Okabe_Ito[2]
)
plot3d(
  balancer_Q3,
  soma = T, lwd = 1.5, add = T,
  alpha = 0.6, col = Okabe_Ito[6]
)
plot3d(
  balancer_Q4,
  soma = T, lwd = 1.5, add = T,
  alpha = 0.6, col = Okabe_Ito[7]
)
plot3d(
  lithocyte,
  soma = T, lwd = 1.5, add = T,
  alpha = 0.9, col = Okabe_Ito[8]
)
#plot3d(
#  with_soma,
#  soma = T, lwd = 1, add = T,
#  alpha = 0.01, col = Okabe_Ito[8]
#)

plot3d(
  outline,
  add = T, alpha = 0.05, col = "grey50"
)

par3d(zoom=0.7)

for(i in 101:120){
  rgl.snapshot(paste("videoframes/Video1_", i, ".png", sep = ""))
  }

next3d(clear=FALSE)
#remove text
rgl.pop(id = as_tibble(ids3d()) |> filter(type =="text") |> pull(id))


plot3d(
  bridge_Q1Q2,
  soma = T, lwd = 1.5, add = T,
  alpha = 0.8, col = Okabe_Ito[2]
)
plot3d(
  bridge_Q3Q4,
  soma = T, lwd = 1.5, add = T,
  alpha = 0.8, col = Okabe_Ito[7]
)
texts3d(
  15000, 32000, 1000, text = "bridge", col='black', cex = 2
)
next3d(clear=FALSE)
plot3d(
  bridge_Q1Q2,
  soma = T, lwd = 1.5, add = T,
  alpha = 0.8, col = Okabe_Ito[2]
)
plot3d(
  bridge_Q3Q4,
  soma = T, lwd = 1.5, add = T,
  alpha = 0.8, col = Okabe_Ito[7]
)
for(i in 121:140){
  rgl.snapshot(paste("videoframes/Video1_", i, ".png", sep = ""))
}

next3d(clear=FALSE)
#remove text
rgl.pop(id = as_tibble(ids3d()) |> filter(type =="text") |> pull(id))

#plot nerve net ---------------
plot_multinucleated_cell(SSN_Q1Q2,
                         lwd = 2, alpha = 1, col = Okabe_Ito[6])
plot_multinucleated_cell(SSN_Q3Q4,
                         lwd = 2, alpha = 1, col = Okabe_Ito[7])
plot_multinucleated_cell(SSN_Q1Q2Q3Q4,
                         lwd = 2, alpha = 1, col = Okabe_Ito[5])
texts3d(
  15000, 32000, 1000, text = "nerve net", col='black', cex = 2
)
# input from SSN to balancer
plot3d(
  postsyn_balancer$x, 
  postsyn_balancer$y, 
  postsyn_balancer$z, 
  size = 0.6, alpha = 1, col = "red", 
  add = TRUE,
  point_antialias = TRUE,
  type = "s"
)


next3d(clear=FALSE)
plot_multinucleated_cell(SSN_Q1Q2,
                         lwd = 2, alpha = 1, col = Okabe_Ito[6])
plot_multinucleated_cell(SSN_Q3Q4,
                         lwd = 2, alpha = 1, col = Okabe_Ito[7])
plot_multinucleated_cell(SSN_Q1Q2Q3Q4,
                         lwd = 2, alpha = 1, col = Okabe_Ito[5])

# input from SSN to balancer
plot3d(
  postsyn_balancer$x, 
  postsyn_balancer$y, 
  postsyn_balancer$z, 
  size = 0.6, alpha = 1, col = "red", 
  add = TRUE,
  point_antialias = TRUE,
  type = "s"
)


for(i in 141:160){
  rgl.snapshot(paste("videoframes/Video1_", i, ".png", sep = ""))
}

next3d(clear=FALSE)
rgl.pop(id = as_tibble(ids3d()) |> filter(type =="text") |> pull(id))


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


