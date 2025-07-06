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
bridge <- read_smooth_neuron("celltype:bridge")
dome <- read_smooth_neuron("celltype:dome")

balancer_Q1 <- read_smooth_neuron(get_skids_with_annot(pid = 35, c("celltype:balancer", "Q1")))
balancer_Q2 <- read_smooth_neuron(get_skids_with_annot(pid = 35, c("celltype:balancer", "Q2")))
balancer_Q3 <- read_smooth_neuron(get_skids_with_annot(pid = 35, c("celltype:balancer", "Q3")))
balancer_Q4 <- read_smooth_neuron(get_skids_with_annot(pid = 35, c("celltype:balancer", "Q4")))

bridge_Q1Q2 <- read_smooth_neuron(get_skids_with_annot(pid = 35, c("celltype:bridge", "Q1Q2")))
bridge_Q3Q4 <- read_smooth_neuron(get_skids_with_annot(pid = 35, c("celltype:bridge", "Q3Q4")))

SSN_Q1Q2 <- read_smooth_neuron("SSN_Q1Q2")[[1]]
SSN_Q3Q4 <- read_smooth_neuron("SSN_Q3Q4")[[1]]
SSN_Q1Q2Q3Q4 <- read_smooth_neuron("SSN_Q1Q2Q3Q4")[[1]]

Q1 <- read_smooth_neuron("Q1")
Q2 <- read_smooth_neuron("Q2")
Q3 <- read_smooth_neuron("Q3")
Q4 <- read_smooth_neuron("Q4")

# retrieve connectors ----------------------------------------------------------

balancer_conn <- connectors(balancer)
presyn_balancer <- subset(balancer_conn, prepost == 0)
postsyn_balancer <- subset(balancer_conn, prepost == 1)

bridge_conn <- connectors(bridge)
presyn_bridge <- subset(bridge_conn, prepost == 0)
postsyn_bridge <- subset(bridge_conn, prepost == 1)



# functions ------


remove_outline <- function(){
  #as_tibble(ids3d()) #check ids in rgl plot
  rgl.pop(id = as_tibble(ids3d()) |> filter(type =="triangles") |> pull(id))
  next3d(clear=F)
  rgl.pop(id = as_tibble(ids3d()) |> filter(type =="triangles") |> pull(id))
}

add_outline <- function(){
  plot3d(
  outline,
  add = T, alpha = 0.05, col = "grey50"
  )
}

remove_text <- function(){
  rgl.pop(id = as_tibble(ids3d()) |> filter(type =="text") |> pull(id))
}

plot_balancer <- function(){
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
    bounding_dots,
    soma = F, lwd = 0, add = T,
    alpha = 0, col = "white"
  )
}

plot_lithocyte <- function(){
  plot3d(
    lithocyte,
    soma = T, lwd = 1.5, add = T,
    alpha = 0.9, col = Okabe_Ito[8]
  )
}

plot_bridge <- function(){
  plot3d(
    bridge_Q1Q2,
    soma = T, lwd = 3, add = T,
    alpha = 0.8, col = Okabe_Ito[2]
  )
  plot3d(
    bridge_Q3Q4,
    soma = T, lwd = 3, add = T,
    alpha = 0.8, col = Okabe_Ito[7]
  )
}

# plot cells -------------------------------------------------------------------

nopen3d() # opens a pannable 3d window
mfrow3d(1, 2)  #defines the two scenes
par3d(windowRect = c(0, 0, 1600, 800)) #to define the size of the rgl window
par3d(zoom=0.9)
#nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
#nview3d("ventral", extramat=rotationMatrix(1.2, 0, 0, 1))
nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))


plot_balancer()
add_outline()
plot3d(
    lithocyte,
    soma = T, lwd = 1, add = T,
    alpha = 0, col = "white"
  )
par3d(zoom=0.6)
texts3d(
  14000, 32000, 1000, text = "balancer", col='black', cex = 2
  )


#go to next scene
next3d(clear=FALSE)
nview3d("posterior", 
        extramat = rotationMatrix(2.54, 0.1, 0, 1)
        )

plot_balancer()
add_outline()
#add
plot3d(
    lithocyte,
    soma = T, lwd = 1, add = T,
    alpha = 0, col = "white"
  )

par3d(zoom=0.7)
texts3d(11000, 52000, 23000, text = "Q1", cex = 2)
texts3d(22000, 30000, 23000, text = "Q2", cex = 2)
texts3d(30000, 68000, 23000, text = "Q3", cex = 2)
texts3d(42000, 43000, 23000, text = "Q4", cex = 2)

for(i in 101:120){
  rgl.snapshot(paste("videoframes/Video1_", i, ".png", sep = ""))
}
#remove lithocytes
id_litho = tail(as_tibble(ids3d()) |> filter(type =="spheres") |> pull(id), 8) 
rgl.pop(id = id_litho)

next3d(clear=FALSE)
#remove text
remove_text()
#remove lithocytes
id_litho = tail(as_tibble(ids3d()) |> filter(type =="spheres") |> pull(id), 8) 
rgl.pop(id = id_litho)

#plot lithocytes ---------

plot_lithocyte()
texts3d(14000, 32000, 1000, text = "lithocytes", cex = 2)
next3d(clear=FALSE)
plot_lithocyte()

for(i in 121:140){
  rgl.snapshot(paste("videoframes/Video1_", i, ".png", sep = ""))
}

next3d(clear=FALSE)
#remove text
remove_text()

#plot bridge -----------

plot_bridge()

texts3d(
  15000, 32000, 1000, text = "bridge", col='black', cex = 2
)
next3d(clear=FALSE)
plot_bridge()

for(i in 141:160){
  rgl.snapshot(paste("videoframes/Video1_", i, ".png", sep = ""))
}

next3d(clear=FALSE)
#remove text
remove_text()

#plot dome -------------------


plot3d(
  dome, soma = T, lwd = 1, add = T,
  alpha = 1, col = "grey90"
)
texts3d(
  14000, 32000, 1000, text = "dome cells", col='black', cex = 2
)

next3d(clear=FALSE)
plot3d(
  dome, soma = T, lwd = 1, add = T,
  alpha = 1, col = "grey90"
)

for(i in 161:180){
  rgl.snapshot(paste("videoframes/Video1_", i, ".png", sep = ""))
}

next3d(clear=FALSE)
#remove text
remove_text()


#plot nerve net ---------------

remove_outline()

plot_multinucleated_cell(
  SSN_Q1Q2, lwd = 2, alpha = 1, col = Okabe_Ito[6]
  )
plot_multinucleated_cell(
  SSN_Q3Q4, lwd = 2, alpha = 1, col = Okabe_Ito[7]
  )

add_outline()

next3d(clear=FALSE)

plot_multinucleated_cell(
  SSN_Q1Q2, lwd = 2, alpha = 1, col = Okabe_Ito[6]
)
plot_multinucleated_cell(
  SSN_Q3Q4, lwd = 2, alpha = 1, col = Okabe_Ito[7]
)

add_outline()
texts3d(
  15000, 32000, 1000, text = "nerve net Q1Q2, Q3Q4", col='black', cex = 2
)
for(i in 181:200){
  rgl.snapshot(paste("videoframes/Video1_", i, ".png", sep = ""))
}

#plot large nerve net -----------
remove_outline()
next3d(clear=FALSE)
plot_multinucleated_cell(SSN_Q1Q2Q3Q4,
                         lwd = 3, alpha = 0.6, col = Okabe_Ito[5])
remove_text()
texts3d(
  15000, 32000, 1000, text = "nerve net Q1-4", col='black', cex = 2
)
add_outline()
next3d(clear=FALSE)
plot_multinucleated_cell(SSN_Q1Q2Q3Q4,
                         lwd = 4, alpha = 0.6, col = Okabe_Ito[5])
add_outline()

for(i in 201:220){
  rgl.snapshot(paste("videoframes/Video1_", i, ".png", sep = ""))
}

# input from SSN to balancer -----------------
next3d(clear=FALSE)

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
remove_text()

texts3d(
  15000, 32000, 1000, text = "synapses to balancer", col='black', cex = 2
)

for(i in 221:240){
  rgl.snapshot(paste("videoframes/Video1_", i, ".png", sep = ""))
}

remove_text()

#plot bridge synapses --------------

# input from SSN to balancer
plot3d(
  postsyn_bridge$x, 
  postsyn_bridge$y, 
  postsyn_bridge$z, 
  size = 0.7, alpha = 1, col = "#A52A2A", 
  add = TRUE,
  point_antialias = TRUE,
  type = "s"
)

plot3d(
  presyn_bridge$x, 
  presyn_bridge$y, 
  presyn_bridge$z, 
  size = 0.8, alpha = 1, col = "black", 
  add = TRUE,
  point_antialias = TRUE,
  type = "s"
)

texts3d(
  15000, 32000, -1000, text = "synapses to bridge", col='#A52A2A', cex = 2
)
texts3d(
  15000, 32000, 1000, text = "synapses from bridge", col='black', cex = 2
)

next3d(clear=F)

plot3d(
  postsyn_bridge$x, 
  postsyn_bridge$y, 
  postsyn_bridge$z, 
  size = 0.7, alpha = 1, col = "#A52A2A", 
  add = TRUE,
  point_antialias = TRUE,
  type = "s"
)

plot3d(
  presyn_bridge$x, 
  presyn_bridge$y, 
  presyn_bridge$z, 
  size = 0.8, alpha = 1, col = "black", 
  add = TRUE,
  point_antialias = TRUE,
  type = "s"
)

for(i in 241:260){
  rgl.snapshot(paste("videoframes/Video1_", i, ".png", sep = ""))
}

remove_text()

# get rotation matrix
next3d(clear=F)
remove_text()

nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))

um1 <- par3d()$userMatrix
next3d(clear=F)
#nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
nview3d("posterior", 
        extramat = rotationMatrix(2.54, 0.1, 0, 1)
)

um2 <- par3d()$userMatrix
next3d(clear=F)

rotation=300

for (l in 1:90){
  #rotate in a loop (with l e.g. 1:90 for a 180 turn)
  nview3d(userMatrix = um1 %*%rotationMatrix(pi*l/90, 0, 0, 1)
          %*%rotationMatrix(-0.3*l/90, 0, 1, 0)
          %*%rotationMatrix(0.2*l/90, 1, 0, 0)
          )
  next3d(clear=F)
  nview3d(userMatrix = um2 %*%rotationMatrix(pi*l/90, 0, 0, 1)
          %*%rotationMatrix(-0.3*l/90, 0, 1, 0)
          %*%rotationMatrix(0.2*l/90, 1, 0, 0)
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


#get coordinates of the actual state of the rotation --------------
um3 <- par3d()$userMatrix
next3d(clear=F)
um4 <- par3d()$userMatrix

# remove outline ----------

remove_outline()

plot3d(Q1, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.1, col = Okabe_Ito[1])
plot3d(Q2, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.1, col = Okabe_Ito[2])
plot3d(Q3, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.1, col = Okabe_Ito[6])
plot3d(Q4, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.1, col = Okabe_Ito[7])
add_outline()

next3d(clear=F)
plot3d(Q1, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.1, col = Okabe_Ito[1])
plot3d(Q2, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.1, col = Okabe_Ito[2])
plot3d(Q3, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.1, col = Okabe_Ito[6])
plot3d(Q4, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.1, col = Okabe_Ito[7])
add_outline()

next3d(clear=F)
rotation=400

for (l in 1:90){
  #rotate in a loop (with l e.g. 1:90 for a 180 turn)
  nview3d(userMatrix = um3 %*%rotationMatrix(pi*l/90, 0, 0, 1)
          %*%rotationMatrix(-0.3*l/90, 0, 1, 0)
          %*%rotationMatrix(0.2*l/90, 1, 0, 0)
  )
  next3d(clear=F)
  nview3d(userMatrix = um4 %*%rotationMatrix(pi*l/90, 0, 0, 1)
          %*%rotationMatrix(-0.3*l/90, 0, 1, 0)
          %*%rotationMatrix(0.2*l/90, 1, 0, 0)
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


#read png files and write video --------------------------------------------------
av::av_encode_video(
  paste('videoframes/', list.files("videoframes/", '*.png'), 
        sep = ""),
  framerate = 10,
  output = 'manuscript/videos/Video1.mp4'
  )

unlink("videoframes", recursive = T)


