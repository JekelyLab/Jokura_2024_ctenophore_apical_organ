#Video1 of the Jokura et al ctenophore AO paper
#Sanja Jasek, Kei Jokura, Gaspar Jekely

#load packages and functions ------------
source("analysis/scripts/packages_and_functions.R")

#create temp dir to store video frames -----------
mainDir = getwd()
dir.create(file.path(mainDir, "videoframes"), showWarnings = FALSE)

# load cell types from CATMAID ---------

#with_soma <- read_smooth_neuron("with_soma")
bounding_dots <- read.neurons.catmaid("bounding dot", pid = pid)
balancer <- read_smooth_neuron("celltype:balancer")
statolith <- read_smooth_neuron("statolith") # these are only the lithocytes on the balancer
bridge <- read_smooth_neuron("celltype:bridge")
dome <- read_smooth_neuron("celltype:dome")

balancer_Q1 <- read_smooth_neuron(get_skids_with_annot(pid = pid, c("celltype:balancer", "Q1")))
balancer_Q2 <- read_smooth_neuron(get_skids_with_annot(pid = pid, c("celltype:balancer", "Q2")))
balancer_Q3 <- read_smooth_neuron(get_skids_with_annot(pid = pid, c("celltype:balancer", "Q3")))
balancer_Q4 <- read_smooth_neuron(get_skids_with_annot(pid = pid, c("celltype:balancer", "Q4")))

bridge_Q1Q2 <- read_smooth_neuron(get_skids_with_annot(pid = pid, c("celltype:bridge", "Q1Q2")))
bridge_Q3Q4 <- read_smooth_neuron(get_skids_with_annot(pid = pid, c("celltype:bridge", "Q3Q4")))

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

draw_everything <- function() {
  add_outline()
  
  plot3d(
    bounding_dots,
    soma = F, lwd = 0, add = T,
    alpha = 0, col = "white"
  )
  
  plot3d(
    statolith,
    soma = T, lwd = 1.5, add = T,
    alpha = 1, col = "grey"
  )
  
  plot3d(
    balancer,
    soma = T, lwd = 1.5, add = T,
    alpha = 0.2, col = Okabe_Ito[1]
  )
  
  plot3d(
    bridge,
    soma = T, lwd = 3, add = T,
    alpha = 0.2, col = Okabe_Ito[2]
  )
  
  plot_multinucleated_cell(
    SSN_Q1Q2, lwd = 2, alpha = 1, col = Okabe_Ito[3]
  )
  plot_multinucleated_cell(
    SSN_Q3Q4, lwd = 2, alpha = 1, col = Okabe_Ito[5]
  )
  
  plot_multinucleated_cell(SSN_Q1Q2Q3Q4,
                           lwd = 3, alpha = 1, col = Okabe_Ito[7])
  
  plot3d(
    postsyn_balancer$x, 
    postsyn_balancer$y, 
    postsyn_balancer$z, 
    size = 0.6, alpha = 1, col = "red", 
    add = TRUE,
    point_antialias = TRUE,
    type = "s"
  )
  
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

  
}


# plot cells -------------------------------------------------------------------

nopen3d() # opens a pannable 3d window
par3d(windowRect = c(0, 0, 800, 800)) #to define the size of the rgl window
par3d(zoom=0.7)

# get views we want for later
nview3d("posterior", extramat = rotationMatrix(2.54, 0.1, 0, 1))
um1 <- par3d("userMatrix")

nview3d("right", extramat = rotationMatrix(-1.7, 190, -120, -140))
um2 <- par3d("userMatrix")


# start actual video

nview3d("posterior", extramat = rotationMatrix(2.54, 0.1, 0, 1))
#nview3d("anterior", extramat = rotationMatrix(2.54, 0.1, 0, 1))

add_outline()

plot3d(
  bounding_dots,
  soma = F, lwd = 0, add = T,
  alpha = 0, col = "white"
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


#texts3d(
#  14000, 32000, 1000, text = "balancer", col='black', cex = 2
#  )

#texts3d(
#  7000, 55000, 30000, text = "balancer", col='black', cex = 2
#)

texts3d(
  5000, 58000, 1000, text = "balancer", col='black', cex = 2
)

texts3d(11000, 52000, 23000, text = "Q1", cex = 2)
texts3d(22000, 30000, 23000, text = "Q2", cex = 2)
texts3d(31000, 67000, 23000, text = "Q3", cex = 2)
texts3d(44000, 45000, 23000, text = "Q4", cex = 2)

# plot invisible lithocyte so that the view doesn't move from this frame to the next
plot3d(
  statolith,
  soma = T, lwd = 1.5, add = T,
  alpha = 0, col = "grey"
)

for(i in 101:120){
  rgl.snapshot(paste("videoframes/Video1_", i, ".png", sep = ""))
}

# remove all to recolor
clear3d()

add_outline()


plot3d(
  bounding_dots,
  soma = F, lwd = 0, add = T,
  alpha = 0, col = "white"
)


plot3d(
  balancer,
  soma = T, lwd = 1.5, add = T,
  alpha = 0.2, col = Okabe_Ito[1]
)

plot3d(
  statolith,
  soma = T, lwd = 1.5, add = T,
  alpha = 0, col = "grey"
)


for(i in 121:125){
  rgl.snapshot(paste("videoframes/Video1_", i, ".png", sep = ""))
}

n_stat= length(statolith)
id_stat = tail(as_tibble(ids3d()) |> filter(type =="spheres" | type == "linestrip") |> pull(id), n_stat*2) 
rgl.pop(id = id_stat)

plot3d(
  statolith,
  soma = T, lwd = 1.5, add = T,
  alpha = 0.9, col = "grey"
)
texts3d(5000, 58000, 1000, text = "statoliths", cex = 2)

for(i in 126:140){
  rgl.snapshot(paste("videoframes/Video1_", i, ".png", sep = ""))
}

#remove text
remove_text()

#plot bridge -----------

texts3d(
  5000, 58000, 1000, text = "bridge", col='black', cex = 2
)

plot3d(
  bridge_Q1Q2,
  soma = T, lwd = 3, add = T,
  alpha = 0.8, col = Okabe_Ito[2]
)
plot3d(
  bridge_Q3Q4,
  soma = T, lwd = 3, add = T,
  alpha = 0.8, col = Okabe_Ito[3]
)

for(i in 141:160){
  rgl.snapshot(paste("videoframes/Video1_", i, ".png", sep = ""))
}


n_bridge= length(bridge)
id_bridge = tail(as_tibble(ids3d()) |> filter(type =="spheres" | type == "linestrip") |> pull(id), n_bridge*2) 
rgl.pop(id = id_bridge)

plot3d(
  bridge,
  soma = T, lwd = 3, add = T,
  alpha = 0.2, col = Okabe_Ito[2]
)

#remove text
remove_text()


#plot nerve net ---------------


plot_multinucleated_cell(
  SSN_Q1Q2, lwd = 2, alpha = 1, col = Okabe_Ito[3]
  )
plot_multinucleated_cell(
  SSN_Q3Q4, lwd = 2, alpha = 1, col = Okabe_Ito[5]
  )


#texts3d(
#  15000, 32000, 1000, text = "nerve net Q1Q2, Q3Q4", col='black', cex = 2
#)
texts3d(
  5000, 58000, 1000, text = "nerve net Q1Q2, Q3Q4", col='black', cex = 2
)

for(i in 161:180){
  rgl.snapshot(paste("videoframes/Video1_", i, ".png", sep = ""))
}

#plot large nerve net -----------

plot_multinucleated_cell(SSN_Q1Q2Q3Q4,
                         lwd = 3, alpha = 1, col = Okabe_Ito[7])
remove_text()
#texts3d(
#  15000, 32000, 1000, text = "nerve net Q1-4", col='black', cex = 2
#)
texts3d(
  5000, 58000, 1000, text = "nerve net Q1-4", col='black', cex = 2
)

for(i in 181:200){
  rgl.snapshot(paste("videoframes/Video1_", i, ".png", sep = ""))
}

remove_text()

# input from SSN to balancer -----------------

plot3d(
  postsyn_balancer$x, 
  postsyn_balancer$y, 
  postsyn_balancer$z, 
  size = 0.6, alpha = 1, col = "red", 
  add = TRUE,
  point_antialias = TRUE,
  type = "s"
)

#texts3d(
#  15000, 32000, 1000, text = "synapses to balancer", col='black', cex = 2
#)
texts3d(
  5000, 58000, 1000, text = "synapses to balancer", col='black', cex = 2
)

for(i in 201:220){
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
  4500, 60000, 1000, text = "synapses to bridge", col='#A52A2A', cex = 2
)
texts3d(
  5000, 58000, 1000, text = "synapses from bridge", col='black', cex = 2
)


for(i in 221:240){
  rgl.snapshot(paste("videoframes/Video1_", i, ".png", sep = ""))
}

remove_text()

for(i in 241:230){
  rgl.snapshot(paste("videoframes/Video1_", i, ".png", sep = ""))
}



# capture full current view
p1 <- par3d()

# insert your starting and ending matrices
p1$userMatrix <- um1
p2 <- p1
p2$userMatrix <- um2


ip <- par3dinterp(times = c(0, 1),
                  userMatrix = list(um1, um2))


rotation=300
# redraw everything in every frame, because it flickers less
for (l in 1:90) {
  nopen3d() 
  par3d(windowRect = c(0, 0, 800, 800)) 
  par3d(zoom=0.7)
  draw_everything()
  
  t <- l / 90
  par3d(userMatrix = ip(t)$userMatrix)

  material3d(depth_mask = FALSE)
  
  
  filename <- paste("./videoframes/Video1_",
                    rotation,
                    formatC(l, digits = 2, flag = "0"),
                    ".png", sep = "")
  
  rgl.snapshot(filename)
  
  print(l)
  
  close3d()
}


nopen3d() 
par3d(windowRect = c(0, 0, 800, 800)) 
par3d(userMatrix = um2)
par3d(zoom=0.7)
draw_everything()

rotation = 400

for (l in 1:90){
  #rotate in a loop (with l e.g. 1:90 for a 180 turn)
  nview3d(userMatrix = um2 %*%rotationMatrix(pi*l/90, 0, 0, 1)
          %*%rotationMatrix(-0.3*l/90, 0, 1, 0)
          %*%rotationMatrix(0.2*l/90, 1, 0, 0)
          )
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

# remove outline ----------

remove_outline()

plot3d(Q1, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.1, col = Okabe_Ito[1])
plot3d(Q2, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.1, col = Okabe_Ito[2])
plot3d(Q3, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.1, col = Okabe_Ito[6])
plot3d(Q4, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.1, col = Okabe_Ito[7])
add_outline()

rotation=500

for (l in 1:90){
  #rotate in a loop (with l e.g. 1:90 for a 180 turn)
  nview3d(userMatrix = um3 %*%rotationMatrix(pi*l/90, 0, 0, 1)
          %*%rotationMatrix(-0.3*l/90, 0, 1, 0)
          %*%rotationMatrix(0.2*l/90, 1, 0, 0)
  )
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
  output = 'manuscript/videos/Video1.3.mp4'
  )

unlink("videoframes", recursive = T)


