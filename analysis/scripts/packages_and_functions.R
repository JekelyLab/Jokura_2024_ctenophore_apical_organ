rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memory and report the memory usage.

# load some packages
library(catmaid)
library(plyr)
library(tidyverse)
library(cowplot)
library(png)
library(igraph)
library(networkD3)
library(visNetwork)
library(webshot2)
library(patchwork)
library(RColorBrewer)
library(igraph)
library(tidygraph)
library(av)
library(jpeg)
library(magick)


# catmaid connection, needs username, password AND token - weird!
{
  # can run this separate file using source function
  conn <- source("~/R/conn.R")
  #for larger calls we need to use http/1, see https://www.gitmemory.com/issue/natverse/rcatmaid/158/641537466
  #for this we configure to http/1.1
  conn_http1 = catmaid_login(conn=conn, config=httr::config(ssl_verifypeer=0, http_version=1))
}

system("git submodule init")
system("git submodule update")
system("cd rcatmaid_functions_library && git pull origin main")
system("cd rcatmaid_functions_library && git checkout 626b317b8adfc27e1aaa0f138b767385a136edd2")
source("rcatmaid_functions_library/functions.R")

#define some colors
#From Color Universal Design (CUD): https://jfly.uni-koeln.de/color/
Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", 
               "#CC79A7", "#000000")

blues <- brewer.pal(9, 'Blues')
bluepurple <- brewer.pal(9, 'BuPu')
oranges <- brewer.pal(9, 'YlOrRd')

#save session info and Rstudio version info for reproducibility
writeLines(capture.output(sessionInfo()), "analysis/scripts/sessionInfo.txt")
writeLines(capture.output(rstudioapi::versionInfo()), "analysis/scripts/versionInfo.txt")

# functions ----------------------

# skids_by_2annotations has been replaced by get_skids_with_annot from rcatmaid_functions_library

plot_background <- function(x){
  nopen3d() # opens a pannable 3d window
  par3d(zoom=0.75)
  nview3d("frontal", extramat=rotationMatrix(1.2, 0, 0, 1))
  par3d(windowRect = c(0, 0, 800, 800)) #resize for frontal view
  nview3d("anterior", extramat = rotationMatrix(1.05, 250, -200, 1000))
  par3d(zoom=0.7)
  #y-axis clip
  clipplanes3d(1, 0, 0, -11500)
  #x-axis clip
  clipplanes3d(0, 1, 0, -24000)
  
}

plot_background_ventral <- function(x){
  nopen3d() # opens a pannable 3d window
  par3d(zoom=0.75)
  nview3d("ventral", extramat=rotationMatrix(0.9, 0.3, -0.1, 1))
  par3d(windowRect = c(0, 0, 800, 800)) #resize for frontal view
}

read_smooth_neuron <- function(annotation){
  nlapply(
  read.neurons.catmaid(
    annotation, pid = 35
  ),
  function(x) smooth_neuron(x, sigma = 1000)
  )
}

outline <- catmaid_get_volume(
  60, rval = c("mesh3d", "catmaidmesh", "raw"),
  invertFaces = T, conn = NULL, pid = 35
  )
dome_cavity <- catmaid_get_volume(
  54, rval = c("mesh3d", "catmaidmesh", "raw"),
  invertFaces = T, conn = NULL, pid = 35
  )

# define views --------------

anterior <- function(){
  nview3d("anterior", extramat = rotationMatrix(1.05, 250, -200, 1000))
}

sagittal <- function(){
  nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
}
tentacular <- function(){
  nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))
}

# crop_catmaid function has been replaced by crop_substack from rcatmaid_functions_library

# get_celltypes function is in rcatmaid_functions_library


#function to add a text label to a neuron and a line to one of the neurons of the neuronlist starting from the soma
add_label_with_line <- function(neuronlist, offsetx, offsety, offsetz, side, label){
  cords_left = unlist(as_tibble(soma(neuronlist)) %>%
                        arrange(X) %>% slice(1) %>% as.list()
  )
  cords_right = unlist(as_tibble(soma(neuronlist)) %>%
                         arrange(desc(X)) %>% slice(1) %>% as.list()
  )
  
  if (side == "left") {
    cords <- cords_left
  } else {
    cords <- cords_right
  }
  
  print(cords)
  lines3d(c(cords[1], cords[1]-offsetx), 
          c(cords[2], cords[2]-offsety), 
          c(cords[3], cords[3]-offsetz), lwd=2)
  adj_x <- if_else(offsetx < 0, 0, 1)
  adj_y <- if_else(offsety < 0, 1, 0)
  
  if (missing(label)) {
    label <- deparse(substitute(neuronlist))
  }
  texts3d(c(cords[1]-offsetx, cords[2]-offsety, 
            cords[3]-offsetz), 
          text = label, 
          adj = c(adj_x, adj_y, 1),
          cex = 1.7)
}
