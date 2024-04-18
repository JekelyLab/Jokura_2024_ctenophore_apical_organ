rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memory and report the memory usage.

# load some packages
library(catmaid)
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

# catmaid connection, needs username, password AND token - weird!
{
  # can run this separate file using source function
  conn <- source("~/R/conn.R")
  #for larger calls we need to use http/1, see https://www.gitmemory.com/issue/natverse/rcatmaid/158/641537466
  #for this we configure to http/1.1
  conn_http1 = catmaid_login(conn=conn, config=httr::config(ssl_verifypeer=0, http_version=1))
}

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

skids_by_2annotations <- function(annotation1,annotation2){
  skids1 <- catmaid_skids(annotation1, pid = 35)
  skids2 <- catmaid_skids(annotation2, pid = 35)
  intersect <- intersect(skids1, skids2)
  return(intersect)
}

plot_background <- function(x){
  nopen3d() # opens a pannable 3d window
  par3d(zoom=0.75)
  nview3d("frontal", extramat=rotationMatrix(1.2, 0, 0, 1))
  par3d(windowRect = c(0, 0, 800, 800)) #resize for frontal view
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
  53, rval = c("mesh3d", "catmaidmesh", "raw"),
  invertFaces = T, conn = NULL, pid = 35
  )
dome_cavity <- catmaid_get_volume(
  54, rval = c("mesh3d", "catmaidmesh", "raw"),
  invertFaces = T, conn = NULL, pid = 35
  )

# define views --------------

anterior <- function(){
  nview3d("anterior", 
        extramat = rotationMatrix(2.54, 0.1, 0, 1)
        )
}

left <- function(){
  nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))
}



# crop from catmaid -----------------------------------------

crop_catmaid <- function(tagname,
                         half_bb_size_x, half_bb_size_y, half_bb_size_z,
                         zoomlevel,
                         dest_dir) {
  # tag search only returns treenodes with tag, not connectors
  # however, we want to get either treenodes or connectors, so we have to use generic search
  tagname_edited <- gsub(" ", "%20", tagname)
  search_results <- catmaid_fetch(path = paste("35/search?substring=", tagname_edited, sep = ""))
  # search doesn't accept regex, so we have to filter results
  ids=c() # keep track of the number of substacks we have to download from the server later
          # use ids rather than simple counter because it helps with debugging
  for (i in seq_along(search_results)) {
    if (search_results[[i]]$name == tagname && search_results[[i]]$class_name == "label") {
      for (j in seq_along(search_results[[i]]$nodes)) {
        #print("adding treenode id")
        #print(search_results[[i]]$nodes[[j]]$id)
        ids <- c(ids, search_results[[i]]$nodes[[j]]$id)
        x <- search_results[[i]]$nodes[[j]]$x
        y <- search_results[[i]]$nodes[[j]]$y
        z <- search_results[[i]]$nodes[[j]]$z
        catmaid_fetch(
          path = "35/crop/",
          body = list(
            stack_ids=28,
            min_x=x-half_bb_size_x,
            min_y=y-half_bb_size_y,
            min_z=z-half_bb_size_z,
            max_x=x+half_bb_size_x,
            max_y=y+half_bb_size_y,
            max_z=z+half_bb_size_z,
            zoom_level=zoomlevel,
            single_channel='true',
            rotationcw=0
          )
        )
      }
      for (j in seq_along(search_results[[i]]$connectors)) {
        #print("adding connector id")
        #print(search_results[[i]]$connector[[j]]$id)
        ids <- c(ids, search_results[[i]]$connectors[[j]]$id)
        x <- search_results[[i]]$connectors[[j]]$x
        y <- search_results[[i]]$connectors[[j]]$y
        z <- search_results[[i]]$connectors[[j]]$z
        catmaid_fetch(
          path = "35/crop/",
          body = list(
            stack_ids=28,
            min_x=x-half_bb_size_x,
            min_y=y-half_bb_size_y,
            min_z=z-half_bb_size_z,
            max_x=x+half_bb_size_x,
            max_y=y+half_bb_size_y,
            max_z=z+half_bb_size_z,
            zoom_level=zoomlevel,
            single_channel='true',
            rotationcw=0
          )
        )
      }
    }
  }
  print("checkpoint 1")
  if (length(ids) == 0) {
    print(paste("No treenodes or connectors tagged with", tagname, "were found"))
    return()
  }
  
  print(paste("Found", length(ids), "treenodes or connectors tagged with", tagname))
  print("Downloading the following ids:")
  print(ids)
  print("This will take some time")
        
  # server needs time to process this
  Sys.sleep(90)
  
  # download stacks ---------------------------------------------------------
  # the exact file names are needed for download,
  # to get them we need to check messages from the server
  # originally I was checking which messages are new, and downloading those
  # but this could lead to false positives
  # so I decided that's it's probably safe to download the last n messages
  # where n is the number of jobs
  
  # It would be cool if I can put the synapse ID as file name,
  # but it's not guaranteed that jobs will return in the same order
  
  cat_messages <- catmaid_fetch(
    path = "messages/list"
  )
  
  for (k in seq_along(ids)) {
    cat_message <- cat_messages[[k]]$text
    link <- regmatches(cat_message, 
                       regexpr("/crop/download/crop_.{6}.tiff", cat_message))
    filename <- regmatches(link, 
                           regexpr("crop_.{6}.tiff", link))
    full_link <- paste("https:/catmaid.ex.ac.uk", link, sep = "")
    destpath <- paste(dest_dir, "/", filename, sep = "")
    download.file(full_link, destfile = destpath)
  }
}


get_celltypes <- function(pid) {
  annotations <- catmaid_get_annotationlist(pid = 35)
  celltypes <- annotations$annotations |> filter(grepl("^celltype:", name)) |>
    select(name) |> pull() |> sub(pattern="^celltype:", replacement="")
  return(celltypes)
}
