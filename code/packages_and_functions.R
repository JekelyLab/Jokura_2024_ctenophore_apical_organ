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
writeLines(capture.output(sessionInfo()), "code/sessionInfo.txt")
writeLines(capture.output(rstudioapi::versionInfo()), "code/versionInfo.txt")

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
