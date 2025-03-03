source("analysis/scripts/packages_and_functions.R")

celltype <- "balancer"
cilium_lengths <- read_csv("analysis/data/cilium_lengths.csv")
organelle_stats <- read_csv("analysis/data/organelle_stats.csv")
stats_master <- read_csv("analysis/data/stats_master.csv")

celltype_flashcard(celltype) {
  # things to show:
  # - location by quadrant
  # - morphology of example cell (with cilia)
  # - EM
  # - number of cilia
  # - number of centrioles
  # - cilium length
  # - ciliary pocket depth
  # - number of mitochondria
  # - number of pre synapses
  # - number of post synapses
  # - synapse polyadicity
  # - pre partner celltypes
  # - post partner celltypes
  
  ### location by quadrant
  # read, smooth and plot skeletons
  cells <- read.neurons.catmaid(paste("celltype:", celltype, sep=""), pid=35)
  skids <- names(cells)
  
  sigma = 2000
  lapply(cells, plot_multinucleated_cell, sigma=sigma, color="darkblue", lwd=1)
  # define view
  anterior()
  par3d(zoom=1.1)
  par3d(windowRect = c(0, 0, 800, 800)) #resize for frontal view
  # draw lines for quadrants
  
  ### morphology for single cell, with cilia
  cell <- cells[[1]]
  cilia <- segments_between_tags(cell, "cilium tip", "basal body")
  # here we also want cut cilia
  cilia_cut <- segments_between_tags(cell, "cilium cut", "basal body")
  # visualize part of cilium in ciliary pocket
  # this is a bit tricky, because it has "exit_ciliary_pocket" if it exits half way,
  # but "cilium tip intracellular" if the whole cilium is inside the cell
  in_pocket <- neuronlist()
  for (cilium in cilia) {
    exit_cil_pocket <- cilium$tags[["exit_ciliary_pocket"]]
    ctip_intra <- cilium$tags[["cilium tip intracellular"]]
    if (length(exit_cil_pocket) > 0) {
      pocket <- segments_between_tags(cilium, "exit_ciliary_pocket", "basal body")
    } else if (length(ctip_intra) > 0) {
      pocket <- cilium
    } else {
      pocket <- neuronlist()
    }
    in_pocket <- c(in_pocket, as.neuronlist(pocket))
  }
  
  cell_smooth <- smooth_neuron(cell, sigma = 500)
  cilia_smooth <- lapply(cilia, smooth_neuron, sigma = 500) %>% as.neuronlist()
  cilia_cut_smooth <- lapply(cilia_cut, smooth_neuron, sigma = 500) %>% as.neuronlist()
  in_pocket_smooth <- lapply(in_pocket, smooth_neuron, sigma = 500) %>% as.neuronlist()
  plot3d(cell_smooth, soma = TRUE, lwd = 1, color = "grey")
  if (length(cilia_smooth) > 0) {
    plot3d(cilia_smooth, soma = TRUE, lwd = 2, color = "orange")
  }
  if (length(cilia_cut_smooth) > 0) {
    plot3d(cilia_cut_smooth, soma = TRUE, lwd = 2, color = "orange")
  }
  if (length(in_pocket_smooth) > 0) {
    plot3d(in_pocket_smooth, soma = TRUE, lwd = 5, color = "purple")
  }
  # define view
  left()
  par3d(zoom=1.1)
  par3d(windowRect = c(0, 0, 800, 800)) #resize for frontal view
  
  ### EM images
  # get locations of "flashcard" tags
  # we don't know which skeleton has flashcard tag
  
  for (cell in cells) {
    if (length(cell$tags$'flashcard:distal') == 1) {
      flashcard_distal_treenode <- cell$tags$'flashcard:distal'
    }
    if (length(cell$tags$'flashcard:proximal') == 1) {
      flashcard_proximal_treenode <- cell$tags$'flashcard:proximal'
    }
    if (length(cell$tags$'flashcard:soma') == 1) {
      flashcard_soma_treenode <- cell$tags$'flashcard:soma'
    }
  }
  
  # crop distal
  pos <- which(cell$d$PointNo == flashcard_distal_treenode)
  x <- cell$d$X[pos]
  y <- cell$d$Y[pos]
  z <- cell$d$Z[pos]
  crop_substack_point(x, y, z, 500, 500, 0, 0,
                      paste("manuscript/pictures/", celltype, "flashcard_distal.tif"),
                      35, 28)
  
  # crop proximal
  pos <- which(cell$d$PointNo == flashcard_proximal_treenode)
  x <- cell$d$X[pos]
  y <- cell$d$Y[pos]
  z <- cell$d$Z[pos]
  crop_substack_point(x, y, z, 500, 500, 0, 0,
                      paste("manuscript/pictures/", celltype, "flashcard_proximal.tif"),
                      35, 28)
  # crop soma
  pos <- which(cell$d$PointNo == flashcard_soma_treenode)
  x <- cell$d$X[pos]
  y <- cell$d$Y[pos]
  z <- cell$d$Z[pos]
  crop_substack_point(x, y, z, 500, 500, 0, 0,
                      paste("manuscript/pictures/", celltype, "flashcard_soma.tif"),
                      35, 28)
  
  
  n_cilia <- organelle_stats %>%
    filter(skid %in% skids) %>%
    select("basal body") %>%
    pull() %>%
    mean()
  skids_with_cilia <- cilium_lengths %>% filter(skid %in% skids)
  # if there is no basal body, report 0
  # cilium_lengths table only tells us how many cilium tips there are
  if (n_cilia == 0) {
    avg_cil_length <- 0
  }
  # if there is no "cilium tip" report NA
  else if (nrow(skid_with_cilia) == 0) {
    avg_cil_length <- NA
  } else {
    avg_cil_length <- cilium_lengths %>%
      filter(skid %in% skids) %>%
      select(cil_length) %>%
      pull() %>% 
      mean()
  }

  if (n_cilia == 0) {
    avg_pocket_length <- 0
  }
  # if there is no "cilium tip" report NA
  else if (nrow(skid_with_cilia) == 0) {
    avg_pocket_length <- NA
  } else {
    avg_cil_length <- cilium_lengths %>%
      filter(skid %in% skids) %>%
      select(pocket_length) %>%
      pull() %>% 
      mean()
  }
  
  ########### mito count
  n_mito <- organelle_stats %>%
    filter(skid %in% skids) %>%
    select(mitochondrion) %>%
    pull %>%
    mean()
  
  ######### polyadicity 
  avg_polyadicity <- stats_master %>%
   filter(skid %>% skids) %>%
   select(polyadicity) %>%
   pull() %>%
   mean()
}




crop_substack_point <- function(x, y, z,
                                half_bb_size_x, half_bb_size_y, half_bb_size_z,
                                zoomlevel,
                                dest_file_path,
                                pid, stackid) {
  # this only has one coordinate is input, so we can define filename
  catmaid_fetch(
    path = paste(pid, "/crop/", sep = ""),
    body = list(
      stack_ids=stackid,
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
  
  print("Downloading, this will take some time")
  
  # server needs time to process this
  Sys.sleep(61)
  
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
    download.file(full_link, destfile = dest_file_path)
  }
}
