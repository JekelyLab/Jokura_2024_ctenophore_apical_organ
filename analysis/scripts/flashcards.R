source("analysis/scripts/packages_and_functions.R")

celltype <- "balancer"

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
  sigma = 2000
  cells <- read.neurons.catmaid(paste("celltype:", celltype, sep=""), pid=35)
  lapply(cells, plot_multinucleated_cell, sigma=sigma, color="darkblue", lwd=1)
  # define view
  anterior()
  par3d(zoom=1.1)
  par3d(windowRect = c(0, 0, 800, 800)) #resize for frontal view
  # draw lines for quadrants
  
  ### morphollogy for single cell, with cilia
  cell <- cells[[1]]
  cilia <- segments_between_tags(cell, "cilium tip", "basal body")
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
  in_pocket_smooth <- lapply(in_pocket, smooth_neuron, sigma = 500) %>% as.neuronlist()
  
  plot3d(cell_smooth, soma = TRUE, lwd = 1, color = "grey")
  if (length(cilia_smooth) > 0) {
    plot3d(cilia_smooth, soma = TRUE, lwd = 2, color = "orange")
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
    sskid <- cell$skid
    
    ########## cilia and pocket lengths
    cilium_lengths <- tibble(cil_length=integer(),
                             pocket_length=integer())
    cilia <- segments_between_tags(cell, "cilium tip", "basal body")
    for (cilium in cilia) {
      cil_length <- smooth_neuron(cilium, sigma=sigma) %>%
        summary() %>% select(cable.length)
      # check if it has exit_ciliary_pocket tag
      #print(cilium$tags)
      exit_cil_pocket <- cilium$tags[["exit_ciliary_pocket"]]
      ctip_intra <- cilium$tags[["cilium tip intracellular"]]
      if (length(exit_cil_pocket) > 0) {
        pocket <- segments_between_tags(cilium, "exit_ciliary_pocket", "basal body")
        pocket_length <- smooth_neuron(pocket[[1]], sigma = sigma) %>%
          summary() %>% select(cable.length)
      }
      else if (length(ctip_intra) > 0) {
        pocket_length <- cil_length
      }
      else {
        pocket_length <- 0
      }
    }
    cilium_lengths <- tibble(cil_length=integer(),
                             pocket_length=integer())
    cilium_lengths <- bind_rows(cilium_lengths, cil_len_tb)
    
    ########## polyadicity 
    polyadicities <- list()
    polyadicity <- get_syn_polyadicity_for_neuron(cell)
    polyadicities <- c(polyadicities, polyadicity)
  }
  avg_celltype_polyadicity <- polyadicities %>% unlist %>% mean
}
