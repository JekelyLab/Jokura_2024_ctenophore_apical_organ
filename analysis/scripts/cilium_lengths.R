# code to measure length of cilia
# WARNING: It is absolutely essential that the root node is proximal to the basal body 
# i.e., the order of nodes needs to be:
# "root"  -> (optional "exit_ciliary_pocket") -> "basal body" -> "cilium tip"

source("analysis/scripts/packages_and_functions.R")

tags <- catmaid_get_label_stats(pid = 35)

ctip_skids <- tags[tags$labelName == "cilium tip", "skeletonID"] %>% unique()
cil_neurons <- read.neurons.catmaid(ctip_skids, pid = 35, fetch.annotations = TRUE)
annotations <- attr(cil_neurons, 'anndf') |> as.tibble()

cilium_lengths <- tibble(celltype=character(),
                         skid=character(),
                         cil_length=integer(),
                         pocket_length=integer(),
                         axoneme=character())

for (cell in cil_neurons) {
  sigma = 2000
  sskid <- cell$skid
  # get celltype annotation
  celltype <- annotations |> 
    filter(skid==sskid) |> 
    filter(grepl("celltype:", annotation)) |> 
    select(annotation) |> 
    pull()
  celltype <- gsub(".*:","", celltype)
  if (length(celltype)==0) {
    celltype <- "NA"
  }
  # get cilium lengths
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
    if (length(cilium$tags[["9+2"]])) {
      axoneme <- "9+2"
    } else if (length(cilium$tags[["9+0"]])) {
      axoneme <- "9+0"
    } else if (length(cilium$tags[["9+unclear"]])) {
      axoneme <- "9+unclear"
    } else {
      axoneme <- "NA"
    }
    cil_len_tb <- tibble(celltype=as.character(celltype),
                         skid=as.character(sskid),
                         cil_length=as.integer(cil_length),
                         pocket_length=as.integer(pocket_length),
                         axoneme=as.character(axoneme))
    
    cilium_lengths <- bind_rows(cilium_lengths, cil_len_tb)
  }
}

write.csv(cilium_lengths, "analysis/data/cilium_lengths.csv")


# visualize lamellate body cells -----------------------------------------------
lamellate_intra <- read.neurons.catmaid("celltype:lamellate-intra", pid = 35)
for (cell in lamellate_intra) {
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
  lamellae <- segments_between_tags(cell, "lamellae end", "lamellae start")
  
  cell_smooth <- smooth_neuron(cell, sigma = 500)
  cilia_smooth <- lapply(cilia, smooth_neuron, sigma = 500) %>% as.neuronlist()
  in_pocket_smooth <- lapply(in_pocket, smooth_neuron, sigma = 500) %>% as.neuronlist()
  lamellae_smooth <- lapply(lamellae, smooth_neuron, sigma = 500) %>% as.neuronlist()
  
  plot3d(cell_smooth, soma = TRUE, lwd = 1, color = "grey")
  plot3d(cilia_smooth, soma = TRUE, lwd = 2, color = "black")
  plot3d(in_pocket_smooth, soma = TRUE, lwd = 5, color = "green")
  plot3d(lamellae_smooth, soma = TRUE, lwd = 10, color = "purple")
}

lamellate_extra <- read.neurons.catmaid("celltype:lamellate-extra", pid = 35)
for (cell in lamellate_extra) {
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
  lamellae <- segments_between_tags(cell, "lamellae end", "lamellae start")
  
  cell_smooth <- smooth_neuron(cell, sigma = 500)
  cilia_smooth <- lapply(cilia, smooth_neuron, sigma = 500) %>% as.neuronlist()
  in_pocket_smooth <- lapply(in_pocket, smooth_neuron, sigma = 500) %>% as.neuronlist()
  lamellae_smooth <- lapply(lamellae, smooth_neuron, sigma = 500) %>% as.neuronlist()
  
  plot3d(cell_smooth, soma = TRUE, lwd = 1, color = "grey")
  plot3d(cilia_smooth, soma = TRUE, lwd = 2, color = "black")
  plot3d(in_pocket_smooth, soma = TRUE, lwd = 5, color = "yellow")
  plot3d(lamellae_smooth, soma = TRUE, lwd = 10, color = "red")
}

# biciliated cell axoneme comparison
cilium_lengths %>%
  filter(celltype=="biciliated") %>%
  ggplot() +
  aes(x=axoneme, y=cil_length) +
    geom_violin() +
    geom_jitter(shape=16, position=position_jitter(0.2)) +
    theme_minimal() +
    theme(
      panel.grid.major = element_line(linewidth = 0.2),
      panel.background = element_blank(),
      axis.ticks = element_line(size = 0.2)
    ) +
    scale_fill_manual(
      values = c("#CC79A7", "#0072B2", "#E69F00", "grey40")
    )

