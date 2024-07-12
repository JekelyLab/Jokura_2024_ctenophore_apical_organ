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
                         pocket_length=integer())

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
    cil_len_tb <- tibble(celltype=as.character(celltype),
                         skid=as.character(sskid),
                         cil_length=as.integer(cil_length),
                         pocket_length=as.integer(pocket_length))
    
    cilium_lengths <- bind_rows(cilium_lengths, cil_len_tb)
  }
}

write.csv(cilium_lengths, "analysis/data/cilium_lengths.csv")
