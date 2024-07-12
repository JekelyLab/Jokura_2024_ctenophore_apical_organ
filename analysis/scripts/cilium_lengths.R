# code to measure length of cilia
# WARNING: It ia absolutely essential that the root node is proximal to the basal body 
# i.e., the order of nodes needs to be "root" -> "basal body" -> "cilium tip"
# it should work fine if skeleton has multiple cilia

library(catmaid)

conn <- source("~/R/conn.R")

tags <- catmaid_get_label_stats(pid = 35)

# using "cilium tip" tag because that shows that cilium is complete, cut off cilia don't have that tag
# measuring from tips because it's easy to get parent node, but harder to get child
ctips <- tags[tags$labelName == "cilium tip", "treenodeID"]

cilium_length <- function(ctip_id, smooth_sigma = 1) {
  # this will only work if root node is proximal to basal body tag
  # treenode compact detail response:
  # [ID, parent ID, x, y, z, confidence, radius, skeleton_id, edition_time, user_id]
  path <- paste("/35/treenodes/", ctip_id, "/compact-detail", sep = "")
  ctip_treenode <- catmaid_fetch(path = path)
  skid <- ctip_treenode[[8]]
  neuron_info <- smooth_neuron(read.neuron.catmaid(skid, pid = 35), sigma=smooth_sigma)
  neuron_name <- catmaid_fetch(path = "/35/skeleton/neuronnames",
                               body = list(skids=skid))[[1]]
  # check which treenodes have basal body annotation, because this is the stop point
  basalb_ids <- neuron_info$tags$`basal body`
  treenode_id <- ctip_id
  distance <- 0
  while (!(treenode_id %in% basalb_ids)) {
    index <- match(treenode_id, neuron_info$d$PointNo)
    parent_id <- neuron_info$d$Parent[[index]]
    if (parent_id == -1) {
      print("root reached before basal body tag, from cilium tip")
      print(ctip_id)
      stop()
    }
    index_parent <- match(parent_id, neuron_info$d$PointNo)
 
    distance_2nodes <- dist(rbind(c(neuron_info$d$X[[index]],
                                neuron_info$d$Y[[index]],
                                neuron_info$d$Z[[index]]),
                              c(neuron_info$d$X[[index_parent]],
                                neuron_info$d$Y[[index_parent]],
                                neuron_info$d$Z[[index_parent]])
                                ))
    distance <- distance + distance_2nodes
    #print(distance)
    treenode_id <- parent_id
  }
  return(c(distance, length(basalb_ids), neuron_name))
}


cilium_lengths_monociliated <- list()
cilium_lengths_biciliated <- list()
cilium_lengths_intramulti <- list()
cilium_lengths_lamellate <- list()
for (ctip_id in ctips) {
  c_length_type <- c(cilium_length(ctip_id, 250))
  cilium_length_new <- c_length_type[[1]]
  nbb <- c_length_type[[2]]
  nname <- c_length_type[[3]]
  names(cilium_length_new) <- ctip_id
  if (nbb == 1) {
    cilium_lengths_monociliated <- c(cilium_lengths_monociliated, cilium_length_new)
  }
  if (nbb == 2) {
    cilium_lengths_biciliated <- c(cilium_lengths_biciliated, cilium_length_new)
  }
  if (nbb > 2) {
    #print(nname)
    #print(ctip_id)
    if (grepl("intracellular multiciliated cell", nname)) {
      cilium_lengths_intramulti <- c(cilium_lengths_intramulti, cilium_length_new)
    }
    if (grepl("lamellate", nname)) {
      cilium_lengths_lamellate <- c(cilium_lengths_lamellate, cilium_length_new)
    }
  }
}

hist(as.numeric(cilium_lengths_monociliated))
hist(as.numeric(cilium_lengths_biciliated))
hist(as.numeric(cilium_lengths_intramulti))
hist(as.numeric(cilium_lengths_lamellate))
t.test(as.numeric(cilium_lengths_monociliated), as.numeric(cilium_lengths_biciliated))
 


# --------------------------------------------------------------
neurons_info <- read.neurons.catmaid("intracellular multiciliated cell", pid = 35)
neurons_info <- read.neurons.catmaid("lamellate body", pid = 35)
n_extra <- list()
for (neuroninf in neurons_info) {
  n_extra <- c(n_extra, length(neuroninf$tags$`cilium tip extracellular`))
}
summary(as.numeric(n_extra))
n_intra <- list()
for (neuroninf in neurons_info) {
  n_intra <- c(n_intra, length(neuroninf$tags$`cilium tip intracellular`))
}
summary(as.numeric(n_intra))


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
    exit_cil_pocket <- cilium$tags[["exit_ciliary_pocket"]]
    ctip_intra <- cilium$tags[["cilium tip instracellular"]]
    if (length(exit_cil_pocket) > 0) {
      pocket <- segments_between_tags(cilium, "exit_ciliary_pocket", "basal body")
      pocket_length <- smooth_neuron(pocket, sigma = sigma) %>%
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
