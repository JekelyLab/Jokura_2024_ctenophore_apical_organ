# Code to generate Figure 1 of the Jokura et al 2024 Ctenophore apical organ connectome paper

# source packages and functions ------------------------------------------------
source("analysis/scripts/packages_and_functions.R")





# load cell types from catmaid -------------------------------------------------

balancer <-
  nlapply(read.neurons.catmaid("celltype:balancer", pid = 35),
          function(x)
            smooth_neuron(x, sigma = 1000))

bridge <-
  nlapply(read.neurons.catmaid("celltype:bridge", pid = 35),
          function(x)
            smooth_neuron(x, sigma = 1000))

bristle <-
  nlapply(read.neurons.catmaid("celltype:bristle", pid = 35),
          function(x)
            smooth_neuron(x, sigma = 1000))

dome <-
  nlapply(read.neurons.catmaid("celltype:dome", pid = 35),
          function(x)
            smooth_neuron(x, sigma = 1000))

groove <-
  nlapply(read.neurons.catmaid("celltype:groove", pid = 35),
          function(x)
            smooth_neuron(x, sigma = 1000))

intramulticilia <-
  nlapply(read.neurons.catmaid("celltype:intra-multi-ciliated", pid = 35),
          function(x)
            smooth_neuron(x, sigma = 1000))

lamellate <-
  nlapply(read.neurons.catmaid("celltype:lamellate", pid = 35),
          function(x)
            smooth_neuron(x, sigma = 1000))

lithocyte <-
  nlapply(read.neurons.catmaid("celltype:lithocyte", pid = 35),
          function(x)
            smooth_neuron(x, sigma = 1000))

neuron <-
  nlapply(read.neurons.catmaid("celltype:neuron", pid = 35),
          function(x)
            smooth_neuron(x, sigma = 1000))

plumose <-
  nlapply(read.neurons.catmaid("celltype:plumose", pid = 35),
          function(x)
            smooth_neuron(x, sigma = 1000))

dense_vesicle <-
  nlapply(read.neurons.catmaid("celltype:dense_vesicle", pid = 35),
          function(x)
            smooth_neuron(x, sigma = 1000))

monocilia <-
  nlapply(read.neurons.catmaid("celltype:monociliated", pid = 35),
          function(x)
            smooth_neuron(x, sigma = 1000))

bicilia <-
  nlapply(read.neurons.catmaid("celltype:biciliated", pid = 35),
          function(x)
            smooth_neuron(x, sigma = 1000))

non_cilia <-
  nlapply(read.neurons.catmaid("celltype:nonciliated", pid = 35),
          function(x)
            smooth_neuron(x, sigma = 1000))



# cell types -------------------------------------------------------------------

celltypes <- list(balancer, bridge, bristle, dome, groove, intramulticilia, 
                     lamellate, lithocyte, neuron, plumose, dense_vesicle, 
                     monocilia, bicilia, non_cilia)

celltype_names <- list("balancer", "bridge", "bristle", "dome", "groove", 
                          "intramulticilia", "lamellate", "lithocyte", "neuron", 
                          "plumose", "dense_vesicle", "monocilia", "bicilia", "non_cilia")


# iterate through cell group neuron lists and get connectivity
# define empty synapse list with the right dimensions

synapse_list <- c()

for (i in 1:length(celltypes)) {
  for (j in 1:length(celltypes)) {
    # get connectors between two cell groups
    presyn_skids <- attr(celltypes[i][[1]], "df")$skid
    postsyn_skids <- attr(celltypes[j][[1]], "df")$skid
    connectivity <- catmaid_get_connectors_between(
      pre = presyn_skids,
      post = postsyn_skids, pid = 35
    )
    # check the number of synapses from group1 -> group2
    N_synapses <- dim(connectivity)[1]
    if(length(connectivity) == 0) {N_synapses = 0}
    # add value to synapse list
    synapse_list <- c(synapse_list, N_synapses)
  }
}
synapse_list

# convert synapse list into a matrix of appropriate dimensions
synapse_matrix <- matrix(
  unlist(synapse_list), byrow = TRUE, 
  nrow = length(celltypes)
)
synapse_matrix

rownames(synapse_matrix) <- as.character(celltype_names)
colnames(synapse_matrix) <- as.character(celltype_names)
synapse_matrix

# with the make_graph function of igraph we turn it into a graph (input is the list of edge pairs)
graph <- graph_from_adjacency_matrix(
  synapse_matrix,
  mode = c("directed"),
  weighted = TRUE, diag = TRUE
)

graph |> as_tbl_graph() 


# calculate node weighted degree ---------------------------------------
degree=degree(
  graph, v = V(graph), mode = c("all"), 
  loops = TRUE, normalized = FALSE
)
degree



# use visNetwork to plot the network --------------------------------------

## convert to VisNetwork-list
graph.visn <- toVisNetworkData(graph)

## copy column "weight" to new column "value" in list "edges"
graph.visn$edges$value <- graph.visn$edges$weight
graph.visn$nodes$value <- degree

#define node color
graph.visn$nodes$color <- Okabe_Ito[1:14]

#hierarchical layout - define level of nodes
graph.visn$nodes$level <- c(1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2)

#hierarchical layout
visNetwork(graph.visn$nodes, graph.visn$edges) %>%
  visIgraphLayout(
    layout = "layout_nicely", physics = TRUE, 
    randomSeed = 42, type="square"
  ) %>%
  visHierarchicalLayout(
    levelSeparation=250, 
    nodeSpacing=10,
    direction='LR',
    sortMethod='hubsize',
    shakeTowards='roots'
  ) %>%
  visEdges(
    smooth = list(type = 'curvedCW', roundness=0.2),
    scaling=list(min=2, max=12),
    color = list(inherit=TRUE, opacity=0.7),
    arrows = list(
      to = list(enabled = TRUE, 
                scaleFactor = 1, type = 'arrow'))
  ) %>%
  visNodes(
    borderWidth=0.3, 
    color = list(background=graph.visn$nodes$color, border='black'),
    opacity=0.9,
    shape='dot', 
    font=list(color='black', size=44),
    scaling = list(label=list(enabled=TRUE, min=48, max=56)),
    level= graph.visn$nodes$level
  )
