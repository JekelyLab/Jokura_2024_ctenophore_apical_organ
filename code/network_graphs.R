#network analysis and plotting

source("code/packages_and_functions.R")

# Create a random graph with 10 nodes
g <- erdos.renyi.game(10, p = 0.3)

# Assign names to the nodes
V(g)$name <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")
#edges
E(g)
#vertices (nodes)
V(g)

g_tb <- g |> 
  as_tbl_graph()

g_tb |>
  activate(nodes) |>
  filter(name == "A")

g_tb |>
  activate(nodes) |>
  mutate(color = "blue")

plot(g_tb)

# load cell types from catmaid

balancer <- nlapply(
  read.neurons.catmaid(
    "celltype:balancer", pid = 35
    ),
  function(x) smooth_neuron(x, sigma = 1000)
)

LB <- nlapply(
  read.neurons.catmaid(
    "celltype:lamellate", pid = 35
    ),
  function(x) smooth_neuron(x, sigma = 1000)
)

syn_neuron <- nlapply(
  read.neurons.catmaid(
    "celltype:neuron", pid = 35
    ),
  function(x) smooth_neuron(x, sigma = 1000)
)

bridge <- nlapply(
  read.neurons.catmaid(
    "celltype:bridge", pid = 35
  ),
  function(x) smooth_neuron(x, sigma = 1000)
)

# cell types
AO_celltypes <- list(balancer, LB, syn_neuron, bridge)
AO_celltype_names <- list(
  "balancer", "LB", "syn_neuron", "bridge"
  )

# iterate through cell group neuron lists and get connectivity
# define empty synapse list with the right dimensions
synapse_list <- c()

for (i in 1:length(AO_celltypes)) {
  for (j in 1:length(AO_celltypes)) {
    # get connectors between two cell groups
    presyn_skids <- attr(AO_celltypes[i][[1]], "df")$skid
    postsyn_skids <- attr(AO_celltypes[j][[1]], "df")$skid
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

# create graph with cell types as nodes
node_IDs <- data.frame(name = cell_types)
graph.funct <- tbl_graph(nodes = node_IDs)

# define sources and targets
sources <- 
targets <- 

# edges 
graph_edges <- data.frame(
    from = rep(sources, each = length(targets)),
    to = rep(targets, length(sources)),
    value = ,
    color = ,
    type = 
  )
  
# add edges to graph
graph.funct.edges <- graph.funct %>%
    bind_edges(graph_edges)