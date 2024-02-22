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

monociliated <- nlapply(
  read.neurons.catmaid(
    "monociliated_cell", pid = 35
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
AO_celltypes <- list(
  balancer, LB, syn_neuron, bridge, monociliated
  )
AO_celltype_names <- list(
  "balancer", "LB", "syn_neuron", "bridge", "monociliated"
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
# convert synapse list into a matrix of appropriate dimensions
synapse_matrix <- matrix(
  unlist(synapse_list), byrow = TRUE, 
  nrow = length(AO_celltypes)
  )
synapse_matrix

rownames(synapse_matrix) <- as.character(AO_celltype_names)
colnames(synapse_matrix) <- as.character(AO_celltype_names)
synapse_matrix

# with the make_graph function of igraph we turn it into a graph (input is the list of edge pairs)
AO_graph <- graph_from_adjacency_matrix(
  synapse_matrix,
  mode = c("directed"),
  weighted = TRUE, diag = TRUE
)

AO_graph |> as_tbl_graph() 


# calculate node weighted degree -------------
degree=degree(
  AO_graph, v = V(AO_graph), mode = c("all"), 
  loops = TRUE, normalized = FALSE
  )
degree

# use visNetwork to plot the network --------------------------------------

## convert to VisNetwork-list
AO_graph.visn <- toVisNetworkData(AO_graph)

## copy column "weight" to new column "value" in list "edges"
Conn_graph.visn$edges$value <- Conn_graph.visn$edges$weight
Conn_graph.visn$nodes$value = degree
Conn_graph.visn$nodes$group <- cell_group_attr$type


#hierarchical layout

#level	: Number. Default to undefined. When using the hierarchical layout, the level determines where the node is going to be positioned.
Conn_graph.visn$nodes$level <- cell_group_attr$level
#hierarchical layout
{
  visNet <- visNetwork(Conn_graph.visn$nodes,Conn_graph.visn$edges) %>% 
    visIgraphLayout(layout = "layout_nicely", physics = TRUE, 
                    randomSeed = 42, type="square") %>%
    visHierarchicalLayout(levelSeparation=250, 
                          nodeSpacing=10,
                          direction='LR',
                          sortMethod='hubsize',
                          shakeTowards='roots') %>%
    visEdges(smooth = list(type = 'curvedCW', roundness=0.1),
             scaling=list(min=2, max=12),
             color = list(inherit=TRUE, opacity=0.7),
             arrows = list(to = list(enabled = TRUE, 
                                     scaleFactor = 1.2, type = 'arrow'))) %>%
    visNodes(borderWidth=0.3, 
             color = list(background=Conn_graph.visn$nodes$color, border='black'),
             opacity=0.9,
             shape='dot', 
             font=list(color='black', size=44),
             scaling = list(label=list(enabled=TRUE, min=48, max=56)),
             level= Conn_graph.visn$nodes$level) %>%
    visOptions(highlightNearest = list(enabled=TRUE, degree=1, algorithm='hierarchical',labelOnly=FALSE)) %>%
    visInteraction(dragNodes = TRUE, dragView = TRUE,
                   zoomView = TRUE, hover=TRUE,
                   multiselect=TRUE) %>%
    visGroups(groupname = "MBSN", shape = "square", 
              opacity=1, color="#E69F00") %>%
    visGroups(groupname = "MBintrIN", shape = "square", 
              opacity=1, color="#CC79A7") %>%
    visGroups(groupname = "MBON", shape = "square", 
              opacity=1, color="#56B4E9") %>%
    visGroups(groupname = "IN", shape = "dot", 
              opacity=1, color="#0072B2") %>%
    visGroups(groupname = "SN", shape = "dot", 
              opacity=1, color="#D55E00") %>%
    visGroups(groupname = "MN", shape = "dot", 
              opacity=1, color="#cccccc")  %>%
    addFontAwesome() %>%
    visLegend(addNodes = list(
      list(label = "66 syn", shape = "icon", 
           icon = list(code = "f2d1", size = 30, color = "#D55E00"))), 
      useGroups = TRUE,  width=0.1,ncol = 1,
      position='right', stepY=70)
  
  
  visNet
}
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