#network analysis and plotting

source("analysis/scripts/packages_and_functions.R")

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
    "celltype:SSN", pid = 35
    ),
  function(x) smooth_neuron(x, sigma = 1000)
)

monociliated <- nlapply(
  read.neurons.catmaid(
    "celltype:monociliated", pid = 35
  ),
  function(x) smooth_neuron(x, sigma = 1000)
)

biciliated <- nlapply(
  read.neurons.catmaid(
    "celltype:biciliated", pid = 35
  ),
  function(x) smooth_neuron(x, sigma = 1000)
)

nonciliated <- nlapply(
  read.neurons.catmaid(
    "celltype:nonciliated", pid = 35
  ),
  function(x) smooth_neuron(x, sigma = 1000)
)

bridge <- nlapply(
  read.neurons.catmaid(
    "celltype:bridge", pid = 35
  ),
  function(x) smooth_neuron(x, sigma = 1000)
)

dense_vesicle <- nlapply(
  read.neurons.catmaid(
    "celltype:dense_vesicle", pid = 35
  ),
  function(x) smooth_neuron(x, sigma = 1000)
)

intra_multiciliated <- nlapply(
  read.neurons.catmaid(
    "celltype:intra-multi-ciliated", pid = 35
  ),
  function(x) smooth_neuron(x, sigma = 1000)
)

plumose <- nlapply(
  read.neurons.catmaid(
    "celltype:plumose", pid = 35
  ),
  function(x) smooth_neuron(x, sigma = 1000)
)



# cell types
AO_celltypes <- list(
  balancer, LB, syn_neuron, bridge, monociliated, biciliated, 
  nonciliated, dense_vesicle, intra_multiciliated, plumose
  )
AO_celltype_names <- list(
  "balancer", "lamellar body", "nerve net", "bridge", "monociliated", "biciliated", "nonciliated", "dense-vesicle", 
  "intra-multiciliated", "plumose"
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

rownames(synapse_matrix) <- as.character(AO_celltype_names)
colnames(synapse_matrix) <- as.character(AO_celltype_names)

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

# plot cells
plot_background()

# "#E69F00" "#56B4E9" "#009E73" "#F0E442" "#0072B2" "#D55E00" "#CC79A7" "#000000"

plot3d(
  balancer, soma = TRUE, color = Okabe_Ito[1], 
  alpha = 0.1, lwd = 3
  )

plot3d(
  LB, soma = TRUE, color = Okabe_Ito[2], 
  alpha = 0.4, lwd = 2
)

plot3d(
  syn_neuron, soma = TRUE, color = Okabe_Ito[c(1,5,8)], 
  alpha = 1, lwd = c(4,3,3)
)

plot3d(
  bridge, soma = TRUE, color = Okabe_Ito[4], 
  alpha = 0.4, lwd = 3
)

plot3d(
  monociliated, soma = TRUE, color = Okabe_Ito[5], 
  alpha = 0.1, lwd = 3
)

plot3d(
  biciliated, soma = TRUE, color = Okabe_Ito[6], 
  alpha = 0.2, lwd = 3
)

plot3d(
  nonciliated, soma = TRUE, color = Okabe_Ito[7], 
  alpha = 0.1, lwd = 3
)

plot3d(
  intra_multiciliated, soma = TRUE, color = Okabe_Ito[8], 
  alpha = 0.2, lwd = 3
)

plot3d(
  plumose, soma = TRUE, color = bluepurple[7], 
  alpha = 0.1, lwd = 3
)

nview3d("ventral", extramat=rotationMatrix(0.9, 0.3, -0.1, 1))
rgl.snapshot("manuscript/pictures/celltypes_in_network_lateral.png")

nview3d("frontal", extramat=rotationMatrix(2.1, 0, 0, 1)
        %*%rotationMatrix(3.14, 1, 0, 0))
rgl.snapshot("manuscript/pictures/celltypes_in_network_oral.png")
close3d()

# use visNetwork to plot the network --------------------------------------

## convert to VisNetwork-list
AO_graph.visn <- toVisNetworkData(AO_graph)

#filter low-weight edges
AO_graph.visn$edges$weight
AO_graph.visn$edges$weight[AO_graph.visn$edges$weight < 3]  <- 0
AO_graph.visn$edges$weight

## copy column "weight" to new column "value" in list "edges"
AO_graph.visn$edges$value <- AO_graph.visn$edges$weight
AO_graph.visn$nodes$value <- degree

#define node color
AO_graph.visn$nodes$color <- Okabe_Ito[1:10]

#hierarchical layout - define level of nodes
AO_graph.visn$nodes$level <- c(2, 4, 1, 1, 3, 3, 3, 3, 3, 4)
#balancer, LB, syn_neuron, bridge, monociliated, biciliated, nonciliated, dense_vesicle

#hierarchical layout
visNet <- visNetwork(AO_graph.visn$nodes, AO_graph.visn$edges) %>%
  visIgraphLayout(
    layout = "layout_nicely", physics = FALSE
    ) %>%
  visHierarchicalLayout(
    levelSeparation=250, 
    nodeSpacing=200,
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
    color = list(background=AO_graph.visn$nodes$color, border='black'),
    opacity=0.9,
    shape='dot', 
    font=list(color='black', size=44),
    scaling = list(label=list(enabled=TRUE, min=22, max=80)),
    level= AO_graph.visn$nodes$level
    )
#%>%
  visOptions(highlightNearest = TRUE, width = 1800, height = 800) %>%
  visInteraction(navigationButtons = FALSE,
           dragNodes = TRUE, dragView = FALSE,
           zoomView = TRUE)
visNet


saveNetwork(visNet, "manuscript/pictures/visNetwork_celltype_connectome.html")
webshot2::webshot(url = "manuscript/pictures/visNetwork_celltype_connectome.html",
                  file = "manuscript/pictures/visNetwork_celltype_connectome.png",
                  vwidth = 800, vheight = 800, #define the size of the browser window
                  cliprect = c(0, 0, 800, 800), zoom = 5, delay = 1)
