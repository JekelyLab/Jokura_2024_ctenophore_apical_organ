# Code to generate Figure 2 of the Jokura et al 2024 Ctenophore apical organ connectome paper

# source packages and functions ------------------------------------------------
source("analysis/scripts/packages_and_functions.R")

# load cell types from catmaid -------------------------------------------------

read_smooth_cell <- function(annotation){
  nlapply(read.neurons.catmaid(annotation, pid = 35),
          function(x)
            smooth_neuron(x, sigma = 1000))
}


bridge_Q1Q2 <-
  nlapply(read.neurons.catmaid(skids_by_2annotations("celltype:bridge", "Q1Q2"),
                               pid = 35),
          function(x)
            smooth_neuron(x, sigma = 1000))

bridge_Q3Q4 <-
  nlapply(read.neurons.catmaid(skids_by_2annotations("celltype:bridge", "Q3Q4"),
                               pid = 35),
          function(x)
            smooth_neuron(x, sigma = 1000))


balancer_Q1 <-
  nlapply(read.neurons.catmaid(skids_by_2annotations("celltype:balancer", "Q2"),
                               pid = 35),
          function(x)
            smooth_neuron(x, sigma = 1000))
balancer_Q2 <-
  nlapply(read.neurons.catmaid(skids_by_2annotations("celltype:balancer", "Q2"),
                               pid = 35),
          function(x)
            smooth_neuron(x, sigma = 1000))
balancer_Q3 <-
  nlapply(read.neurons.catmaid(skids_by_2annotations("celltype:balancer", "Q3"),
                               pid = 35),
          function(x)
            smooth_neuron(x, sigma = 1000))
balancer_Q4 <-
  nlapply(read.neurons.catmaid(skids_by_2annotations("celltype:balancer", "Q4"),
                               pid = 35),
          function(x)
            smooth_neuron(x, sigma = 1000))


SSN_Q1Q2 <-
  nlapply(read.neurons.catmaid(skids_by_2annotations("celltype:SSN", "Q1Q2"),
                               pid = 35),
          function(x)
            smooth_neuron(x, sigma = 1000))
SSN_Q3Q4 <-
  nlapply(read.neurons.catmaid(skids_by_2annotations("celltype:SSN", "Q3Q4"),
                               pid = 35),
          function(x)
            smooth_neuron(x, sigma = 1000))
SSN_Q1Q2Q3Q4 <-
  nlapply(read.neurons.catmaid(skids_by_2annotations("celltype:SSN", "Q1Q2Q3Q4"),
                               pid = 35),
          function(x)
            smooth_neuron(x, sigma = 1000))




# circuit analysis and plotting -------------------------------------------

# get connectivity between mech cell clusters and their output clusters
cell_groups <- list(bridge_Q1Q2, bridge_Q3Q4, 
                    balancer_Q1, balancer_Q2, balancer_Q3, balancer_Q4, 
                    SSN_Q1Q2, SSN_Q3Q4, SSN_Q1Q2Q3Q4)


N_cell_groups <- length(cell_groups)
N_cell_groups


cell_group_attr <- data.frame(
  cell_group_names = c(
    "bridge_Q1Q2",
    "bridge_Q3Q4",
    "balancer_Q1",
    "balancer_Q2",
    "balancer_Q3",
    "balancer_Q4",
    "SSN_Q1Q2",
    "SSN_Q3Q4",
    "SSN_Q1Q2Q3Q4"
    ),
  type = c(
    "bridge",
    "bridge",
    "balancer",
    "balancer",
    "balancer",
    "balancer",
    "neuron",
    "neuron",
    "neuron"
  ),
  level = c("2", "2", "3", "1", "1", "3", "2", "2", "2")
)

dim(cell_group_attr)



# iterate through cell group neuron lists and get connectivity for all against all
{
  # define empty synapse list with the right dimensions
  synapse_list <- vector("list", N_cell_groups * N_cell_groups)
  for (i in 1:N_cell_groups) {
    for (j in 1:N_cell_groups) {
      # get connectors between two cell groups
      presyn_skids <- attr(cell_groups[i][[1]], "df")$skid
      postsyn_skids <- attr(cell_groups[j][[1]], "df")$skid
      connectivity <- catmaid_get_connectors_between(
        pre = presyn_skids,
        post = postsyn_skids, pid = 35
      )
      # check the number of synapses from group1 -> group2
      N_synapses <- dim(connectivity)[1]
      # change "NULL" to 0
      if (is.null(N_synapses)) {
        N_synapses <- 0
      }
      print((i * N_cell_groups - N_cell_groups) + j)
      print(N_synapses)
      # add value to synapse list
      synapse_list[[(i * N_cell_groups - N_cell_groups) + j]] <- N_synapses
    }
  }
}




# convert synapse list into a matrix of appropriate dimensions
synapse_matrix <- matrix(unlist(synapse_list), byrow = TRUE, nrow = N_cell_groups)
rownames(synapse_matrix) <- cell_group_attr$cell_group_names
colnames(synapse_matrix) <- cell_group_attr$cell_group_names
synapse_matrix



# plot with ggplot---------------------------------------------------------

as.data.frame((synapse_matrix)) %>%
  rownames_to_column(var = "presyn_cell_group") %>%
  pivot_longer(-presyn_cell_group,
               names_to = "postsyn_cell_group",
               values_to = "synapses") %>%
  group_by(postsyn_cell_group) %>%
  mutate(synapse_fraction = synapses / sum(synapses, na.rm = TRUE)) %>%
  ggplot(aes(x = postsyn_cell_group, y = presyn_cell_group)) +
  geom_raster(aes(fill = sqrt(synapses))) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 10
    ),
    axis.text.y = element_text(
      angle = 0,
      hjust = 1,
      vjust = 0.5,
      size = 10
    ),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15)
  ) +
  labs(x = "postsynaptic cell groups", y = "presynaptic cell groups", title = "") +
  scale_x_discrete(limits = cell_group_attr$cell_group_names) +
  scale_y_discrete(limits = cell_group_attr$cell_group_names) +
  scale_fill_gradientn(colours = c("white", "#0072B2")) +
  geom_text(aes(label = synapses, size = synapses / (synapses + 0.5))) +
  scale_radius(range = c(0, 2)) +
  guides(size = "none")

  

# Saving R ggplot with R ggsave Function
ggsave(
  "manuscript/pictures/mech_girdle_chaeMech_syn_matrix.png",
  width = 1800,
  height = 1200,
  limitsize = TRUE,
  units = c("px")
)
























































# network---------------------------------------------------------

# cell types
cell_groups <- list(bridge_Q1Q2, bridge_Q3Q4, 
                    balancer_Q1, balancer_Q2, balancer_Q3, balancer_Q4, 
                    SSN_Q1Q2, SSN_Q3Q4, SSN_Q1Q2Q3Q4)

cell_groups_names <- list("bridge_Q1Q2", "bridge_Q3Q4", 
                          "balancer_Q1", "balancer_Q2", "balancer_Q3", "balancer_Q4", 
                          "SSN_Q1Q2", "SSN_Q3Q4", "SSN_Q1Q2Q3Q4"
)

# iterate through cell group neuron lists and get connectivity
# define empty synapse list with the right dimensions
synapse_list <- c()

for (i in 1:length(cell_groups)) {
  for (j in 1:length(cell_groups)) {
    # get connectors between two cell groups
    presyn_skids <- attr(cell_groups[i][[1]], "df")$skid
    postsyn_skids <- attr(cell_groups[j][[1]], "df")$skid
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
  nrow = length(cell_groups)
)

rownames(synapse_matrix) <- as.character(cell_groups_names)
colnames(synapse_matrix) <- as.character(cell_groups_names)

# with the make_graph function of igraph we turn it into a graph (input is the list of edge pairs)
cell_groups_graph <- graph_from_adjacency_matrix(
  synapse_matrix,
  mode = c("directed"),
  weighted = TRUE, diag = TRUE
)

#cell_groups_graph |> as_tbl_graph() 

# calculate node weighted degree -------------
degree=degree(
  cell_groups_graph, v = V(cell_groups_graph), mode = c("all"), 
  loops = TRUE, normalized = FALSE
)

# use visNetwork to plot the network --------------------------------------

write_rds(cell_groups_graph, "manuscript/source_data/cell_groups_synapse_matrix.rds")
cell_groups_graph <- read_rds("manuscript/source_data/cell_groups_synapse_matrix.rds")

## convert to VisNetwork-list
cell_groups_graph.visn <- toVisNetworkData(cell_groups_graph)

#filter low-weight edges
cell_groups_graph.visn$edges$weight
#cell_groups_graph.visn$edges$weight[cell_groups_graph.visn$edges$weight < 3]  <- 0
cell_groups_graph.visn$edges$weight <- sqrt(cell_groups_graph.visn$edges$weight)
cell_groups_graph.visn$edges$weight

## copy column "weight" to new column "value" in list "edges"
cell_groups_graph.visn$edges$value <- cell_groups_graph.visn$edges$weight
cell_groups_graph.visn$nodes$value <- degree

#define node color
cell_groups_graph.visn$nodes$color <- Okabe_Ito[1:9]

#hierarchical layout - define level of nodes
cell_groups_graph.visn$nodes$level <- c(2, 2, 3, 1, 1, 3, 2, 2, 2)
#bridge_Q1Q2, bridge_Q3Q4,balancer_Q1, balancer_Q2, balancer_Q3, balancer_Q4,SSN_Q1Q2, SSN_Q3Q4, SSN_Q1Q2Q3Q4

#hierarchical layout
visNet <- visNetwork(cell_groups_graph.visn$nodes, cell_groups_graph.visn$edges) %>%
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
    color = list(background=cell_groups_graph.visn$nodes$color, border='black'),
    opacity=0.9,
    shape='dot', 
    font=list(color='black', size=44),
    scaling = list(label=list(enabled=TRUE, min=22, max=80)),
    level= cell_groups_graph.visn$nodes$level
  ) %>%
  visOptions(highlightNearest = TRUE, width = 1800, height = 1800)
visNet


saveNetwork(visNet, "manuscript/pictures/visNetwork_bridge_balancer_connectome.html")
webshot2::webshot(url = "manuscript/pictures/visNetwork_bridge_balancer_connectome.html",
                  file = "manuscript/pictures/visNetwork_bridge_balancer_connectome.png",
                  vwidth = 1800, vheight = 1800, #define the size of the browser window
                  cliprect = c(390, 250, 880, 1000), zoom = 1, delay = 1)



