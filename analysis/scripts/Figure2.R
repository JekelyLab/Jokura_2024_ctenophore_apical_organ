# Code to generate Figure 2 of the Jokura et al 2024 Ctenophore apical organ connectome paper

# source packages and functions ------------------------------------------------
source("analysis/scripts/packages_and_functions.R")

# load cell types from catmaid -------------------------------------------------

#read_smooth_cell <- function(annotation){
#  nlapply(read.cells.catmaid(annotation, pid = 35),
#          function(x)
#            smooth_cell(x, sigma = 1000))
#}

#balancer <- read_smooth_cell("celltype:balancer")
#bridge <- read_smooth_cell("celltype:bridge")
#bristle <- read_smooth_cell("celltype:bristle")
#dome <- read_smooth_cell("celltype:dome")
#groove <- read_smooth_cell("celltype:groove")
#intramulticilia <- read_smooth_cell("celltype:intra-multi-ciliated")
#lamellate <- read_smooth_cell("celltype:lamellate")
#lithocyte <- read_smooth_cell("celltype:lithocyte")
#cell <- read_smooth_cell("celltype:neuron")
#plumose <- read_smooth_cell("celltype:plumose")
#dense_vesicle <- read_smooth_cell("celltype:dense_vesicle")
#monocilia <- read_smooth_cell("celltype:monociliated")
#bicilia <- read_smooth_cell("celltype:biciliated")
#non_cilia <- read_smooth_cell("celltype:nonciliated")



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



syncytial_Q1Q2 <-
  nlapply(read.neurons.catmaid(skids_by_2annotations("syncytial neuron", "Q1Q2"),
                               pid = 35),
          function(x)
            smooth_neuron(x, sigma = 1000))
syncytial_Q3Q4 <-
  nlapply(read.neurons.catmaid(skids_by_2annotations("syncytial neuron", "Q3Q4"),
                               pid = 35),
          function(x)
            smooth_neuron(x, sigma = 1000))
syncytial_Q1Q2Q3Q4 <-
  nlapply(read.neurons.catmaid(skids_by_2annotations("syncytial neuron", "Q1Q2Q3Q4"),
                               pid = 35),
          function(x)
            smooth_neuron(x, sigma = 1000))




# circuit analysis and plotting -------------------------------------------

# get connectivity between mech cell clusters and their output clusters
cell_groups <- list(balancer_Q1, balancer_Q2, balancer_Q3, balancer_Q4, 
                    syncytial_Q1Q2, syncytial_Q3Q4, syncytial_Q1Q2Q3Q4)

N_cell_groups <- length(cell_groups)
N_cell_groups


cell_group_attr <- data.frame(
  cell_group_names = c(
    "balancer_Q1",
    "balancer_Q2",
    "balancer_Q3",
    "balancer_Q4",
    "syncytial_Q1Q2",
    "syncytial_Q3Q4",
    "syncytial_Q1Q2Q3Q4"
  ),
  type = c(
    "balancer",
    "balancer",
    "balancer",
    "balancer",
    "neuron",
    "neuron",
    "neuron"
  ),
  level = c("3", "1", "1", "3", "2", "2", "2")
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
  labs(x = "postsynaptic cell groups", y = "presynaptic cell groups", title = " ") +
  scale_x_discrete(limits = cell_group_attr$cell_group_names) +
  scale_y_discrete(limits = cell_group_attr$cell_group_names) +
  scale_fill_gradientn(colours = c("white", "#0072B2")) +
  geom_text(aes(label = synapses, size = synapses / (synapses + 0.1))) +
  scale_radius(range = c(0, 2)) +
  guides(size = "none")

  

# Saving R ggplot with R ggsave Function
ggsave(
  "pictures/mech_girdle_chaeMech_syn_matrix.png",
  width = 1700,
  height = 1300,
  limitsize = TRUE,
  units = c("px")
)























# visNetwork plotting -----------------------------------------------------

## convert to VisNetwork-list
Conn_graph.visn <- toVisNetworkData(Conn_graph)
## copy column "weight" to new column "value" in list "edges"
Conn_graph.visn$edges$value <- Conn_graph.visn$edges$weight
Conn_graph.visn$nodes$value <- degree
Conn_graph.visn$nodes$group <- cell_group_attr$type










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
graph.visn$nodes$color <- Okabe_Ito[1:9]

#hierarchical layout - define level of nodes
graph.visn$nodes$level <- c(1, 2, 3, 4, 1, 2, 3, 4, 1)

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






visNetwork(graph.visn$nodes, graph.visn$edges) %>%
  visIgraphLayout(
    layout = "layout_nicely", physics = TRUE, 
    randomSeed = 42, type="square"
  ) %>%
#  visHierarchicalLayout(
#    levelSeparation=250, 
#    nodeSpacing=10,
#    direction='LR',
#    sortMethod='hubsize',
#    shakeTowards='roots'
#  ) %>%
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
    font=list(color='black', size=1),
    scaling = list(label=list(enabled=TRUE, min=10, max=12)),
    level= graph.visn$nodes$level
  )



results <- sapply(annot_to_search, function(ann1) {
  sapply(annot_to_search, function(ann2) {
    length(skids_by_2annotations(ann1, ann2))
  })
})

# add row and col names
colnames(results) <- annot_to_search
rownames(results) <- annot_to_search
results

# convert to tibble
results.tb <- results %>%
  as_tibble(rownames = "annotation1") %>%
  pivot_longer(
    -annotation1,
    names_to = "annotation2",
    values_to = "number"
  )

results.tb

write.table(results.tb, "source_data/Figure3_fig_suppl3_source_data1.txt", sep = "\t")
results.tb <- read.table("source_data/Figure3_fig_suppl3_source_data1.txt", sep = "\t")


results.tb %>%
  ggplot(aes(annotation2, annotation1)) +
  geom_raster(aes(fill = sqrt(number))) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
    axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 10),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15)
  ) +
  labs(x = "annotation 2", y = "annotation 1", title = " ") +
  scale_x_discrete(limits = unlist(annot_to_search)) +
  scale_y_discrete(limits = rev(unlist(annot_to_search))) +
  scale_fill_gradientn(colours = c("white", "#0072B2")) +
  geom_text(aes(label = number, size = number / (number + 0.1))) +
  scale_radius(range = c(0, 2)) +
  guides(size = "none")
