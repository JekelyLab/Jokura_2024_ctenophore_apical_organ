

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