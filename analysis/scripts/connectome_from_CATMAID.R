# code to generate the ctenophore apical organ connectome graph based on the CATMAID database
# Gaspar Jekely 2023-2025

# source packages and functions ------------------------------------------------
source("analysis/scripts/packages_and_functions.R")

# get all synapses from CATMAID
all_syn_connectors <- catmaid_fetch(
  path = "35/connectors/",
  body = list(
    relation_type = "presynaptic_to",
    relation_type = "postsynaptic_to",
    with_partners = "true"
  )
)
length(all_syn_connectors$connectors)

connector_ids <- unlist(unname(sapply(all_syn_connectors$connectors, "[[", 1)))
connectivity <- catmaid_get_connectors(connector_ids, pid = 35)
length(connector_ids)
dim(connectivity)[1]
connectivity

#create graph
edges <- tibble(
  pre = as.character(connectivity |>
            select(pre) |>
            pull()),
  post = as.character(connectivity |>
             select(post) |>
             pull()
           )
)
  
conn.tb <- tbl_graph(edges = edges, directed = TRUE)
conn.tb

# retrieve node names, simplify and filter the graph --------------------------------------------------------

# add weighted degree to nodes
conn.tb.str <- conn.tb %>%
  activate(nodes) %>%
  mutate(strength = centrality_degree(
    weights = NULL,
    mode = "all",
    loops = TRUE,
    normalized = FALSE
  ))

# test if graph is simple (no multi-edges)
is_simple(as.igraph(conn.tb.str))

# add edge weight of 1
conn.tb.str <- set_edge_attr(conn.tb.str, "weight", value = 1)
conn.tb.str

# merge edges by converting to simple and summing edge weights
conn.tb.str.sum <- conn.tb.str %>%
  activate(edges) %>%
  convert(to_simple) %>%
  mutate(weight = map_dbl(.orig_data, ~ sum(.x$weight))) %>%
  select(-.orig_data, -.tidygraph_edge_index) %>%
  activate(nodes) %>%
  select(-.tidygraph_node_index)

conn.tb.str.sum
is_simple(as.igraph(conn.tb.str.sum))
components(conn.tb.str.sum)

# filter by node strength
conn.tb.str.sum.filt <- conn.tb.str.sum %>%
  activate(nodes) %>%
  filter(strength > 0)

conn.tb.str.sum.filt

# check connected components  ---------------------------

# check connected components, the connectome should only have one component
cl <- components(conn.tb.str.sum.filt)
# size of the largest subnetwork
length(which(cl$membership == 1))

# generate subgraph
conn.tb.str.sum.filt.cl <- subgraph(conn.tb.str.sum.filt, which(cl$membership == 1))

# check if graph is directed
igraph::is_directed(conn.tb.str.sum.filt.cl)

# shorten name, convert to tbl
conn.tb <- conn.tb.str.sum.filt.cl %>% as_tbl_graph()

conn.tb

# check names, search for annotations --------------------------------------------------

# list skids
skids <- conn.tb %>%
  activate(nodes) %>%
  select(name) %>%
  pull()

# get neuron names
names <- catmaid_get_neuronnames(skids, pid = 35)

conn.tb.skids.names <- conn.tb %>%
  mutate(names) %>%
  rename(skids = name)

# list neuron names
cell_names <- conn.tb.skids.names %>%
  activate(nodes) %>%
  select(names) %>%
  pull() %>%
  unname()

# check duplicated names (there should be none)
cell_names[duplicated(cell_names)]
length(cell_names)
cell_names

# centrality measures --------------------------------------

conn.tb.skids.names <- conn.tb.skids.names %>%
  activate(nodes) %>%
  mutate("betweenness" = centrality_betweenness(
    weights = E(conn.tb.str.sum.filt.cl)$weight,
    directed = TRUE
  )) %>%
  mutate("authority" = centrality_authority(weights = E(conn.tb.str.sum.filt.cl)$weight)) %>%
  mutate("pagerank" = centrality_pagerank(
    weights = E(conn.tb.str.sum.filt.cl)$weight,
    directed = TRUE
  )) %>%
  mutate("closeness" = centrality_closeness(
    weights = E(conn.tb.str.sum.filt.cl)$weight,
    mode = "all"
  )) %>%
  mutate("eigen" = centrality_eigen(directed = TRUE)) %>%
  mutate("hub" = centrality_hub(weights = E(conn.tb.str.sum.filt.cl)$weight)) %>%
  mutate("node_eccentricity_out" = node_eccentricity(mode = "out")) %>%
  mutate("node_eccentricity_in" = node_eccentricity(mode = "in")) %>%
  mutate("node_is_sink" = node_is_sink()) %>%
  mutate("node_is_source" = node_is_source()) %>%
  mutate("node_is_cut" = node_is_cut()) %>%
  mutate("local_triangles" = local_triangles())

conn.tb.skids.names
# save conn.tb
saveRDS(conn.tb.skids.names, "manuscript/source_data/connectome_graph_tibble.rds")
conn.tb.skids.names <- readRDS("manuscript/source_data/connectome_graph_tibble.rds")

# VisNetwork conversion ---------------------------------------------------

# convert to visNetwork graph
conn_graph.visn <- toVisNetworkData(conn.tb.skids.names)

## copy column "weight" to new column "value" in list "edges"
conn_graph.visn$edges$value <- conn_graph.visn$edges$weight
conn_graph.visn$nodes$value <- conn_graph.visn$nodes$strength

