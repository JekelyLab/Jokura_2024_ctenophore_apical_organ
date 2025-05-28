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

# create graph
edges <- tibble(
  pre = as.character(connectivity |>
    select(pre) |>
    pull()),
  post = as.character(connectivity |>
    select(post) |>
    pull())
)

conn.tb <- tbl_graph(edges = edges, directed = TRUE)
conn.tb

# retrieve node names, simplify and filter the graph --------------------------------------------------------

# add weighted degree to nodes
add_node_wdegree <- function(graph.tb){
  graph.tb|>
    activate(nodes) |>
    mutate(strength = centrality_degree(
      weights = NULL,
      mode = "all",
      loops = TRUE,
      normalized = FALSE
    ))
}
  
conn.tb.str <- add_node_wdegree(conn.tb)

# test if graph is simple (no multi-edges)
is_simple(as.igraph(conn.tb.str))

# add edge weight of 1
conn.tb.str <- set_edge_attr(conn.tb.str, "weight", value = 1)
conn.tb.str

# merge edges by converting to simple and summing edge weights

conn.tb.str.sum <- conn.tb.str |>
  activate(edges) |>
  group_by(from, to) |>
  mutate(weight = sum(weight)) |>
  distinct() |>
  activate(nodes)

# merge edges by converting to simple and summing edge weights
#conn.tb.str.sum <- conn.tb.str %>%
#  activate(edges) %>%
#  convert(to_simple) %>%
#  mutate(weight = map_dbl(.orig_data, ~ sum(.x$weight))) %>%
#  select(-.orig_data, -.tidygraph_edge_index) %>%
#  activate(nodes) %>%
#  select(-.tidygraph_node_index)

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
conn.tb.str.sum.filt.cl <- subgraph(
  conn.tb.str.sum.filt, 
  which(cl$membership == 1)
  )

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

# get neuron names from CATMAID
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

# truncate names for plot ---------------------

conn.tb.skids.names <- conn.tb.skids.names |>
  mutate(names_short = sub("_.*", "", names)) |>
  mutate(names_short = sub("\\s.*", "", names_short))

# add annotations to nodes -------------------

# annotations to search for
annot_to_search <- c(
  "Q1", "Q2", "Q3",
  "Q4", "Q1Q2", "Q3Q4",
  "Q1Q2Q3Q4"
)

#use purrr:map to get all annotations for all skids
annotations <- map(skids, pid = 35, catmaid_get_annotations_for_skeletons)

#function to match terms to the annotations
find_annotation <- function(annots, term){
  if(sum(unlist(annots) %in% term) == 1){
    return(term)
  }
  else {
    return(NA)
  }
}

#search the annotations with the function
annot.tb <- map_dfc(annot_to_search, ~map_chr(annotations, find_annotation, .x))
annot.tb
#generate list wilag()#generate list with quadra...7#generate list with quadrant annotations
quadrant_of_cell <- annot.tb |>
  mutate(cons = case_when(
    !is.na(...1) ~ ...1,
    !is.na(...2) ~ ...2,
    !is.na(...3) ~ ...3,
    !is.na(...4) ~ ...4,
    !is.na(...5) ~ ...5,
    !is.na(...6) ~ ...6,
    !is.na(...7) ~ ...7)
  ) |>
  select(cons) |>
  pull()


conn.tb.skids.names <- conn.tb.skids.names |>
  mutate(quadrant = quadrant_of_cell)

# color nodes by quadrant, matching in the name
conn.tb.skids.names.col <- conn.tb.skids.names |>
  mutate(color = case_when(
    quadrant == "Q1Q2Q3Q4" ~ "#DDDDDD",
    quadrant == "Q1Q2" ~ "#EEEEEE",
    quadrant == "Q3Q4" ~ "#EEEEEE",
    quadrant == "else" ~ "#AAAAAA",
    quadrant == "Q1" ~ Okabe_Ito[1],
    quadrant == "Q2" ~ Okabe_Ito[2],
    quadrant == "Q3" ~ Okabe_Ito[6],
    quadrant == "Q4" ~ Okabe_Ito[7]
  )) |>
  mutate(shape = case_when(
    quadrant == "Q1Q2Q3Q4" ~ "circle",
    quadrant == "Q1Q2" ~ "circle",
    quadrant == "Q3Q4" ~ "circle",
    quadrant == "else" ~ "dot",
    quadrant == "Q1" ~ "dot",
    quadrant == "Q2" ~ "dot",
    quadrant == "Q3" ~ "dot",
    quadrant == "Q4" ~ "dot"
  ))

conn.tb.skids.names.col |>
  select(color) |>
  pull()

conn.tb.skids.names.col |>
  select(names_short) |>
  pull()

# save conn.tb
saveRDS(conn.tb.skids.names.col, "manuscript/source_data/connectome_graph_tibble.rds")
conn.tb.skids.names.col <- readRDS("manuscript/source_data/connectome_graph_tibble.rds")

# VisNetwork conversion ---------------------------------------------------

# convert to visNetwork graph
conn_graph.visn <- toVisNetworkData(conn.tb.skids.names.col)

## copy column "weight" to new column "value" in list "edges"
conn_graph.visn$edges$value <- conn_graph.visn$edges$weight
conn_graph.visn$nodes$value <- conn_graph.visn$nodes$strength

conn_graph.visn$nodes$label <- conn_graph.visn$nodes$names_short

# hierarchical layout
visNet <- visNetwork(conn_graph.visn$nodes, conn_graph.visn$edges) %>%
  visIgraphLayout(
    layout = "layout_nicely", physics = FALSE
  ) %>%
  visEdges(
    smooth = list(type = "curvedCW", roundness = 0.2),
    scaling = list(min = 2, max = 12),
    color = list(inherit = TRUE, opacity = 0.7),
    font = list(color = "black", size = 22),
    arrows = list(
      to = list(
        enabled = TRUE,
        scaleFactor = 1, type = "arrow"
      )
    )
  ) %>%
  visNodes(
    borderWidth = 0.3,
    color = list(background = conn_graph.visn$nodes$color, border = "black"),
    opacity = 0.9,
    #  shape='circle',
    font = list(color = "black", size = 32),
    scaling = list(label = list(enabled = TRUE, min = 24, max = 50))
  ) %>%
  visOptions(highlightNearest = TRUE, width = 1800, height = 1800)

visNet


saveNetwork(visNet, "manuscript/pictures/cell_level_connectome.html")
webshot2::webshot(
  url = "manuscript/pictures/cell_level_connectome.html",
  file = "manuscript/pictures/cell_level_connectome.png",
  vwidth = 1800, vheight = 1800, # define the size of the browser window
  cliprect = c(350, 400, 1400, 700), zoom = 5, delay = 1
)
