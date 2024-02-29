source("analysis/scripts/packages_and_functions.R")

# get all abutting contacts from CATMAID ------------

contacts <- catmaid_fetch(
  path = "35/connectors/",
  body = list(
    relation_type = "abutting",
    with_partners = "true"
  )
)

# method1 -------------------------------
start.time <- Sys.time()

authtoken <- conn$value$token

has_enter_cil_nocil_tag <- function(treenode_id) {
  # catmaid fetch doesn't work here for some reason?
  tags_of_node <- GET(paste(
    "https://catmaid.ex.ac.uk/35/labels/treenode/", treenode_id, sep = ""
    ), 
    add_headers(Authorization = paste("Token", authtoken))) |>
    content()
  enter_cil <- "enter_cil" %in% tags_of_node
  enter_nocil <- "enter_nocil" %in% tags_of_node
  return(c(enter_cil, enter_nocil))
}

enter_cell_info <- function(connector_id) {
  connector_info <- GET(paste(
     "https://catmaid.ex.ac.uk/35/connectors/", connector_id, sep = ""
  ), 
  add_headers(Authorization = paste("Token", authtoken))) |>
    content()
  partner1_treenode <- connector_info$partners[[1]][[2]]
  partner2_treenode <- connector_info$partners[[2]][[2]]
  partner1_skid <- connector_info$partners[[1]][[4]]
  partner2_skid <- connector_info$partners[[2]][[4]]
  newlist <- c(
    connector_id,
    partner1_skid,
    partner2_skid,
    has_enter_cil_nocil_tag(partner1_treenode),
    has_enter_cil_nocil_tag(partner2_treenode)
  )
  names(newlist) <- c("connector_id",
                      "partner1_skid",
                      "partner2_skid",
                      "partner1_enter_cil",
                      "partner1_enter_nocil",
                      "partner2_enter_cil",
                      "partner2_enter_nocil"
                      )
  return(newlist)
}
          
# get connector ids
connector_ids <- lapply(contacts$connectors, "[[", 1) |> unlist()

connector_prepost_info <- lapply(connector_ids, enter_cell_info)

connector_prepost_info_table <- bind_rows(connector_prepost_info, .id = "fitname")

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

start.time <- Sys.time()

# method 2 -------------------------------------
# get connector ids
connector_ids <- c()
for(i in 1:length(contacts$connectors)){
  connector_ids[i] <- contacts$connectors[[i]][[1]]
  }
connector_ids

partner1_skid <- c()
for(i in 1:length(contacts$connectors)){
  partner1_skid[i] <- contacts$partners[[i]][[1]][[3]]
}
partner1_skid

partner2_skid <- c()
for(i in 1:length(contacts$connectors)){
  partner2_skid[i] <- contacts$partners[[i]][[2]][[3]]
}
partner2_skid

node_id_of_contact_partner1 <- c()
for(i in 1:length(contacts$connectors)){
  node_id_of_contact_partner1[i] <- contacts$partners[[i]][[1]][[2]]
}
node_id_of_contact_partner1

node_id_of_contact_partner2 <- c()
for(i in 1:length(contacts$connectors)){
  node_id_of_contact_partner2[i] <- contacts$partners[[i]][[2]][[2]]
}
node_id_of_contact_partner2


#check for tags enter_cil and enter_nocil
partner1_pre_cil <- c()
partner1_pre_nocil <- c()
for(i in 1:length(contacts$connectors)){
  neuron <- read.neuron.catmaid(partner1_skid[i], pid = 35)
  partner1_pre_cil[i] <- node_id_of_contact_partner1[i] %in% neuron$tags$enter_cil
  partner1_pre_nocil[i] <- node_id_of_contact_partner1[i] %in% neuron$tags$enter_nocil
}
partner1_pre_cil
partner1_pre_nocil

partner2_pre_cil <- c()
partner2_pre_nocil <- c()
for(i in 1:length(contacts$connectors)){
  neuron <- read.neuron.catmaid(partner2_skid[i], pid = 35)
  partner2_pre_cil[i] <- node_id_of_contact_partner2[i] %in% neuron$tags$enter_cil
  partner2_pre_nocil[i] <- node_id_of_contact_partner2[i] %in% neuron$tags$enter_nocil
}
partner2_pre_cil
partner2_pre_nocil

# define connectivity tibble

connector_tb <- tibble(
  connectors = connector_ids,
  skid1 = partner1_skid,
  skid2 = partner2_skid,
  partner1_pre_cil = partner1_pre_cil,
  partner1_pre_nocil = partner1_pre_nocil,
  partner2_pre_cil = partner2_pre_cil,
  partner2_pre_nocil = partner2_pre_nocil
)

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

connector_tb


# derive network from table -----------------


connector_prepost_info_table

#get all nodes (all skids)
nodes1 <- connector_prepost_info_table |>
  pull(partner1_skid)
nodes2 <- connector_prepost_info_table |>
  pull(partner2_skid)
all_nodes <- unique(c(nodes1, nodes2))

#edge lists
edges1_to_2_from <- connector_prepost_info_table |>
  filter(partner1_enter_cil == 1) |>
  pull(partner1_skid) |> 
  as.character()
edges1_to_2_to <- connector_prepost_info_table |>
  filter(partner1_enter_cil == 1) |>
  pull(partner2_skid) |> 
  as.character()
edges2_to_1_from <- connector_prepost_info_table |>
  filter(partner2_enter_cil == 1) |>
  pull(partner2_skid) |> 
  as.character()
edges2_to_1_to <- connector_prepost_info_table |>
  filter(partner2_enter_cil == 1) |>
  pull(partner1_skid) |> 
  as.character()
edges2_to_1_from
edges2_to_1_to

# make graph
cil_graph <- create_empty(0)

#add nodes
cil_graph <- cil_graph |> 
  bind_nodes(as.data.frame(all_nodes))

#add edges
cil_graph <- cil_graph |> bind_edges(data.frame(from = edges1_to_2_from, to = edges1_to_2_to))
cil_graph <- cil_graph |> bind_edges(data.frame(from = edges2_to_1_from, to = edges2_to_1_to))
cil_graph

# derive network method 2 ----------------

connector_tb

#get all nodes (all skids)
nodes1 <- connector_tb |>
  pull(skid1)
nodes2 <- connector_tb |>
  pull(skid2)
all_nodes <- unique(c(nodes1, nodes2))

#edge lists
edges1_to_2_from <- connector_tb |>
  filter(partner1_pre_cil == TRUE) |>
  pull(skid1) |> 
  as.character()
edges1_to_2_from
edges1_to_2_to <- connector_tb |>
  filter(partner1_pre_cil == TRUE) |>
  pull(skid2) |> 
  as.character()
edges1_to_2_to
edges2_to_1_from <- connector_tb |>
  filter(partner2_pre_cil == TRUE) |>
  pull(skid2) |> 
  as.character()
edges2_to_1_from
edges2_to_1_to <- connector_tb |>
  filter(partner2_pre_cil == TRUE) |>
  pull(skid1) |> 
  as.character()
edges2_to_1_to

# make graph
cil_graph2 <- create_empty(0, directed = TRUE)

#add nodes
cil_graph2 <- cil_graph |> 
  bind_nodes(as.data.frame(all_nodes))

#add edges
cil_graph2 <- cil_graph |> bind_edges(data.frame(from = edges1_to_2_from, to = edges1_to_2_to))
cil_graph2 <- cil_graph |> bind_edges(data.frame(from = edges2_to_1_from, to = edges2_to_1_to))
cil_graph2

# nocil graph -----------------------


connector_prepost_info_table

#get all nodes (all skids)
nodes1 <- connector_prepost_info_table |>
  pull(partner1_skid)
nodes2 <- connector_prepost_info_table |>
  pull(partner2_skid)
all_nodes <- unique(c(nodes1, nodes2))

neuron_names <- catmaid_get_neuronnames(all_nodes, pid = 35)

#edge lists
edges1_to_2_from <- connector_prepost_info_table |>
  filter(partner1_enter_nocil == 1) |>
  pull(partner1_skid) |> 
  as.character()
edges1_to_2_to <- connector_prepost_info_table |>
  filter(partner1_enter_nocil == 1) |>
  pull(partner2_skid) |> 
  as.character()
edges2_to_1_from <- connector_prepost_info_table |>
  filter(partner2_enter_nocil == 1) |>
  pull(partner2_skid) |> 
  as.character()
edges2_to_1_to <- connector_prepost_info_table |>
  filter(partner2_enter_nocil == 1) |>
  pull(partner1_skid) |> 
  as.character()
edges2_to_1_from
edges2_to_1_to

# make graph
nocil_graph <- create_empty(0)

#add nodes
nocil_graph <- nocil_graph |> 
  bind_nodes(as.data.frame(all_nodes))

#add edges
nocil_graph <- nocil_graph |> bind_edges(data.frame(from = edges1_to_2_from, to = edges1_to_2_to))
nocil_graph <- nocil_graph |> bind_edges(data.frame(from = edges2_to_1_from, to = edges2_to_1_to))
nocil_graph <- nocil_graph |>
  mutate(name = unlist(neuron_names))
nocil_graph <- nocil_graph |> activate(edges) |>
  mutate(value = 1)

nocil_graph

as.igraph(nocil_graph)
# graph plot --------------------------------------------------------------


## convert to VisNetwork-list
nocil_graph.visn <- toVisNetworkData(nocil_graph)
nocil_graph.visn

#define node color
nocil_graph.visn$nodes$color <- sample(Okabe_Ito[1:8], 102, replace = TRUE)
nocil_graph.visn$nodes$value <- 1

#hierarchical layout
visNetwork(nocil_graph.visn$nodes, nocil_graph.visn$edges) %>%
  visIgraphLayout(
    layout = "layout_nicely", physics = TRUE, 
    randomSeed = 42, type="square"
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
    color = list(background=nocil_graph.visn$nodes$color, border='black'),
    opacity=0.9,
    shape='dot', 
    font=list(color='black', size=44),
    scaling = list(label=list(enabled=TRUE, min=48, max=56))
    )
