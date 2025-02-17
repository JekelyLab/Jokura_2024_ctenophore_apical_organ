
#Synaptic connection from SSN to balancer

# source packages and functions -----------------
source("analysis/scripts/packages_and_functions.R")

# load cell type ----------------------------------

balancer <- read_smooth_neuron("celltype:balancer")
bridge <- read_smooth_neuron("celltype:bridge")
bristle <- read_smooth_neuron("celltype:bristle")
Cgroove_sag <- read_smooth_neuron("celltype:Cgroove-sag")
Cgroove_tag <- read_smooth_neuron("celltype:Cgroove-tag")
dense_vesicle <- read_smooth_neuron("celltype:dense_vesicle")
dome <- read_smooth_neuron("celltype:dome")
epithelial_floor <- read_smooth_neuron("celltype:epithelial_floor")
intra_multi_ciliated <- read_smooth_neuron("celltype:intra-multi-ciliated")
lamellate <- read_smooth_neuron("celltype:lamellate")
lithocyte <- read_smooth_neuron("celltype:lithocyte")
plumose <- read_smooth_neuron("celltype:plumose")
SSN <- read_smooth_neuron("celltype:SSN")

monociliated <- read_smooth_neuron("celltype:monociliated")
biciliated <- read_smooth_neuron("celltype:biciliated")
multiciliated <- read_smooth_neuron("celltype:multiciliated")
nonciliated <- read_smooth_neuron("celltype:nonciliated")

Q1 <- read_smooth_neuron("Q1")
Q2 <- read_smooth_neuron("Q2")
Q3 <- read_smooth_neuron("Q3")
Q4 <- read_smooth_neuron("Q4")

all_celltypes <- list(balancer,
                      bridge,
                      bristle,
                      Cgroove_sag,
                      Cgroove_tag,
                      dense_vesicle,
                      dome,
                      epithelial_floor,
                      intra_multi_ciliated,
                      lamellate,
                      lithocyte,
                      plumose,
                      SSN,
                      monociliated,
                      biciliated,
                      multiciliated,
                      nonciliated)


# retrieve connectors ----------------

conn_syn <- connectors(SSN)
presyn_syn <- subset(conn_syn, prepost == 0)
postsyn_syn <- subset(conn_syn, prepost == 1)

# plot syn ----------------------
plot_background()

plot3d(
  SSN, soma = FALSE, 
  color = c("#1f77b4", "#ff7f0e", "#2ca02c"), 
  alpha = 0.5, lwd = c(2,3,3)
)


for (i in 1:length(all_celltypes)) {
  print(i)
  plot3d(
    all_celltypes[[i]], soma = TRUE, lwd = 1, add = TRUE, 
    alpha = 0.05, col = Okabe_Ito[8]
  )
}



#plot3d(
#  outline, WithNodes = F,
#  add=T, alpha=0.07, col="#E2E2E2"
#) 


# plot presynapses
plot3d(
  presyn_syn$x, 
  presyn_syn$y, 
  presyn_syn$z, 
  size = 7, alpha = 1, col = "red", 
  add = T
)

# plot postsynapses
plot3d(
  postsyn_syn$x, 
  postsyn_syn$y, 
  postsyn_syn$z, 
  size = 7, alpha = 1, col = "blue", 
  add = T
)

anterior()
rgl.snapshot("manuscript/pictures/syn_prepost_ant.png")

left()
rgl.snapshot("manuscript/pictures/syn_prepost_left.png")

close3d()


# plot syn ----------------------
plot_background()

plot3d(
  Q1, soma = TRUE, 
  color = Okabe_Ito[1], 
  alpha = 0.7, lwd = c(1,2,3,4)
)
plot3d(
  Q2, soma = TRUE, 
  color = Okabe_Ito[2], 
  alpha = 0.7, lwd = c(1,2,3,4)
)
plot3d(
  Q3, soma = TRUE, 
  color = Okabe_Ito[3], 
  alpha = 0.7, lwd = c(1,2,3,4)
)
plot3d(
  Q4, soma = TRUE, 
  color = Okabe_Ito[4], 
  alpha = 0.7, lwd = c(1,2,3,4)
)



rgl.snapshot("manuscript/pictures/quadrants.png")

close3d()

# assemble figure ---------------------------------------------------------

panel_syn_oral <- ggdraw() + draw_image(readPNG("manuscript/pictures/syn_prepost_ant.png")) +
  draw_label("oral view", x = 0.75, y = 0.95, size = 10) +
  draw_label("nerve net", x = 0.1, y = 0.95, size = 10, hjust = 0)
panel_syn_side <- ggdraw() + draw_image(readPNG("manuscript/pictures/syn_prepost_left.png")) +
  draw_label("side view", x = 0.75, y = 0.95, size = 10)


layout <- "
AB"

Figure_nerve_net <-
  panel_syn_oral + panel_syn_side +
  plot_layout(design = layout) +
  patchwork::plot_layout(design = layout, 
                         heights = c(1,1,1,1,1), 
                         widths = c(1,1,1)) +
  patchwork::plot_annotation(tag_levels = c("A")) &
  ggplot2::theme(plot.tag = element_text(size = 12, face='plain', color='black'))


ggsave(
  "manuscript/figures/Figure_nerve_net.png", 
  limitsize = FALSE, units = c("px"), 
  Figure_nerve_net, width = 1600, height = 800, bg = "white"
)

ggsave(
  "manuscript/figures/Figure_nerve_net.pdf", limitsize = FALSE, 
  units = c("px"), Figure_nerve_net, width = 1600, height = 800
)






























#network analysis and plotting

source("analysis/scripts/packages_and_functions.R")

# load cell type ---------------------------------------------------------------

balancer <- read_smooth_neuron("celltype:balancer")
bridge <- read_smooth_neuron("celltype:bridge")
bristle <- read_smooth_neuron("celltype:bristle")
Cgroove_sag <- read_smooth_neuron("celltype:Cgroove-sag")
Cgroove_tag <- read_smooth_neuron("celltype:Cgroove-tag")
dense_vesicle <- read_smooth_neuron("celltype:dense_vesicle")
dome <- read_smooth_neuron("celltype:dome")
epithelial_floor <- read_smooth_neuron("celltype:epithelial_floor")
intra_multi_ciliated <- read_smooth_neuron("celltype:intra-multi-ciliated")
lamellate <- read_smooth_neuron("celltype:lamellate")
lithocyte <- read_smooth_neuron("celltype:lithocyte")
plumose <- read_smooth_neuron("celltype:plumose")
SSN <- read_smooth_neuron("celltype:SSN")

monociliated <- read_smooth_neuron("celltype:monociliated")
biciliated <- read_smooth_neuron("celltype:biciliated")
multiciliated <- read_smooth_neuron("celltype:multiciliated")
nonciliated <- read_smooth_neuron("celltype:nonciliated")


# cell types
all_celltypes <- list(balancer,
                      bridge,
                      bristle,
                      Cgroove_sag,
                      Cgroove_tag,
                      dense_vesicle,
                      dome,
                      epithelial_floor,
                      intra_multi_ciliated,
                      lamellate,
                      lithocyte,
                      plumose,
                      SSN,
                      monociliated,
                      biciliated,
                      multiciliated,
                      nonciliated)

all_celltypes_names <- list("balancer",
                      "bridge",
                      "bristle",
                      "Cgroove_sag",
                      "Cgroove_tag",
                      "dense_vesicle",
                      "dome",
                      "epithelial_floor",
                      "intra_multi_ciliated",
                      "lamellate",
                      "lithocyte",
                      "plumose",
                      "SSN",
                      "monociliated",
                      "biciliated",
                      "multiciliated",
                      "nonciliated")







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

#AO_graph |> as_tbl_graph() 

# calculate node weighted degree -------------
degree=degree(
  AO_graph, v = V(AO_graph), mode = c("all"), 
  loops = TRUE, normalized = FALSE
  )

# retrieve connectors ----------------

mono_conn <- connectors(monociliated)
presyn_mono <- subset(mono_conn, prepost == 0)
postsyn_mono <- subset(mono_conn, prepost == 1)

bridge_conn <- connectors(bridge)
presyn_bridge <- subset(bridge_conn, prepost == 0)
postsyn_bridge <- subset(bridge_conn, prepost == 1)

bal_conn <- connectors(balancer)
presyn_bal <- subset(bal_conn, prepost == 0)
postsyn_bal <- subset(bal_conn, prepost == 1)

dense_conn <- connectors(dense_vesicle)
presyn_dense <- subset(dense_conn, prepost == 0)
postsyn_dense <- subset(dense_conn, prepost == 1)

dense_noncil <- connectors(nonciliated)
presyn_noncil <- subset(dense_noncil, prepost == 0)
postsyn_noncil <- subset(dense_noncil, prepost == 1)

dense_bicil <- connectors(biciliated)
presyn_bicil <- subset(dense_bicil, prepost == 0)
postsyn_bicil <- subset(dense_bicil, prepost == 1)

# plot monociliated syn ----------------------
plot_background()

plot3d(
  syn_neuron, soma = TRUE, color = Okabe_Ito[c(1,5,8)], 
  alpha = 0.3, lwd = c(4,3,3)
)
plot3d(
  monociliated, soma = TRUE, color = Okabe_Ito[8], 
  alpha = 0.05, lwd = 3
)
# plot postsynapses
plot3d(
  postsyn_mono$x, 
  postsyn_mono$y, 
  postsyn_mono$z, 
  size = 10, alpha = 1, col = "red", 
  add = T
  )

anterior()
rgl.snapshot("manuscript/pictures/mono_syn_ant.png")

left()
rgl.snapshot("manuscript/pictures/mono_syn_left.png")

close3d()

# plot bridge -------------------

plot_background()
plot3d(
  syn_neuron, soma = TRUE, color = Okabe_Ito[c(1,5,8)], 
  alpha = 0.3, lwd = c(4,3,3)
)
plot3d(
  bridge, soma = TRUE, color = Okabe_Ito[8], 
  alpha = 0.1, lwd = 3
)
# plot postsynapses
plot3d(
  postsyn_bridge$x, 
  postsyn_bridge$y, 
  postsyn_bridge$z, 
  size = 10, alpha = 1, col = "red", 
  add = T
  )
plot3d(
  presyn_bridge$x, 
  presyn_bridge$y, 
  presyn_bridge$z, 
  size = 10, alpha = 1, col = "cyan", 
  add = T
  )

anterior()
rgl.snapshot("manuscript/pictures/bridge_syn_ant.png")

left()
rgl.snapshot("manuscript/pictures/bridge_syn_left.png")

close3d()


# plot balancer syn -------------------------------------------------------

plot_background()

plot3d(
  syn_neuron, soma = TRUE, color = Okabe_Ito[c(1,5,8)], 
  alpha = 0.3, lwd = c(4,3,3)
)

plot3d(
  balancer, soma = TRUE, color = Okabe_Ito[8], 
  alpha = 0.05, lwd = 3
  )

# plot postsynapses
plot3d(
  postsyn_bal$x, 
  postsyn_bal$y, 
  postsyn_bal$z, 
  size = 10, alpha = 1, col = "red", 
  add = T
  )

anterior()
rgl.snapshot("manuscript/pictures/bal_syn_ant.png")

left()
rgl.snapshot("manuscript/pictures/bal_syn_left.png")

close3d()


# plot dense vesicle syn --------------------------------------------------

plot_background()
plot3d(
  syn_neuron, soma = TRUE, color = Okabe_Ito[c(1,5,8)], 
  alpha = 0.3, lwd = c(4,3,3)
)
plot3d(
  dense_vesicle, soma = TRUE, color = Okabe_Ito[8], 
  alpha = 0.1, lwd = 3
  )
# plot postsynapses
plot3d(
  postsyn_dense$x, 
  postsyn_dense$y, 
  postsyn_dense$z, 
  size = 10, alpha = 1, col = "red", 
  add = T
  )

anterior()
rgl.snapshot("manuscript/pictures/dense_syn_ant.png")

left()
rgl.snapshot("manuscript/pictures/dense_syn_left.png")

close3d()


# plot nonciliated syn ----------------------------------------------------


plot_background()
plot3d(
  syn_neuron, soma = TRUE, color = Okabe_Ito[c(1,5,8)], 
  alpha = 0.3, lwd = c(4,3,3)
)
plot3d(
  nonciliated, soma = TRUE, color = Okabe_Ito[8], 
  alpha = 0.1, lwd = 3
)
# plot postsynapses
plot3d(
  postsyn_noncil$x, 
  postsyn_noncil$y, 
  postsyn_noncil$z, 
  size = 10, alpha = 1, col = "red", 
  add = T
  )

anterior()
rgl.snapshot("manuscript/pictures/noncil_syn_ant.png")

left()
rgl.snapshot("manuscript/pictures/noncil_syn_left.png")

close3d()


# plot biciliated syn -----------------------------------------------------

plot_background()
plot3d(
  syn_neuron, soma = TRUE, color = Okabe_Ito[c(1,5,8)], 
  alpha = 0.3, lwd = c(4,3,3)
)
plot3d(
  biciliated, soma = TRUE, color = Okabe_Ito[8], 
  alpha = 0.1, lwd = 3
)
# plot postsynapses
plot3d(
  postsyn_bicil$x, 
  postsyn_bicil$y, 
  postsyn_bicil$z, 
  size = 10, alpha = 1, col = "red", 
  add = T
  )

anterior()
rgl.snapshot("manuscript/pictures/bicil_syn_ant.png")

left()
rgl.snapshot("manuscript/pictures/bicil_syn_left.png")

close3d()




# plot all cells in network -----------------------------------------------

plot_background()
plot3d(
  syn_neuron, soma = TRUE, color = Okabe_Ito[3], 
  alpha = 0.3, lwd = c(4,3,3)
)
plot3d(
  bridge, soma = TRUE, color = Okabe_Ito[4], 
  alpha = 0.9, lwd = 5
)
plot3d(
  monociliated, soma = TRUE, color = Okabe_Ito[5], 
  alpha = 0.4, lwd = 3
)
plot3d(
  biciliated, soma = TRUE, color = Okabe_Ito[6], 
  alpha = 0.6, lwd = 3
)
plot3d(
  balancer, soma = TRUE, color = Okabe_Ito[1], 
  alpha = 0.3, lwd = 3
)
plot3d(
  dense_vesicle, soma = TRUE, color = Okabe_Ito[8], 
  alpha = 0.6, lwd = 2
)
plot3d(
  nonciliated, soma = TRUE, color = Okabe_Ito[7], 
  alpha = 0.2, lwd = 3
)
plot3d(
  intra_multiciliated, soma = TRUE, color = Okabe_Ito[2], 
  alpha = 0.8, lwd = 6
)


anterior()
rgl.snapshot("manuscript/pictures/celltypes_in_network_ant.png")
left()
rgl.snapshot("manuscript/pictures/celltypes_in_network_left.png")
close3d()

# use visNetwork to plot the network --------------------------------------

write_rds(AO_graph, "manuscript/source_data/celltype_synapse_matrix.rds")
AO_graph <- read_rds("manuscript/source_data/celltype_synapse_matrix.rds")

## convert to VisNetwork-list
AO_graph.visn <- toVisNetworkData(AO_graph)

#filter low-weight edges
AO_graph.visn$edges$weight
#AO_graph.visn$edges$weight[AO_graph.visn$edges$weight < 3]  <- 0
AO_graph.visn$edges$weight <- sqrt(AO_graph.visn$edges$weight)
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
    ) %>%
  visOptions(highlightNearest = TRUE, width = 1800, height = 1800)
visNet


saveNetwork(visNet, "manuscript/pictures/visNetwork_celltype_connectome.html")
webshot2::webshot(url = "manuscript/pictures/visNetwork_celltype_connectome.html",
                  file = "manuscript/pictures/visNetwork_celltype_connectome.png",
                  vwidth = 1800, vheight = 1800, #define the size of the browser window
                  cliprect = c(390, 250, 880, 1000), zoom = 1, delay = 1)


# assemble figure ---------------------------------------------------------

panel_mono_oral <- ggdraw() + draw_image(readPNG("manuscript/pictures/mono_syn_ant.png")) +
  draw_label("oral view", x = 0.75, y = 0.95, size = 10) +
  draw_label("monociliated", x = 0.1, y = 0.95, size = 10, hjust = 0)
panel_mono_side <- ggdraw() + draw_image(readPNG("manuscript/pictures/mono_syn_left.png")) +
  draw_label("side view", x = 0.75, y = 0.95, size = 10)

panel_bridge_oral <- ggdraw() + draw_image(readPNG("manuscript/pictures/bridge_syn_ant.png")) +
  draw_label("oral view", x = 0.75, y = 0.95, size = 10) +
  draw_label("bridge", x = 0.1, y = 0.95, size = 10, hjust = 0)
panel_bridge_side <- ggdraw() + draw_image(readPNG("manuscript/pictures/bridge_syn_left.png")) +
  draw_label("side view", x = 0.75, y = 0.95, size = 10)

panel_bal_oral <- ggdraw() + draw_image(readPNG("manuscript/pictures/bal_syn_ant.png")) +
  draw_label("oral view", x = 0.75, y = 0.95, size = 10) +
  draw_label("balancer", x = 0.1, y = 0.95, size = 10, hjust = 0)
panel_bal_side <- ggdraw() + draw_image(readPNG("manuscript/pictures/bal_syn_left.png")) +
  draw_label("side view", x = 0.75, y = 0.95, size = 10)

panel_noncil_oral <- ggdraw() + draw_image(readPNG("manuscript/pictures/noncil_syn_ant.png")) +
  draw_label("oral view", x = 0.75, y = 0.95, size = 10) +
  draw_label("nonciliated", x = 0.1, y = 0.95, size = 10, hjust = 0)
panel_noncil_side <- ggdraw() + draw_image(readPNG("manuscript/pictures/noncil_syn_left.png")) +
  draw_label("side view", x = 0.75, y = 0.95, size = 10)

panel_dense_oral <- ggdraw() + draw_image(readPNG("manuscript/pictures/dense_syn_ant.png")) +
  draw_label("oral view", x = 0.75, y = 0.95, size = 10) +
  draw_label("dense vesicle", x = 0.1, y = 0.95, size = 10, hjust = 0)
panel_dense_side <- ggdraw() + draw_image(readPNG("manuscript/pictures/dense_syn_left.png")) +
  draw_label("side view", x = 0.75, y = 0.95, size = 10)

panel_bicil_oral <- ggdraw() + draw_image(readPNG("manuscript/pictures/bicil_syn_ant.png")) +
  draw_label("oral view", x = 0.75, y = 0.95, size = 10) +
  draw_label("biciliated", x = 0.1, y = 0.95, size = 10, hjust = 0)
panel_bicil_side <- ggdraw() + draw_image(readPNG("manuscript/pictures/bicil_syn_left.png")) +
  draw_label("side view", x = 0.75, y = 0.95, size = 10)

panel_all_oral <- ggdraw() + draw_image(readPNG("manuscript/pictures/celltypes_in_network_ant.png")) +
  draw_label("oral view", x = 0.75, y = 0.95, size = 10) +
  draw_label("all cells", x = 0.1, y = 0.95, size = 10, hjust = 0)
panel_all_side<- ggdraw() + draw_image(readPNG("manuscript/pictures/celltypes_in_network_left.png")) +
  draw_label("side view", x = 0.75, y = 0.95, size = 10)

panel_conn <- ggdraw() + draw_image(readPNG("manuscript/pictures/visNetwork_celltype_connectome.png"))

layout <- "
ABCDEF
GGHIJK
GGLMNO"

Figure_network <-
  panel_all_oral + panel_all_side +
  panel_mono_oral + panel_mono_side + 
  panel_bridge_oral + panel_bridge_side +
  panel_conn +
  panel_bal_oral + panel_bal_side + panel_noncil_oral + panel_noncil_side +
  panel_dense_oral + panel_dense_side + panel_bicil_oral + panel_bicil_side +
  plot_layout(design = layout) +
  patchwork::plot_layout(design = layout, 
                         heights = c(1,1,1,1,1), 
                         widths = c(1,1,1)) +
  patchwork::plot_annotation(tag_levels = list(c(
    "A", "B", "C", "", "D", "", "E", "F", "", "G", "", "H", "", "I"))) &
  ggplot2::theme(plot.tag = element_text(size = 12, face='plain', color='black'))


ggsave("manuscript/figures/Figure_network.png", limitsize = FALSE, 
       units = c("px"), Figure_network, width = 3800, height = 2000, bg = "white")

ggsave("manuscript/figures/Figure_network.pdf", limitsize = FALSE, 
       units = c("px"), Figure_network, width = 3800, height = 2000)


