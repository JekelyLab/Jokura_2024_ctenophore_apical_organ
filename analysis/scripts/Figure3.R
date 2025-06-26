# Code to generate Figure 3 of the Jokura et al 2024 Ctenophore apical organ connectome paper

# source packages and functions ------------------------------------------------
source("analysis/scripts/packages_and_functions.R")


# load cell type ---------------------------------------------------------------

with_soma <- read_smooth_neuron("with_soma")

balancer <- read_smooth_neuron("celltype:balancer")
bridge <- read_smooth_neuron("celltype:bridge")

balancer_Q1 <- read_smooth_neuron(get_skids_with_annot(pid = 35, c("celltype:balancer", "Q1")))
balancer_Q2 <- read_smooth_neuron(get_skids_with_annot(pid = 35, c("celltype:balancer", "Q2")))
balancer_Q3 <- read_smooth_neuron(get_skids_with_annot(pid = 35, c("celltype:balancer", "Q3")))
balancer_Q4 <- read_smooth_neuron(get_skids_with_annot(pid = 35, c("celltype:balancer", "Q4")))

bridge_Q1Q2 <- read_smooth_neuron(get_skids_with_annot(pid = 35, c("celltype:bridge", "Q1Q2")))
bridge_Q3Q4 <- read_smooth_neuron(get_skids_with_annot(pid = 35, c("celltype:bridge", "Q3Q4")))


balancer_skids <- names(balancer)


# bar graph of outputs from SNNs -----------------------------------------------

SSN_Q1Q2 <- read_smooth_neuron("SSN_Q1Q2")[[1]]
SSN_Q3Q4 <- read_smooth_neuron("SSN_Q3Q4")[[1]]
SSN_Q1Q2Q3Q4 <- read_smooth_neuron("SSN_Q1Q2Q3Q4")[[1]]

SSN_Q1Q2_skid <- SSN_Q1Q2$skid
SSN_Q3Q4_skid <- SSN_Q3Q4$skid
SSN_Q1Q2Q3Q4_skid <- SSN_Q1Q2Q3Q4$skid


# "presynapse (output)"

stats_synapse <- read.csv("analysis/data/stats_synapse.csv")

syn_cent_out_connectors <- stats_synapse %>%
  filter(skid==SSN_Q1Q2Q3Q4_skid) %>%
  filter(prepost==0) %>%
  select(connector_id) %>%
  pull()

syn_lateral_out_connectors  <- stats_synapse %>%
  filter(skid==SSN_Q1Q2_skid | skid==SSN_Q3Q4_skid) %>%
  filter(prepost==0) %>%
  select(connector_id) %>%
  pull()


# "postsynapse (input)"

syn_cent_to_lateral <- stats_synapse %>%
  filter(skid==SSN_Q1Q2_skid | skid==SSN_Q3Q4_skid) %>%
  filter(connector_id %in% syn_cent_out_connectors) %>%
  filter(prepost==1)

syn_cent_to_cent <- stats_synapse %>%
  filter(skid==SSN_Q1Q2Q3Q4_skid) %>%
  filter(connector_id %in% syn_cent_out_connectors) %>%
  filter(prepost==1)

syn_lateral_to_cent <- stats_synapse %>%
  filter(skid==SSN_Q1Q2Q3Q4_skid) %>%
  filter(connector_id %in% syn_lateral_out_connectors) %>%
  filter(prepost==1)

syn_lateral_to_lateral <- stats_synapse %>%
  filter(skid==SSN_Q1Q2_skid | skid==SSN_Q3Q4_skid) %>%
  filter(connector_id %in% syn_lateral_out_connectors) %>%
  filter(prepost==1)


stats_master <- read.csv("analysis/data/stats_master.csv")

# get outgoing syn from SSNs
conn_from_SSNs <- stats_synapse %>%
  filter(prepost == 0) %>%
  filter(skid == SSN_Q1Q2_skid | skid == SSN_Q3Q4_skid | skid == SSN_Q1Q2Q3Q4_skid) %>%
  select(connector_id) %>% 
  pull()

# skeletons that receive synapses from SSNs
SSN_downstream <- stats_synapse %>%
  filter(prepost==1) %>%
  filter(connector_id %in% conn_from_SSNs)

# in SSN downstream skid is skid of neuron postsynaptic to SSNs
# use get_celltype_annot_for_skid function from cell_statistics_master.R to get celltypes
#SSN_downstream <- SSN_downstream %>%
#  mutate(celltype = map_chr(skid, get_celltype_annot_for_skid))

# or get celltype from stats_master.csv
stats_master <- read.csv("analysis/data/stats_master.csv")
SSN_downstream <- SSN_downstream %>%
  left_join(stats_master %>% select(skid, celltype), by = "skid")



# bar plot
all_celltypes <- c("balancer", "bridge", "large_glanular_cell", "Cgroove", "dense_vesicle", "dome", 
                   "intra-multi-ciliated", "lamellate", "lithocyte", "plumose", "SSN", 
                   "epithelial_floor")

label_mapping <- c(
  "balancer" = "bal", "bridge" = "brg", "large_glanular_cell" = "lgc", 
  "Cgroove" = "cg", 
  "dense_vesicle" = "dv", "dome" = "do", 
  "intra-multi-ciliated" = "imc", "lamellate" = "la", "lithocyte" = "li", 
  "plumose" = "pl", "SSN" = "ANN", "epithelial_floor" = "ef"
)


# pre-synapse (take output-side connector)
syn_out_SSN_Q1Q2 <- stats_synapse %>%
  filter(prepost == 0, skid == SSN_Q1Q2_skid) %>%
  select(connector_id) %>%
  mutate(SSN_source = "SSN_Q1Q2")

syn_out_SSN_Q3Q4 <- stats_synapse %>%
  filter(prepost == 0, skid == SSN_Q3Q4_skid) %>%
  select(connector_id) %>%
  mutate(SSN_source = "SSN_Q3Q4")

syn_out_SSN_Q1Q2Q3Q4 <- stats_synapse %>%
  filter(prepost == 0, skid == SSN_Q1Q2Q3Q4_skid) %>%
  select(connector_id) %>%
  mutate(SSN_source = "SSN_Q1Q2Q3Q4")

# Putting together all the connector_id and SSN_source
connector_source_table <- bind_rows(
  syn_out_SSN_Q1Q2,
  syn_out_SSN_Q3Q4,
  syn_out_SSN_Q1Q2Q3Q4
)

# join to connector_id when taking postsynapse
SSN_downstream <- stats_synapse %>%
  filter(prepost == 1) %>%
  inner_join(connector_source_table, by = "connector_id") %>%
  left_join(stats_master %>% select(skid, celltype), by = "skid")



# plot
plot_output_number <- 
  SSN_downstream %>%
  filter(!celltype %in% c("monociliated", "biciliated", "multiciliated", "nonciliated", NA)) %>%
  ggplot(aes(x = factor(celltype, levels = all_celltypes), 
             fill = factor(SSN_source, levels = c("SSN_Q1Q2Q3Q4", "SSN_Q1Q2", "SSN_Q3Q4")))) +
  geom_bar(position = "stack", alpha = 0.75) +
  theme_minimal() +
  labs(
    x = "Postsynaptic cell types",
    y = "# of synapses from ANNs",
    fill = NULL  
  ) +
  theme(
    axis.text.x = element_text(size = 17, angle = -70, vjust = 0.5, hjust = 0, margin = margin(t = -7)),
    axis.text.y = element_text(size = 17),
    axis.title = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 13)
  ) +
  scale_x_discrete(limits = all_celltypes, labels = label_mapping) +
  scale_fill_manual(
    values = c(
      "SSN_Q1Q2Q3Q4" = Okabe_Ito[5],   
      "SSN_Q1Q2" = Okabe_Ito[6],
      "SSN_Q3Q4" = Okabe_Ito[7]
    ),
    labels = c(
      "SSN_Q1Q2Q3Q4" = "ANN Q1-4",
      "SSN_Q1Q2" = "ANN Q1Q2",
      "SSN_Q3Q4" = "ANN Q3Q4"
    )
  )

plot_output_number

ggsave(
  filename = "manuscript/pictures/output_from_SSNs.png",
  plot = plot_output_number,
  width = 2585,
  height = 1100,
  units = "px",
  dpi = 300,
  bg='white'
)

# save to source data (bar graph of outputs from SNN) -------------------------------------------------------

SSN_downstream %>%
  filter(!celltype %in% c("monociliated", "biciliated", "multiciliated", "nonciliated", NA)) %>%
  write_csv("manuscript/source_data/Figure3_source_data1.csv")

# plot balancer -------------------------------------------------------------


plot_balancer <- function() {
  plot3d(
    balancer_Q1,
    soma = T, lwd = 1.5, add = T,
    alpha = 0.5, col = Okabe_Ito[1]
  )
  plot3d(
    balancer_Q2,
    soma = T, lwd = 1.5, add = T,
    alpha = 0.5, col = Okabe_Ito[2]
  )
  plot3d(
    balancer_Q3,
    soma = T, lwd = 1.5, add = T,
    alpha = 0.5, col = Okabe_Ito[6]
  )
  plot3d(
    balancer_Q4,
    soma = T, lwd = 1.5, add = T,
    alpha = 0.5, col = Okabe_Ito[7]
  )
  
  plot3d(
    with_soma,
    soma = T, lwd = 1, add = T,
    alpha = 0.01, col = Okabe_Ito[8]
  )
  
  plot3d(
    outline,
    add = T, alpha = 0.05, col = "grey50"
  )
  par3d(zoom = 0.61)
}


close3d()

# 3d plotting of cells
nopen3d()
mfrow3d(1, 3) 
# define the size of the rgl window, the view and zoom
par3d(windowRect = c(0, 0, 1200, 350))
#par3d(windowRect = c(0, 0, 2400, 700))

# plot aboral view
plot_balancer()
aboral()


# move to next panel in rgl window
next3d(clear = F)

# plot lateral view of Sagittal plane
plot_balancer()
sagittal()


# move to next panel in rgl window
next3d(clear = F)

# plot lateral view of Tentacular plane
plot_balancer()
tentacular()

# make a snapshot to the working directory
rgl.snapshot("manuscript/pictures/balancer.png")
close3d()





# plot bridge -----------------------------------------------

plot_bridge <- function() {
  plot3d(
    bridge_Q1Q2,
    soma = T, lwd = 1.5, add = T,
    alpha = 0.75, col = Okabe_Ito[2]
  )
  plot3d(
    bridge_Q3Q4,
    soma = T, lwd = 1.5, add = T,
    alpha = 0.75, col = Okabe_Ito[7]
  )
  
  plot3d(
    with_soma,
    soma = T, lwd = 1, add = T,
    alpha = 0.01, col = Okabe_Ito[8]
  )
  
  plot3d(
    outline,
    add = T, alpha = 0.05, col = "grey50"
  )
  par3d(zoom = 0.61)
}

close3d()

# 3d plotting of cells
nopen3d()
mfrow3d(1, 3) 
# define the size of the rgl window, the view and zoom
par3d(windowRect = c(0, 0, 1200, 350))
#par3d(windowRect = c(0, 0, 2400, 700))

# plot aboral view
plot_bridge()
aboral()

# move to next panel in rgl window
next3d(clear = F)

# plot lateral view of Sagittal plane
plot_bridge()
sagittal()

# move to next panel in rgl window
next3d(clear = F)

# plot lateral view of Tentacular plane
plot_bridge()
tentacular()

#make a snapshot to the working directory
rgl.snapshot("manuscript/pictures/bridge.png")


close3d()


# plot bridge and balancer -----------------------------------------------------

plot_bridge_and_balancer <- function() {
  plot3d(balancer,
         soma = T, lwd = 1, add = T, 
         alpha = 0.1, col = Okabe_Ito[8],
         WithConnectors = F, WithNodes = F)
  
  plot3d(
    bridge_Q1Q2,
    soma = T, lwd = 1.5, add = T,
    alpha = 0.75, col = Okabe_Ito[2]
  )
  plot3d(
    bridge_Q3Q4,
    soma = T, lwd = 1.5, add = T,
    alpha = 0.75, col = Okabe_Ito[7]
  )
  
  plot3d(
    with_soma,
    soma = T, lwd = 1, add = T,
    alpha = 0.01, col = Okabe_Ito[8]
  )
  
  plot3d(
    outline,
    add = T, alpha = 0.05, col = "grey50"
  )
  par3d(zoom = 0.61)
}


close3d()

# 3d plotting of cells
nopen3d()
mfrow3d(1, 3) 
# define the size of the rgl window, the view and zoom
par3d(windowRect = c(0, 0, 1200, 350))

# plot aboral view
plot_bridge_and_balancer()
aboral()

# move to next panel in rgl window
next3d(clear = F)

# plot lateral view of Sagittal plane
plot_bridge_and_balancer()
sagittal()

# move to next panel in rgl window
next3d(clear = F)

# plot lateral view of Tentacular plane
plot_bridge_and_balancer()
tentacular()

#make a snapshot to the working directory
rgl.snapshot("manuscript/pictures/bridge_and_balancer.png")

close3d()



# load of mitochondrial location information-------------------------------

mito_vesicle_info <- read.csv("analysis/data/mito_vesicle_info.csv")

# 3D plot mitochondria positions in bridge ----------------------------------------

pos_ves_syn <- mito_vesicle_info |>
  filter(celltype == "bridge") |>
  filter(mito_type == "vesicles_syn") |>
  select(x, y, z)

pos_ves_no_syn <- mito_vesicle_info |>
  filter(celltype == "bridge") |>
  filter(mito_type != "vesicles_syn") |>
  select(x, y, z)


plot_mito_bridge <- function() {
  plot3d(bridge_Q1Q2,
         soma = FALSE, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[2],
         WithConnectors = F, WithNodes = F)
  
  plot3d(bridge_Q3Q4,
         soma = FALSE, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[7],
         WithConnectors = F, WithNodes = F)
  
  plot3d(outline,
         add = T, alpha = 0.05, col = "grey50"
  )
  
  # plot mitochondria
  plot3d(pos_ves_syn, 
         add = TRUE,
         col = "green3",
         size = 1.1, 
         alpha = 1,
         point_antialias = TRUE,
         type = "s"
  )
  
  plot3d(pos_ves_no_syn, 
         add = TRUE,
         col = "black",
         size = 0.9, 
         alpha = 0.37,
         point_antialias = TRUE,
         type = "s"
  )
  par3d(zoom=0.61)
}

close3d()
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  #defines the two scenes
#define the size of the rgl window, the view and zoom
par3d(windowRect = c(0, 0, 1200, 350))
#par3d(windowRect = c(0, 0, 2400, 700))

#plot aboral view
plot_mito_bridge()
aboral()

#move to next panel in rgl window
next3d(clear=F)

#plot sagittal view
plot_mito_bridge()
sagittal()

#move to next panel in rgl window
next3d(clear=F)

#plot tentacular view
plot_mito_bridge()
tentacular()

#make a snapshot to the working directory
rgl.snapshot("manuscript/pictures/mito_pos_bridge.png")

close3d()



# retrieve connectors ----------------------------------------------------------

# starting from balancers doesn't give me accurate synapse count?

outgoing_SSN_Q1Q2 <- stats_synapse %>%
  filter(skid == SSN_Q1Q2_skid) %>%
  filter(prepost == 0) %>%
  select(connector_id) %>%
  pull()

SSN_Q1Q2_to_balancer <- stats_synapse %>%
  filter(connector_id %in% outgoing_SSN_Q1Q2) %>%
  filter(prepost == 1) %>%
  filter(skid %in% balancer_skids)

outgoing_SSN_Q3Q4 <- stats_synapse %>%
  filter(skid == SSN_Q3Q4_skid) %>%
  filter(prepost == 0) %>%
  select(connector_id) %>%
  pull()

SSN_Q3Q4_to_balancer <- stats_synapse %>%
  filter(connector_id %in% outgoing_SSN_Q3Q4) %>%
  filter(prepost == 1) %>%
  filter(skid %in% balancer_skids)

outgoing_SSN_Q1Q2Q3Q4 <- stats_synapse %>%
  filter(skid == SSN_Q1Q2Q3Q4_skid) %>%
  filter(prepost == 0) %>%
  select(connector_id) %>%
  pull()

SSN_Q1Q2Q3Q4_to_balancer <- stats_synapse %>%
  filter(connector_id %in% outgoing_SSN_Q1Q2Q3Q4) %>%
  filter(prepost == 1) %>%
  filter(skid %in% balancer_skids)


bridge_conn <- connectors(bridge)
presyn_bridge <- subset(bridge_conn, prepost == 0)
postsyn_bridge <- subset(bridge_conn, prepost == 1)


# 3D plot of synaptic input from SSN to balancer -------------------------------

plot_SNN_to_balancer <- function() {
  plot3d(balancer_Q1,
         soma = FALSE, lwd = 1, add = T, alpha = 0.1, col = Okabe_Ito[8],
         WithConnectors = F, WithNodes = F)
  
  plot3d(balancer_Q2,
         soma = FALSE, lwd = 1, add = T, alpha = 0.1, col = Okabe_Ito[8],
         WithConnectors = F, WithNodes = F)
  
  plot3d(balancer_Q3,
         soma = FALSE, lwd = 1, add = T, alpha = 0.1, col = Okabe_Ito[8],
         WithConnectors = F, WithNodes = F)
  
  plot3d(balancer_Q4,
         soma = FALSE, lwd = 1, add = T, alpha = 0.1, col = Okabe_Ito[8],
         WithConnectors = F, WithNodes = F)
  
  plot3d(SSN_Q1Q2Q3Q4,
         soma = FALSE, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
         WithConnectors = F, WithNodes = F)
  
  plot3d(SSN_Q1Q2,
         soma = FALSE, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[6],
         WithConnectors = F, WithNodes = F)
  
  plot3d(SSN_Q3Q4,
         soma = FALSE, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[7],
         WithConnectors = F, WithNodes = F)
  
  plot3d(outline,
         add = T, alpha = 0.05, col = "grey50"
  )
  
  # input from SSN to balancer
  plot3d(
    SSN_Q1Q2_to_balancer$x, 
    SSN_Q1Q2_to_balancer$y,
    SSN_Q1Q2_to_balancer$z,
    size = 0.8, alpha = 0.7, col = "orange", 
    add = TRUE,
    point_antialias = TRUE,
    type = "s"
  )
  
  plot3d(
    SSN_Q3Q4_to_balancer$x, 
    SSN_Q3Q4_to_balancer$y,
    SSN_Q3Q4_to_balancer$z,
    size = 0.8, alpha = 0.7, col = "magenta", 
    add = TRUE,
    point_antialias = TRUE,
    type = "s"
  )
  
  plot3d(
    SSN_Q1Q2Q3Q4_to_balancer$x, 
    SSN_Q1Q2Q3Q4_to_balancer$y,
    SSN_Q1Q2Q3Q4_to_balancer$z,
    size = 0.8, alpha = 0.7, col = "dodgerblue2", 
    add = TRUE,
    point_antialias = TRUE,
    type = "s"
  )
  par3d(zoom=0.61)
}


close3d()
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  
#define the size of the rgl window, the view and zoom
par3d(windowRect = c(0, 0, 1200, 350))
#par3d(windowRect = c(0, 0, 2400, 700))

#plot aboral view
plot_SNN_to_balancer()
aboral()

#move to next panel in rgl window
next3d(clear=F)

#plot lateral view of Sagittal plane
plot_SNN_to_balancer()
sagittal()

next3d(clear=F)

#plot lateral view of Tentacular plane
plot_SNN_to_balancer()
tentacular()

#make a snapshot to the working directory
rgl.snapshot("manuscript/pictures/SSN_prepost_synapse_balancer.png")

close3d()



# 3D plot of reciprocal synaptic inputs between SSN and bridge-----------------------------

plot_reciprocal_SSN_syn <- function() {
  plot3d(bridge_Q1Q2,
         soma = FALSE, lwd = 1, add = T, alpha = 0.2, col = Okabe_Ito[8],
         WithConnectors = F, WithNodes = F)
  
  plot3d(bridge_Q3Q4,
         soma = FALSE, lwd = 1, add = T, alpha = 0.2, col = Okabe_Ito[8],
         WithConnectors = F, WithNodes = F)
  
  plot3d(SSN_Q1Q2Q3Q4,
         soma = FALSE, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
         WithConnectors = F, WithNodes = F)
  
  plot3d(SSN_Q1Q2,
         soma = FALSE, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[6],
         WithConnectors = F, WithNodes = F)
  
  plot3d(SSN_Q3Q4,
         soma = FALSE, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[7],
         WithConnectors = F, WithNodes = F)
  
  plot3d(outline,
         add = T, alpha = 0.05, col = "grey50"
  )
  
  # inputs from SSN to bridge
  plot3d(
    postsyn_bridge$x, 
    postsyn_bridge$y, 
    postsyn_bridge$z, 
    size = 0.8, alpha = 0.7, col = "magenta2", 
    add = TRUE,
    point_antialias = TRUE,
    type = "s"
  )
  
  # inputs from bridge to SNN
  plot3d(
    presyn_bridge$x, 
    presyn_bridge$y, 
    presyn_bridge$z, 
    size = 0.8, alpha = 1, col = "cyan2", 
    add = TRUE,
    point_antialias = TRUE,
    type = "s"
  )
  par3d(zoom=0.61)
}

close3d()
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  #defines the two scenes
#define the size of the rgl window, the view and zoom
par3d(windowRect = c(0, 0, 1200, 350))
#par3d(windowRect = c(0, 0, 2400, 700))

#plot aboral view
plot_reciprocal_SSN_syn()
aboral()

next3d(clear=F)

#plot lateral view of Sagittal plane
plot_reciprocal_SSN_syn()
sagittal()
next3d(clear=F)


#plot lateral view of Tentacular plane
plot_reciprocal_SSN_syn()
tentacular()
next3d(clear=F)


#make a snapshot to the working directory
rgl.snapshot("manuscript/pictures/SSN_prepost_synapse_bridge.png")

close3d()





# get connectivity for matrix plot ---------------------------------------------

SSN_Q1Q2 <- read_smooth_neuron("SSN_Q1Q2")
SSN_Q3Q4 <- read_smooth_neuron("SSN_Q3Q4")
SSN_Q1Q2Q3Q4 <- read_smooth_neuron("SSN_Q1Q2Q3Q4")


matrix_celltypes <- list(balancer_Q1,
                         balancer_Q2,
                         balancer_Q3,
                         balancer_Q4,
                         bridge_Q1Q2,
                         bridge_Q3Q4,
                         SSN_Q1Q2,
                         SSN_Q3Q4,
                         SSN_Q1Q2Q3Q4)

matrix_celltypes_names <- list("balancer_Q1",
                               "balancer_Q2",
                               "balancer_Q3",
                               "balancer_Q4",
                               "bridge_Q1Q2",
                               "bridge_Q3Q4",
                               "SSN_Q1Q2",
                               "SSN_Q3Q4",
                               "SSN_Q1Q2Q3Q4")

# define empty synapse list with the right dimensions
synapse_list <- c()

for (i in 1:length(matrix_celltypes)) {
  for (j in 1:length(matrix_celltypes)) {
    # get connectors between two cell groups
    presyn_skids <- attr(matrix_celltypes[i][[1]], "df")$skid
    postsyn_skids <- attr(matrix_celltypes[j][[1]], "df")$skid
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
  nrow = length(matrix_celltypes)
)

rownames(synapse_matrix) <- as.character(matrix_celltypes_names)
colnames(synapse_matrix) <- as.character(matrix_celltypes_names)

syn_df <- as.data.frame(as.table(synapse_matrix))

# matrix plot ------------------------------------------------------------------

rename_map <- c("balancer_Q1" = "bal Q1", 
                "balancer_Q2" = "bal Q2", 
                "balancer_Q3" = "bal Q3", 
                "balancer_Q4" = "bal Q4",
                "bridge_Q1Q2" = "brg Q1Q2",
                "bridge_Q3Q4" = "brg Q3Q4",
                "SSN_Q1Q2" = "ANN Q1Q2", 
                "SSN_Q3Q4" = "ANN Q3Q4", 
                "SSN_Q1Q2Q3Q4" = "ANN Q1-4")

custom_order <- c("bal Q1", "bal Q2", "brg Q1Q2", "ANN Q1Q2", 
                  "bal Q3", "bal Q4", "brg Q3Q4", "ANN Q3Q4", 
                  "ANN Q1-4")

syn_df <- syn_df %>%
  mutate(
    Var1 = factor(recode(Var1, !!!rename_map), levels = rev(custom_order)),
    Var2 = factor(recode(Var2, !!!rename_map), levels = custom_order)
  )


plot_syn_matrix <- ggplot(syn_df, aes(x = Var2, y = Var1, fill = Freq)) +
  geom_tile(color = "lightgrey") + 
  geom_text(aes(label = ifelse(Freq == 0, "", Freq)), 
            color = "black",
            size = 5) +
  scale_fill_gradient(name = "Number",
                      low = "white", high = "#0072b2", 
                      limits = c(0, 30), 
                      oob = scales::squish,
                      labels = function(x) ifelse(x == 30, "30>", as.character(x))) + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 16, angle = -45, hjust = 0),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18, margin = margin(t = 10)),
    axis.title.y = element_text(size = 18, margin = margin(r = 10)),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15, margin = margin(b = 12))
  ) +
  labs(
    x = "Postsynaptic",
    y = "Presynaptic"
    )


plot_syn_matrix


ggsave(
  filename = "manuscript/pictures/balancers_bridges_matrix.png",
  plot = plot_syn_matrix,
  width = 2000,
  height = 1700,
  units = "px",
  dpi = 300,
  bg='white'
)

# save to source data (connectivity for matrix plot) -------------------------------------------------------

syn_df %>%
  write_csv("manuscript/source_data/Figure3_source_data2.csv")


# assemble figure -------------------------------------------------------------

panel_output <- ggdraw() + draw_image(readPNG("manuscript/pictures/output_from_SSNs.png"))

panel_bal <- ggdraw() + draw_image(readPNG("manuscript/pictures/balancer.png")) +
  draw_label("Balancers", x = 0.5, y = 1.05, color = "black", size = 10, hjust = 0.5) +
  draw_label("Q1", x = 0.25, y = 0.84, color = Okabe_Ito[1], size = 9, hjust = 0, alpha = 1) +
  draw_label("Q2", x = 0.1, y = 0.77, color = Okabe_Ito[2], size = 9, hjust = 1, alpha = 1) +
  draw_label("Q3", x = 0.075, y = 0.13, color = Okabe_Ito[6], size = 9, hjust = 1, alpha = 1) +
  draw_label("Q4", x = 0.22, y = 0.2, color = Okabe_Ito[7], size = 9, hjust = 0, alpha = 1) +
  draw_label("aboral view", x = 0.18, y = 0.95, color = "black", size = 8, hjust = 0.5) +
  draw_label("lateral view", x = 0.67, y = 0.95, color = "black", size = 8, hjust = 0.5) +
  draw_label("sagittal plane", x = 0.5, y = 0.86, color = "black", size = 7.5, hjust = 0.5) +
  draw_label("tentacular plane", x = 0.83, y = 0.86, color = "black", size = 7.5, hjust = 0.5) +
  draw_label("A", x = 0.67, y = 0.85, size = 7, color = "black", hjust = 0.5) +
  draw_label("O", x = 0.67, y = 0.65, size = 7, color = "black", hjust = 0.5) +
  draw_line(x = c(0.67, 0.67), y = c(0.69, 0.81), color = "black", linewidth = 0.65) +
  draw_label(expression(paste("25 ", mu, "m")), x = 0.95, y = 0.1, color = "black", size = 7, hjust = 0.5) +
  draw_line(x = c(0.91, 0.99), y = c(0.05, 0.05), color = "black", linewidth = 0.7)

panel_bri <- ggdraw() + draw_image(readPNG("manuscript/pictures/bridge.png")) +
  draw_label("Bridges", x = 0.5, y = 0.95, color = "black", size = 10, hjust = 0.5) +
  draw_label("Q1Q2", x = 0.2, y = 0.86, color = Okabe_Ito[2], size = 9, hjust = 0, alpha = 1) +
  draw_label("Q3Q4", x = 0.19, y = 0.18, color = Okabe_Ito[7], size = 9, hjust = 0, alpha = 1)

circle_g <- ggdraw() + draw_image(readPNG("manuscript/pictures/synapse_green.png"))
circle_gy <- ggdraw() + draw_image(readPNG("manuscript/pictures/synapse_gray.png"))

panel_bri_mito <- ggdraw() + draw_image(readPNG("manuscript/pictures/mito_pos_bridge.png")) +
  draw_plot(circle_g, x = 0.06, y = 0.995, width = 0.05, height = 0.05) +
  draw_label("With synapses", x = 0.1, y = 1.02, color = "green4", size = 8.5, hjust = 0) +
  draw_plot(circle_gy, x = 0.06, y = 0.915, width = 0.05, height = 0.05) +
  draw_label("Not forming synapses", x = 0.1, y = 0.94, color = "black", size = 8.5, hjust = 0, alpha = 0.7) +
  draw_label("bridge_Q1Q2", x = 0.95, y = 1.015, color = Okabe_Ito[2], size = 7.5, hjust = 1) +
  draw_label("bridge_Q3Q4", x = 0.95, y = 0.92, color = Okabe_Ito[7], size = 7.5, hjust = 1)

circle_dob <- ggdraw() + draw_image(readPNG("manuscript/pictures/synapse_dodgerblue.png"))
circle_or <- ggdraw() + draw_image(readPNG("manuscript/pictures/synapse_orange.png"))
circle_m <- ggdraw() + draw_image(readPNG("manuscript/pictures/synapse_magenta.png"))

panel_SSN_bal <- ggdraw() + draw_image(readPNG("manuscript/pictures/SSN_prepost_synapse_balancer.png")) +
  draw_label("→ Balancers", x = 0.505, y = 1.025, color = "gray30", size = 8, hjust = 0) +
  draw_plot(circle_dob, x = 0.065, y = 1, width = 0.05, height = 0.05) +
  draw_plot(circle_or, x = 0.2, y = 1, width = 0.05, height = 0.05) +
  draw_plot(circle_m, x = 0.345, y = 1, width = 0.05, height = 0.05) +
  draw_label("ANN Q1-4", x = 0.1, y = 1.025, color = Okabe_Ito[5], size = 8, hjust = 0) +
  draw_label("ANN Q1Q2", x = 0.235, y = 1.025, color = Okabe_Ito[6], size = 8, hjust = 0) +
  draw_label("ANN Q3Q4", x = 0.38, y = 1.025, color = Okabe_Ito[7], size = 8, hjust = 0)

circle_m <- ggdraw() + draw_image(readPNG("manuscript/pictures/synapse_magenta.png"))
circle_lb <- ggdraw() + draw_image(readPNG("manuscript/pictures/synapse_lightblue.png"))

panel_SSN_bri <- ggdraw() + draw_image(readPNG("manuscript/pictures/SSN_prepost_synapse_bridge.png")) +
  draw_plot(circle_m, x = 0.06, y = 0.975, width = 0.05, height = 0.05) +
  draw_label("ANN → Bridges", x = 0.1, y = 1, color = "gray30", size = 8, hjust = 0) +
  draw_plot(circle_lb, x = 0.26, y = 0.975, width = 0.05, height = 0.05) +
  draw_label("Bridges → ANN", x = 0.3, y = 1, color = "gray30", size = 8, hjust = 0) +
  draw_label("ANN Q1-4", x = 0.85, y = 1.13, color = Okabe_Ito[5], size = 7, hjust = 0) +
  draw_label("ANN Q1Q2", x = 0.85, y = 1.06, color = Okabe_Ito[6], size = 7, hjust = 0) +
  draw_label("ANN Q3Q4", x = 0.85, y = 0.99, color = Okabe_Ito[7], size = 7, hjust = 0)

panel_matrix <- ggdraw() + draw_image(readPNG("manuscript/pictures/balancers_bridges_matrix.png"))
panel_map_graph <- ggdraw() + draw_image(readPNG("manuscript/pictures/map_balancer_bridge_SSN_1.png"))



layout <- "
####
A#EE
####
B#FF
C###
C#GH
##GH
D#GH
"

Figure3 <- panel_output + panel_bal + panel_bri + panel_bri_mito +
  panel_SSN_bal + panel_SSN_bri + 
  panel_matrix + panel_map_graph +
  plot_layout(design = layout,
              heights = c(0.05, 1.1, 0.2, 1.1, 0.1, 0.9, 0.1, 0.9),
              widths = c(1.6, 0.1, 1.1, 0.9)) + 
  patchwork::plot_annotation(tag_levels = "A") &  
  ggplot2::theme(plot.tag = element_text(size = 12, 
                                         face='plain', color='black'))


ggsave("manuscript/figures/Figure3.png", limitsize = FALSE, 
       units = c("px"), Figure3, width = 3000, height = 2100, bg='white')  


ggsave("manuscript/figures/Figure3.pdf", limitsize = FALSE, 
       units = c("px"), Figure3, width = 3000, height = 2100) 


