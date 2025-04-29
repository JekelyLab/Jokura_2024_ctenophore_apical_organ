# Code to generate Figure 2 of the Jokura et al 2024 Ctenophore apical organ connectome paper

# source packages and functions ------------------------------------------------
source("analysis/scripts/packages_and_functions.R")


# load cell type ---------------------------------------------------------------

#balancer <- read_smooth_neuron("celltype:balancer")
#bridge <- read_smooth_neuron("celltype:bridge")
#bristle <- read_smooth_neuron("celltype:bristle")
#Cgroove <- read_smooth_neuron("celltype:Cgroove")
#dense_vesicle <- read_smooth_neuron("celltype:dense_vesicle")
#dome <- read_smooth_neuron("celltype:dome")
#epithelial_floor <- read_smooth_neuron("celltype:epithelial_floor")
#intra_multi_ciliated <- read_smooth_neuron("celltype:intra-multi-ciliated")
#lamellate <- read_smooth_neuron("celltype:lamellate")
#lithocyte <- read_smooth_neuron("celltype:lithocyte")
#plumose <- read_smooth_neuron("celltype:plumose")
#SSN <- read_smooth_neuron("celltype:SSN")

#monociliated <- read_smooth_neuron("celltype:monociliated")
#biciliated <- read_smooth_neuron("celltype:biciliated")
#multiciliated <- read_smooth_neuron("celltype:multiciliated")
#nonciliated <- read_smooth_neuron("celltype:nonciliated")

#Q1 <- read_smooth_neuron("Q1")
#Q2 <- read_smooth_neuron("Q2")
#Q3 <- read_smooth_neuron("Q3")
#Q4 <- read_smooth_neuron("Q4")

SSN_Q1Q2 <- read_smooth_neuron("SSN_Q1Q2")[[1]]
SSN_Q3Q4 <- read_smooth_neuron("SSN_Q3Q4")[[1]]
SSN_Q1Q2Q3Q4 <- read_smooth_neuron("SSN_Q1Q2Q3Q4")[[1]]

with_soma <- read_smooth_neuron("with_soma")

#all_celltypes <- list(balancer,bridge,bristle,Cgroove,
#                      dense_vesicle,dome,epithelial_floor,
#                      intra_multi_ciliated,
#                      lamellate,lithocyte,plumose,SSN,
#                      monociliated,biciliated,multiciliated,nonciliated)



# 3d plot SSN Q1Q2 & Q3Q4 neuron -----------------------------------------------

close3d()
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  
#define the size of the rgl window, the view and zoom
par3d(windowRect = c(0, 0, 1200, 350))
#par3d(windowRect = c(0, 0, 2400, 700))

#plot aboral view
plot_multinucleated_cell(SSN_Q1Q2,
                         lwd = 1, alpha = 1, col = Okabe_Ito[6])

plot_multinucleated_cell(SSN_Q3Q4,
                         lwd = 1, alpha = 1, col = Okabe_Ito[7])

plot3d(outline,
       add = T, alpha = 0.05, col = "grey50"
)

#texts3d(58000,51000,5000, "INRGW", cex = 3, col = "#56B4E9")

#aboral view
nview3d("anterior", extramat = rotationMatrix(1.05, 250, -200, 1000))
#rgl.snapshot("manuscript/pictures/SSN_aboral_view.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Sagittal plane
plot_multinucleated_cell(SSN_Q1Q2,
                         lwd = 1, alpha = 1, col = Okabe_Ito[6])

plot_multinucleated_cell(SSN_Q3Q4,
                         lwd = 1, alpha = 1, col = Okabe_Ito[7])

plot3d(with_soma,
       soma = T, lwd = 0.5, add = T, alpha = 0.025, col = Okabe_Ito[8],
       WithConnectors = F, WithNodes = F)

plot3d(outline,
       add = T, alpha = 0.05, col = "grey50"
)

#texts3d(58000,51000,5000, "INRGW", cex = 3, col = "#56B4E9")

#lateral view of Sagittal plane
nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
#rgl.snapshot("manuscript/pictures/SSN_sagittal_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Tentacular plane
plot_multinucleated_cell(SSN_Q1Q2,
                         lwd = 1, alpha = 1, col = Okabe_Ito[6])

plot_multinucleated_cell(SSN_Q3Q4,
                         lwd = 1, alpha = 1, col = Okabe_Ito[7])

plot3d(with_soma,
       soma = T, lwd = 0.5, add = T, alpha = 0.025, col = Okabe_Ito[8],
       WithConnectors = F, WithNodes = F)

plot3d(outline,
       add = T, alpha = 0.05, col = "grey50"
)


#lateral view of Tentacular plane
nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))
#rgl.snapshot("manuscript/pictures/SSN_tentacular_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)


#make a snapshot to the working directory
rgl.snapshot("manuscript/pictures/SSN_Q12_Q34.png")


close3d()

# 3d plot SSN Q1Q2Q3Q4 neuron -----------------------------------------------


close3d()
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  
#define the size of the rgl window, the view and zoom
par3d(windowRect = c(0, 0, 1200, 350))
#par3d(windowRect = c(0, 0, 2400, 700))

#plot aboral view
plot_multinucleated_cell(SSN_Q1Q2Q3Q4,
       lwd = 1, alpha = 1, col = Okabe_Ito[5])

plot3d(with_soma,
       soma = T, lwd = 0.5, add = T, alpha = 0.025, col = Okabe_Ito[8],
       WithConnectors = F, WithNodes = F)

plot3d(outline,
       add = T, alpha = 0.05, col = "grey50"
)

#texts3d(58000,51000,5000, "INRGW", cex = 3, col = "#56B4E9")

#aboral view
nview3d("anterior", extramat = rotationMatrix(1.05, 250, -200, 1000))
#rgl.snapshot("manuscript/pictures/SSN_aboral_view.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Sagittal plane
plot_multinucleated_cell(SSN_Q1Q2Q3Q4,
                         lwd = 1, alpha = 1, col = Okabe_Ito[5])

plot3d(with_soma,
       soma = T, lwd = 0.5, add = T, alpha = 0.025, col = Okabe_Ito[8],
       WithConnectors = F, WithNodes = F)

plot3d(outline,
       add = T, alpha = 0.05, col = "grey50"
)

#texts3d(58000,51000,5000, "INRGW", cex = 3, col = "#56B4E9")

#lateral view of Sagittal plane
nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
#rgl.snapshot("manuscript/pictures/SSN_sagittal_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Tentacular plane
plot_multinucleated_cell(SSN_Q1Q2Q3Q4,
                         lwd = 1, alpha = 1, col = Okabe_Ito[5])

plot3d(with_soma,
       soma = T, lwd = 0.5, add = T, alpha = 0.025, col = Okabe_Ito[8],
       WithConnectors = F, WithNodes = F)

plot3d(outline,
       add = T, alpha = 0.05, col = "grey50"
)


#lateral view of Tentacular plane
nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))
#rgl.snapshot("manuscript/pictures/SSN_tentacular_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)


#make a snapshot to the working directory
rgl.snapshot("manuscript/pictures/SSN_Q1234.png")


close3d()




# calculate percentage of mitochondria which have vesicles or synapses associated with them -----

mito_done <- read.neurons.catmaid("mitochondria done", pid = 35)

get_pos_of_tags_in_neuron <- function(neur, tagname) {
  # argument is neuron (data type), not skid
  tag_positions <- tibble(
    treenodeid = character(),
    x = double(),
    y = double(),
    z = double()
  )
  tags <- neur$tags[[tagname]]
  for (tag in tags) {
    pos <- neur$d |>
      filter(PointNo == tag) |>
      select(X, Y, Z)
    tag <- as.character(tag)
    tag_positions <- tag_positions |> bind_rows(list(
      treenodeid = tag,
      x = pos$X, y = pos$Y, z = pos$Z
    ))
  }
  return(tag_positions)
}

get_mito_pos <- function(neur) {
  # accepts neuron as input
  sskid <- neur$NeuronName
  celltype <- catmaid_get_annotations_for_skeletons(sskid, pid = 35) |>
    select(annotation) |>
    filter(grepl("celltype:", annotation)) |>
    pull()
  celltype <- gsub(".*:", "", celltype)
  if (length(celltype) == 0) {
    celltype <- "NA"
  }
  
  ves_none <- get_pos_of_tags_in_neuron(neur, "mitochondrion no vesicles")
  ves_none <- bind_cols(mito_type = "vesicles_none", ves_none)
  
  if (nrow(ves_none) == 0) {
    ves_none <- get_pos_of_tags_in_neuron(neur, "mitochondrion")
    ves_none <- bind_cols(mito_type = "vesicles_none", ves_none)
  }
  
  ves_unc <- get_pos_of_tags_in_neuron(neur, "mitochondrion unclear vesicles")
  ves_unc <- bind_cols(mito_type = "vesicles_unc", ves_unc)
  
  ves_yes <- get_pos_of_tags_in_neuron(neur, "mitochondrion vesicles")
  ves_yes <- bind_cols(mito_type = "vesicles", ves_yes)
  # treenodes with outgoing synapses
  connectors <- neur$connectors
  if (!is.null(connectors)) {
    syn_out_treenodes <- connectors |>
      filter(prepost == 0) |>
      select(treenode_id) |>
      pull()
  } else {
    syn_out_treenodes <- integer()
  }
  # check if the mitochondrion is associated with a synapse
  for (treenode in ves_yes$treenodeid) {
    if (treenode %in% syn_out_treenodes) {
      mitotype <- "vesicles_syn"
    } else {
      mitotype <- "vesicles_no_syn"
    }
    ves_yes <- ves_yes |>
      mutate(mito_type = replace(mito_type, treenodeid == treenode, mitotype))
  }
  mito_pos <- bind_rows(ves_none, ves_unc, ves_yes)
  mito_pos <- bind_cols(
    celltype = celltype,
    skid = sskid,
    mito_pos
  )
  return(mito_pos)
}

mito_vesicle_info <- lapply(mito_done, get_mito_pos) |>
  bind_rows()

write.csv(mito_vesicle_info, "analysis/data/mito_vesicle_info.csv")

# load of mitochondrial data csv file ------------------------------------------

df <- read_csv("analysis/data/mito_vesicle_info.csv")

# mitochondrial data csv file ------------------------------------------

df <- df %>%
  mutate(synapse_related = ifelse(mito_type == "vesicles_syn", "mean_vesicles_syn", "mean_vesicles_no_syn"))

exclude_types <- c("biciliated", "monociliated", "multiciliated", "nonciliated", "SNN", NA)

df <- df %>%
  filter(!(celltype %in% exclude_types))

rename_map <- c(
  "balancer" = "bal", 
  "bridge" = "brg", 
  "bristle" = "bsl", 
  "Cgroove" = "cg",
  "dense_vesicle" = "dv", 
  "dome" = "do", 
  "intra-multi-ciliated" = "imc", 
  "lamellate" = "la", 
  "lithocyte" = "li", 
  "plumose" = "pl", 
  "SSN" = "ANN", 
  "epithelial_floor" = "ef"
)

df <- df %>%
  mutate(celltype = recode(celltype, !!!rename_map))

# calculate averages using mitochondrial data and format the data --------------

# Aggregation of individual cell data
summary_df <- df %>%
  group_by(celltype, skid, synapse_related) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = synapse_related, values_from = count, values_fill = 0)

# Creation of average data (for bar graphs)
mean_summary <- summary_df %>%
  group_by(celltype) %>%
  summarise(
    mean_vesicles_syn = mean(mean_vesicles_syn),
    mean_vesicles_no_syn = mean(mean_vesicles_no_syn),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(mean_vesicles_syn, mean_vesicles_no_syn), names_to = "characteristic", values_to = "value")


# Align the number of cells (n) of each celltype with celltype_order
celltype_order <- c("bal", "brg", "bsl", "cg", "dv", "do", "imc", "la", "li", "pl", "ANN", "ef")

n_counts <- summary_df %>%
  group_by(celltype) %>%
  summarise(n = n()) %>%
  arrange(factor(celltype, levels = celltype_order))

n_labels <- paste0(n_counts$celltype, " (", n_counts$n, ")")

mean_summary$celltype <- factor(mean_summary$celltype, levels = celltype_order)
summary_df$celltype <- factor(summary_df$celltype, levels = celltype_order)


# plot Average number of mitochondria per cell -----------------------------------------------------------

plot_mito_stats <- 
  ggplot() +
  geom_bar(data = mean_summary, aes(
    x = celltype, 
    y = value, 
    fill = factor(characteristic, 
                  levels = c("mean_vesicles_syn", "mean_vesicles_no_syn"))), 
    stat = "identity", position = "stack", width = 0.7) +
  geom_jitter(data = summary_df, 
              aes(x = celltype, y = mean_vesicles_no_syn), 
              width = 0.2, color = "black", size = 1, alpha = 0.5) +
  geom_jitter(data = summary_df %>% filter(mean_vesicles_syn > 0), 
              aes(x = celltype, y = mean_vesicles_syn + mean_vesicles_no_syn), 
              width = 0.2, color = "green3", size = 1.5, alpha = 0.75) +
  scale_fill_manual(
    values = c("mean_vesicles_syn" = "green4", "mean_vesicles_no_syn" = "darkgrey"),
    labels = c("mean_vesicles_syn" = "mitochondria with synapses", 
               "mean_vesicles_no_syn" = "mitochondria not forming synapses")) +
  theme_bw() +
#  scale_x_discrete(labels = n_labels) +
  theme(
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(size = 15, angle = -70, vjust = 0.5, hjust = 0, margin = margin(t = -7)),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 13, margin = margin(r = 15)),
    legend.title = element_blank(),
    legend.position = c(0.3, 0.75),
    text = element_text(size = 15)
  ) +
  ylab("Mean mitochondria count per cell") +
  xlab("Cell types")

ggsave("manuscript/pictures/plot_mito_stats.png", 
       plot = plot_mito_stats,
       width = 2400,
       height = 1100,
       units = "px",
       dpi = 300)


# load of mitochondrial location information-------------------------------

#mito_done <- read.neurons.catmaid("mitochondria done", pid = 35)
mito_vesicle_info <- read.csv("analysis/data/mito_vesicle_info.csv")

# 3d plot mitochondria positions in SSN ----------------------------------------

close3d()
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  
#define the size of the rgl window, the view and zoom
par3d(windowRect = c(0, 0, 1200, 350))
#par3d(windowRect = c(0, 0, 2400, 700))

#plot aboral view
plot3d(SSN_Q1Q2,
       soma = FALSE, lwd = 1, add = T, alpha = 0.25, col = Okabe_Ito[6],
       WithConnectors = F, WithNodes = F)

plot3d(SSN_Q3Q4,
       soma = FALSE, lwd = 1, add = T, alpha = 0.25, col = Okabe_Ito[7],
       WithConnectors = F, WithNodes = F)

plot3d(SSN_Q1Q2Q3Q4,
       soma = FALSE, lwd = 1, add = T, alpha = 0.25, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)

plot3d(outline,
       add = T, alpha = 0.05, col = "grey50"
)

# plot mitochondria

pos_ves_syn <- mito_vesicle_info |>
  filter(celltype == "SSN") |>
  filter(mito_type == "vesicles_syn") |>
  select(x, y, z)

plot3d(pos_ves_syn, 
       add = TRUE,
       col = "green3",
       size = 0.6, 
       alpha = 0.5,
       point_antialias = TRUE,
       type = "s"
)

pos_ves_no_syn <- mito_vesicle_info |>
  filter(celltype == "SSN") |>
  filter(mito_type != "vesicles_syn") |>
  select(x, y, z)

plot3d(pos_ves_no_syn, 
       add = TRUE,
       col = "black",
       size = 0.6, 
       alpha = 0.5,
       point_antialias = TRUE,
       type = "s"
)

#aboral view
aboral()

par3d(zoom=0.61)


#move to next panel in rgl window
next3d(clear=F)

#plot sagittal view
plot3d(SSN_Q1Q2,
       soma = FALSE, lwd = 1, add = T, alpha = 0.25, col = Okabe_Ito[6],
       WithConnectors = F, WithNodes = F)

plot3d(SSN_Q3Q4,
       soma = FALSE, lwd = 1, add = T, alpha = 0.25, col = Okabe_Ito[7],
       WithConnectors = F, WithNodes = F)

plot3d(SSN_Q1Q2Q3Q4,
       soma = FALSE, lwd = 1, add = T, alpha = 0.25, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)

plot3d(outline,
       add = T, alpha = 0.05, col = "grey50"
)

# plot mitochondria

pos_ves_syn <- mito_vesicle_info |>
  filter(celltype == "SSN") |>
  filter(mito_type == "vesicles_syn") |>
  select(x, y, z)

plot3d(pos_ves_syn, 
       add = TRUE,
       col = "green3",
       size = 0.6, 
       alpha = 0.5,
       point_antialias = TRUE,
       type = "s"
)

pos_ves_no_syn <- mito_vesicle_info |>
  filter(celltype == "SSN") |>
  filter(mito_type != "vesicles_syn") |>
  select(x, y, z)

plot3d(pos_ves_no_syn, 
       add = TRUE,
       col = "black",
       size = 0.6, 
       alpha = 0.5,
       point_antialias = TRUE,
       type = "s"
)

#sagittal view
sagittal()

par3d(zoom=0.61)


#move to next panel in rgl window
next3d(clear=F)

#plot tentacular view
plot3d(SSN_Q1Q2,
       soma = FALSE, lwd = 1, add = T, alpha = 0.25, col = Okabe_Ito[6],
       WithConnectors = F, WithNodes = F)

plot3d(SSN_Q3Q4,
       soma = FALSE, lwd = 1, add = T, alpha = 0.25, col = Okabe_Ito[7],
       WithConnectors = F, WithNodes = F)

plot3d(SSN_Q1Q2Q3Q4,
       soma = FALSE, lwd = 1, add = T, alpha = 0.25, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)

plot3d(outline,
       add = T, alpha = 0.05, col = "grey50"
)

# plot mitochondria

pos_ves_syn <- mito_vesicle_info |>
  filter(celltype == "SSN") |>
  filter(mito_type == "vesicles_syn") |>
  select(x, y, z)

plot3d(pos_ves_syn, 
       add = TRUE,
       col = "green3",
       size = 0.6, 
       alpha = 0.5,
       point_antialias = TRUE,
       type = "s"
)

pos_ves_no_syn <- mito_vesicle_info |>
  filter(celltype == "SSN") |>
  filter(mito_type != "vesicles_syn") |>
  select(x, y, z)

plot3d(pos_ves_no_syn, 
       add = TRUE,
       col = "black",
       size = 0.6, 
       alpha = 0.5,
       point_antialias = TRUE,
       type = "s"
)

#tentacular view
tentacular()

par3d(zoom=0.61)


#make a snapshot to the working directory
rgl.snapshot("manuscript/pictures/mito_pos_SSN.png")

close3d()



# retrieve connectors ----------------------------------------------------------

#conn_syn <- connectors(SSN)
#presyn_syn <- subset(conn_syn, prepost == 0)
#postsyn_syn <- subset(conn_syn, prepost == 1)


# don't use many_to_many connectors, because getting cords looks a bit complicated
#syn_Q1Q2Q3Q4_to_Q1Q2 <- catmaid_fetch(path = "/35/connector/list/many_to_many",
#                                      body = list(skids1=SSN_Q1Q2Q3Q4_skid,
#                                                  skids2=SSN_Q1Q2_skid,
#                                                  relation="presynaptic_to"))
#syn_Q1Q2Q3Q4_to_Q3Q4 <- catmaid_fetch(path = "/35/connector/list/many_to_many",
#                                      body = list(skids1=SSN_Q1Q2Q3Q4_skid,
#                                                  skids2=SSN_Q3Q4_skid,
#                                                  relation="presynaptic_to"))

SSN_Q1Q2_skid <- SSN_Q1Q2$skid
SSN_Q3Q4_skid <- SSN_Q3Q4$skid
SSN_Q1Q2Q3Q4_skid <- SSN_Q1Q2Q3Q4$skid


# "presynapse (output)"

stats_synapse <- read.csv("analysis/data/stats_synapse.csv")

syn_big_out_connectors <- stats_synapse %>%
  filter(skid==SSN_Q1Q2Q3Q4_skid) %>%
  filter(prepost==0) %>%
  select(connector_id) %>%
  pull()

syn_small_out_connectors  <- stats_synapse %>%
  filter(skid==SSN_Q1Q2_skid | skid==SSN_Q3Q4_skid) %>%
  filter(prepost==0) %>%
  select(connector_id) %>%
  pull()


# "postsynapse (input)"

syn_big_to_small <- stats_synapse %>%
  filter(skid==SSN_Q1Q2_skid | skid==SSN_Q3Q4_skid) %>%
  filter(connector_id %in% syn_big_out_connectors) %>%
  filter(prepost==1)

syn_big_to_big <- stats_synapse %>%
  filter(skid==SSN_Q1Q2Q3Q4_skid) %>%
  filter(connector_id %in% syn_big_out_connectors) %>%
  filter(prepost==1)

syn_small_to_big <- stats_synapse %>%
  filter(skid==SSN_Q1Q2Q3Q4_skid) %>%
  filter(connector_id %in% syn_small_out_connectors) %>%
  filter(prepost==1)

syn_small_to_small <- stats_synapse %>%
  filter(skid==SSN_Q1Q2_skid | skid==SSN_Q3Q4_skid) %>%
  filter(connector_id %in% syn_small_out_connectors) %>%
  filter(prepost==1)


# 3d plot SSN synapses ---------------------------------------------------------


close3d()
# 3d plotting of cells
nopen3d() 
mfrow3d(1, 3)  
#define the size of the rgl window, the view and zoom
par3d(windowRect = c(0, 0, 1200, 350))
#par3d(windowRect = c(0, 0, 2400, 700))


#plot aboral view
plot3d(SSN_Q1Q2,
       soma = FALSE, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[6],
       WithConnectors = F, WithNodes = F)

plot3d(SSN_Q3Q4,
       soma = FALSE, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[7],
       WithConnectors = F, WithNodes = F)

plot3d(SSN_Q1Q2Q3Q4,
       soma = FALSE, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)

plot3d(outline,
       add = T, alpha = 0.05, col = "grey50"
)


# plot synapses

pos_big_to_small <- syn_big_to_small %>%
  select(x, y, z)
plot3d(pos_big_to_small,
       size = 1, alpha = 0.75, col = "magenta2", 
       add = TRUE,
       point_antialias = TRUE,
       type = "s"
)

pos_small_to_big <- syn_small_to_big %>%
  select(x, y, z)
plot3d(pos_small_to_big,
       size = 1, alpha = 0.75, col = "cyan2", 
       add = TRUE,
       point_antialias = TRUE,
       type = "s"
)


pos_big_to_big <- syn_big_to_big %>%
  select(x, y, z)
plot3d(pos_big_to_big,
       size = 1, alpha = 0.75, col = "#4477AA", 
       add = TRUE,
       point_antialias = TRUE,
       type = "s"
)

#pos_small_to_small <- syn_small_to_small %>%
#  select(x, y, z)
#plot3d(pos_small_to_small,
#       size = 0.6, alpha = 0.5, col = "#4477AA", 
#       add = TRUE,
#       point_antialias = TRUE,
#       type = "s"
#)


#aboral view
nview3d("anterior", extramat = rotationMatrix(1.05, 250, -200, 1000))
#rgl.snapshot("manuscript/pictures/balancer_aboral_view.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)




#plot lateral view of Sagittal plane
plot3d(SSN_Q1Q2,
       soma = FALSE, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[6],
       WithConnectors = F, WithNodes = F)

plot3d(SSN_Q3Q4,
       soma = FALSE, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[7],
       WithConnectors = F, WithNodes = F)

plot3d(SSN_Q1Q2Q3Q4,
       soma = FALSE, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)

plot3d(outline,
       add = T, alpha = 0.05, col = "grey50"
)

# plot synapses

pos_big_to_small <- syn_big_to_small %>%
  select(x, y, z)
plot3d(pos_big_to_small,
       size = 1, alpha = 0.75, col = "magenta2", 
       add = TRUE,
       point_antialias = TRUE,
       type = "s"
)

pos_small_to_big <- syn_small_to_big %>%
  select(x, y, z)
plot3d(pos_small_to_big,
       size = 1, alpha = 0.75, col = "cyan2", 
       add = TRUE,
       point_antialias = TRUE,
       type = "s"
)


pos_big_to_big <- syn_big_to_big %>%
  select(x, y, z)
plot3d(pos_big_to_big,
       size = 1, alpha = 0.75, col = "#4477AA", 
       add = TRUE,
       point_antialias = TRUE,
       type = "s"
)



nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
#rgl.snapshot("manuscript/pictures/balancer_sagittal_plane.png")
par3d(zoom=0.61)

#move to next panel in rgl window
next3d(clear=F)



#plot lateral view of Tentacular plane
plot3d(SSN_Q1Q2,
       soma = FALSE, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[6],
       WithConnectors = F, WithNodes = F)

plot3d(SSN_Q3Q4,
       soma = FALSE, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[7],
       WithConnectors = F, WithNodes = F)

plot3d(SSN_Q1Q2Q3Q4,
       soma = FALSE, lwd = 1, add = T, alpha = 0.5, col = Okabe_Ito[5],
       WithConnectors = F, WithNodes = F)

plot3d(outline,
       add = T, alpha = 0.05, col = "grey50"
)

# plot synapses

pos_big_to_small <- syn_big_to_small %>%
  select(x, y, z)
plot3d(pos_big_to_small,
       size = 1, alpha = 0.75, col = "magenta2", 
       add = TRUE,
       point_antialias = TRUE,
       type = "s"
)

pos_small_to_big <- syn_small_to_big %>%
  select(x, y, z)
plot3d(pos_small_to_big,
       size = 1, alpha = 0.75, col = "cyan2", 
       add = TRUE,
       point_antialias = TRUE,
       type = "s"
)


pos_big_to_big <- syn_big_to_big %>%
  select(x, y, z)
plot3d(pos_big_to_big,
       size = 1, alpha = 0.75, col = "#4477AA", 
       add = TRUE,
       point_antialias = TRUE,
       type = "s"
)


nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))
#rgl.snapshot("manuscript/pictures/balancer_tentacular_plane.png")
par3d(zoom=0.61)


#make a snapshot to the working directory
rgl.snapshot("manuscript/pictures/SSN_prepost_synapse.png")

close3d()



# assemble figure -------------------------------------------------------------

panel_SSN_Q1234 <- ggdraw() + draw_image(readPNG("manuscript/pictures/SSN_Q1234.png")) +
  draw_label("Aboral nerve net (ANN) Q1-4", x = 0.5, y = 1.05, color = Okabe_Ito[5], size = 10, hjust = 0.5) +
  draw_label("aboral view", x = 0.18, y = 0.95, color = "black", size = 8, hjust = 0.5) +
  draw_label("lateral view", x = 0.67, y = 0.95, color = "black", size = 8, hjust = 0.5) +
  draw_label("Q1", x = 0.3, y = 0.9, color = "black", size = 9, hjust = 0.5, alpha = 0.5) +
  draw_label("Q2", x = 0.06, y = 0.79, color = "black", size = 9, hjust = 0.5, alpha = 0.5) +
  draw_label("Q3", x = 0.03, y = 0.08, color = "black", size = 9, hjust = 0.5, alpha = 0.5) +
  draw_label("Q4", x = 0.275, y = 0.11, color = "black", size = 9, hjust = 0.5, alpha = 0.5) +
  draw_label("sagittal plane", x = 0.5, y = 0.86, color = "black", size = 7.5, hjust = 0.5) +
  draw_label("tentacular plane", x = 0.83, y = 0.86, color = "black", size = 7.5, hjust = 0.5) +
  draw_label("A", x = 0.67, y = 0.85, size = 7, color = "black", hjust = 0.5) +
  draw_label("O", x = 0.67, y = 0.65, size = 7, color = "black", hjust = 0.5) +
  draw_line(x = c(0.67, 0.67), y = c(0.69, 0.81), color = "black", linewidth = 0.65) +
  draw_label(expression(paste("25 ", mu, "m")), x = 0.95, y = 0.1, color = "black", size = 7, hjust = 0.5) +
  draw_line(x = c(0.91, 0.99), y = c(0.05, 0.05), color = "black", linewidth = 0.7)

panel_SSN_Q12_Q34 <- ggdraw() + draw_image(readPNG("manuscript/pictures/SSN_Q12_Q34.png")) +
#  draw_label("ANN Q1Q2, ANN Q3Q4", x = 0.5, y = 1, color = "gray", size = 10, hjust = 0.5) +
  draw_label("ANN Q1Q2", x = 0.347, y = 1, color = Okabe_Ito[6], size = 10, hjust = 0) +
  draw_label(",", x = 0.492, y = 1, color = "black", size = 10, hjust = 0) +
  draw_label("ANN Q3Q4", x = 0.507, y = 1, color = Okabe_Ito[7], size = 10, hjust = 0) +
  draw_label("Q1", x = 0.3, y = 0.9, color = "black", size = 9, hjust = 0.5, alpha = 0.5) +
  draw_label("Q2", x = 0.06, y = 0.79, color = "black", size = 9, hjust = 0.5, alpha = 0.5) +
  draw_label("Q3", x = 0.03, y = 0.08, color = "black", size = 9, hjust = 0.5, alpha = 0.5) +
  draw_label("Q4", x = 0.275, y = 0.11, color = "black", size = 9, hjust = 0.5, alpha = 0.5)

panel_mito_bar_graph <- ggdraw() + draw_image(readPNG("manuscript/pictures/plot_mito_stats.png"))

circle_g <- ggdraw() + draw_image(readPNG("manuscript/pictures/synapse_green.png"))
circle_gy <- ggdraw() + draw_image(readPNG("manuscript/pictures/synapse_gray.png"))

panel_mito_pos <- ggdraw() + draw_image(readPNG("manuscript/pictures/mito_pos_SSN.png")) +
  draw_plot(circle_g, x = 0.06, y = 0.955, width = 0.05, height = 0.05) +
  draw_label("Mitochondria with synapses", x = 0.1, y = 0.98, color = "green4", size = 10.5, hjust = 0) +
  draw_plot(circle_gy, x = 0.06, y = 0.875, width = 0.05, height = 0.05) +
  draw_label("Not forming synapses", x = 0.1, y = 0.9, color = "black", size = 10.5, hjust = 0, alpha = 0.7) +
  draw_label("ANN Q1-4", x = 0.85, y = 1, color = Okabe_Ito[5], size = 8, hjust = 0) +
  draw_label("ANN Q1Q2", x = 0.85, y = 0.94, color = Okabe_Ito[6], size = 8, hjust = 0) +
  draw_label("ANN Q3Q4", x = 0.85, y = 0.88, color = Okabe_Ito[7], size = 8, hjust = 0)

panel_EM_SSN <- ggdraw() + draw_image(readPNG("manuscript/pictures/SSN_EM_schematic.png"))

circle_m <- ggdraw() + draw_image(readPNG("manuscript/pictures/synapse_magenta.png"))
circle_lb <- ggdraw() + draw_image(readPNG("manuscript/pictures/synapse_lightblue.png"))
circle_db <- ggdraw() + draw_image(readPNG("manuscript/pictures/synapse_darkblue.png"))

panel_SSN_prepost_synapse <- ggdraw() + draw_image(readPNG("manuscript/pictures/SSN_prepost_synapse.png")) +
  #  draw_label("ANN Q1-4 → ANN Q1Q2 or Q3Q4", x = 0.06, y = 1.06, color = "gray", size = 8.5, hjust = 0) +
  draw_plot(circle_m, x = 0.02, y = 1.035, width = 0.05, height = 0.05) +
  draw_label("ANN Q1-4", x = 0.06, y = 1.06, color = Okabe_Ito[5], size = 8.5, hjust = 0) +
  draw_label(" → ", x = 0.175, y = 1.06, color = "black", size = 8.5, hjust = 0) +
  draw_label("ANN Q1Q2", x = 0.21, y = 1.06, color = Okabe_Ito[6], size = 8.5, hjust = 0) +
  draw_label(" or ", x = 0.335, y = 1.06, color = "black", size = 8.5, hjust = 0) +
  draw_label("Q3Q4", x = 0.37, y = 1.06, color = Okabe_Ito[7], size = 8.5, hjust = 0) +
  #  draw_label("ANN Q1Q2 or Q3Q4 → ANN Q1-4", x = 0.06, y = 0.96, color = "gray", size = 8.5, hjust = 0) +
  draw_plot(circle_lb, x = 0.02, y = 0.935, width = 0.05, height = 0.05) +
  draw_label("ANN Q1Q2", x = 0.06, y = 0.96, color = Okabe_Ito[6], size = 8.5, hjust = 0) +
  draw_label(" or ", x = 0.185, y = 0.96, color = "black", size = 8.5, hjust = 0) +
  draw_label("Q3Q4", x = 0.22, y = 0.96, color = Okabe_Ito[7], size = 8.5, hjust = 0) +
  draw_label(" → ", x = 0.285, y = 0.96, color = "black", size = 8.5, hjust = 0) +
  draw_label("ANN Q1-4", x = 0.32, y = 0.96, color = Okabe_Ito[5], size = 8.5, hjust = 0) +
  #  draw_label("ANN Q1-4 → ANN Q1-4", x = 0.56, y = 1.06, color = "gray", size = 8.5, hjust = 0) +
  draw_plot(circle_db, x = 0.52, y = 1.035, width = 0.05, height = 0.05) +
  draw_label("ANN Q1-4", x = 0.56, y = 1.06, color = Okabe_Ito[5], size = 8.5, hjust = 0) +
  draw_label(" → ", x = 0.675, y = 1.06, color = "black", size = 8.5, hjust = 0) +
  draw_label("ANN Q1-4", x = 0.71, y = 1.06, color = Okabe_Ito[5], size = 8.5, hjust = 0)

#panel_SSN_graph <- ggdraw() + draw_image(readPNG("manuscript/pictures/SSN_graph.png"))


layout <- "
#########
AAAA#BBBB
#########
CCCC#DDDD
#########
EEEE#FFFF
#########
"

Figure2 <- panel_SSN_Q1234 + panel_SSN_Q12_Q34 +
  panel_mito_bar_graph + panel_mito_pos + 
  panel_EM_SSN + panel_SSN_prepost_synapse +
  plot_layout(design = layout,
              heights = c(0.05, 1, 0.15, 1.3, 0.15, 1, 0.05),
              widths = c(1, 1, 1, 1, 0.2, 1, 1, 1, 1)) + 
  patchwork::plot_annotation(tag_levels = "A") &  
  ggplot2::theme(plot.tag = element_text(size = 12, 
                                         face='plain', color='black'))


ggsave("manuscript/figures/Figure2.png", limitsize = FALSE, 
       units = c("px"), Figure2, width = 3000, height = 1800, bg='white')  


ggsave("manuscript/figures/Figure2.pdf", limitsize = FALSE, 
       units = c("px"), Figure2, width = 3000, height = 1800) 







