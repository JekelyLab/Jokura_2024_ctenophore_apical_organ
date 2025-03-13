# Code to generate Figure 2 of the Jokura et al 2024 Ctenophore apical organ connectome paper

# source packages and functions ------------------------------------------------
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

Q1 <- read_smooth_neuron("Q1")
Q2 <- read_smooth_neuron("Q2")
Q3 <- read_smooth_neuron("Q3")
Q4 <- read_smooth_neuron("Q4")

SSN_Q1Q2 <- read_smooth_neuron("SSN_Q1Q2")[[1]]
SSN_Q3Q4 <- read_smooth_neuron("SSN_Q3Q4")[[1]]
SSN_Q1Q2Q3Q4 <- read_smooth_neuron("SSN_Q1Q2Q3Q4")[[1]]



with_soma <- read_smooth_neuron("with_soma")

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

plot_multinucleated_cell(SSN_Q1Q2Q3Q4,
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





# load of mitochondrial location information-------------------------------

mito_done <- read.neurons.catmaid("mitochondria done", pid = 35)
mito_vesicle_info <- read.csv("analysis/data/mito_vesicle_info.csv")
mito_means_tidy <- read.csv("analysis/data/mito_means_tidy.csv")

# graph plot mito vesicles by celltype -----------------------------------------------

mito_means_tidy <- mito_means_tidy %>%
  mutate(
    characteristic = ifelse(characteristic == "mean_vesicles_syn", 
                            "mean_vesicles_syn", 
                            "mean_vesicles_no_syn"),
    celltype = ifelse(celltype %in% c("Cgroove-sag", "Cgroove-tag"), 
                      "ciliated_groove", 
                      celltype)
  ) %>%
  mutate(celltype = factor(celltype, levels = c(
    "balancer", "bridge", "bristle", "ciliated_groove", "dense_vesicle", "dome", 
    "intra-multi-ciliated", "lamellate", "lithocyte", "plumose", "SSN", 
    "epithelial_floor", "monociliated", "biciliated", "multiciliated", "nonciliated"
  )))

plot_mito_stats <- ggplot(mito_means_tidy, aes(
  fill = factor(characteristic, levels = c("mean_vesicles_syn", "mean_vesicles_no_syn")),
  y = value,
  x = celltype
)) +
  geom_bar(position = "stack", stat = "identity") +
  theme_bw() +
  scale_fill_manual(
    values = c("mean_vesicles_syn" = "green4", "mean_vesicles_no_syn" = "darkgrey"),
    labels = c("mean_vesicles_syn" = "mitochondria with synapses", 
               "mean_vesicles_no_syn" = "mitochondria not forming synapses")
  ) +
  theme(
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(angle = -70, vjust = 0.5, hjust = 0, size = 9, margin = margin(t = -7)),
    axis.text.y = element_text(size = 10), 
    axis.title = element_text(size = 10),
    legend.title = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.3, 0.75),
    text = element_text(size = 11)
  ) +
  ylab("Average number of mitochondria per cell") +
  xlab("Cell types")

plot_mito_stats

ggsave(
  filename = "manuscript/pictures/mito_in_syn.png",
  plot = plot_mito_stats,
  width = 1600,
  height = 1100,
  units = "px",
  dpi = 300
)


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




# outputs from SNN -------------------------------------------------------------
# bar graph: y-axis - number of inputs from SSNs, x-axis - cell type
# exclude monoC, biC, multiC, and nonC.
 
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


SSN_downstream <- SSN_downstream %>%
  mutate(celltype = get_celltype_annot_for_skid(skid))



#stats_master <- read.csv("analysis/data/stats_master.csv") 
#SSN_downstream <- SSN_downstream %>%
#  left_join(stats_master %>% select(skid, celltype), by = "skid")



# use get_celltype_annot_for_skid function from cell_statistics_master.R to get celltypes
# or get celltype from stats_master.csv
# add celltype as column to SSN_downstream

# filter out monoC, etc
#SNN_downstream %>%
#  filter(celltype != "monoC" | celltype != "biC")


SSN_downstream_filtered <- SSN_downstream %>%
  filter(!celltype %in% c("monoC", "biC", "multiC", "nonC"))




ggplot(SSN_downstream_filtered, aes(x = celltype)) +
  geom_bar(fill = "steelblue") +
  theme_minimal() +
  labs(
    title = "Number of Inputs from SSNs by Cell Type",
    x = "Cell Type",
    y = "Number of Inputs"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))









# assemble figure -------------------------------------------------------------

panel_SSN_Q1234 <- ggdraw() + draw_image(readPNG("manuscript/pictures/SSN_Q1234.png")) +
  draw_label("Q1", x = 0.3, y = 0.9, color = "black", size = 9, hjust = 0.5, alpha = 0.5) +
  draw_label("Q2", x = 0.06, y = 0.79, color = "black", size = 9, hjust = 0.5, alpha = 0.5) +
  draw_label("Q3", x = 0.03, y = 0.08, color = "black", size = 9, hjust = 0.5, alpha = 0.5) +
  draw_label("Q4", x = 0.275, y = 0.11, color = "black", size = 9, hjust = 0.5, alpha = 0.5) +
  draw_label("Syncytial synaptic neuron (SSN) Q1-4", x = 0.05, y = 1, color = Okabe_Ito[5], size = 10, hjust = 0)

panel_SSN_Q12_Q34 <- ggdraw() + draw_image(readPNG("manuscript/pictures/SSN_Q12_Q34.png")) +
  draw_label("Q1", x = 0.3, y = 0.9, color = "black", size = 9, hjust = 0.5, alpha = 0.5) +
  draw_label("Q2", x = 0.06, y = 0.79, color = "black", size = 9, hjust = 0.5, alpha = 0.5) +
  draw_label("Q3", x = 0.03, y = 0.08, color = "black", size = 9, hjust = 0.5, alpha = 0.5) +
  draw_label("Q4", x = 0.275, y = 0.11, color = "black", size = 9, hjust = 0.5, alpha = 0.5) +
  draw_label("SSN Q1Q2", x = 0.05, y = 1, color = Okabe_Ito[6], size = 10, hjust = 0) +
  draw_label("SSN Q3Q4", x = 0.225, y = 1, color = Okabe_Ito[7], size = 10, hjust = 0)

panel_mito_bar_graph <- ggdraw() + draw_image(readPNG("manuscript/pictures/mito_in_syn.png"))

panel_mito_pos <- ggdraw() + draw_image(readPNG("manuscript/pictures/mito_pos_SSN.png")) +
  draw_label("mitochondria with synapses", x = 0.05, y = 0.98, color = "green4", size = 10.5, hjust = 0) +
  draw_label("mitochondria not forming synapses", x = 0.05, y = 0.9, color = "black", size = 10.5, hjust = 0, alpha = 0.8) +
  draw_label("SSN Q1-4", x = 0.85, y = 1, color = Okabe_Ito[5], size = 8, hjust = 0) +
  draw_label("SSN Q1Q2", x = 0.85, y = 0.94, color = Okabe_Ito[6], size = 8, hjust = 0) +
  draw_label("SSN Q3Q4", x = 0.85, y = 0.88, color = Okabe_Ito[7], size = 8, hjust = 0)

panel_EM_SSN <- ggdraw() + draw_image(readPNG("manuscript/pictures/SSN_EM_schematic.png"))

panel_SSN_prepost_synapse <- ggdraw() + draw_image(readPNG("manuscript/pictures/SSN_prepost_synapse.png")) +
  draw_label("Synapses from SSN Q1-4 to SSN Q1Q2 or Q3Q4", x = 0.05, y = 1.06, color = "#ff00ff", size = 10, hjust = 0) +
  draw_label("Synapses from SSN Q1Q2 or Q3Q4 to SSN Q1-4", x = 0.05, y = 0.96, color = "#00c9ff", size = 10, hjust = 0)


panel_SSN_graph <- ggdraw() + draw_image(readPNG("manuscript/pictures/SSN_graph.png"))


layout <- "
AAAABBBB
########
CCCDDDDG
########
EEEEFFFF
"

Figure2 <- panel_SSN_Q1234 + panel_SSN_Q12_Q34 +
  panel_mito_bar_graph + panel_mito_pos + 
  panel_EM_SSN + panel_SSN_prepost_synapse + 
  panel_SSN_graph +
  plot_layout(design = layout,
              heights = c(1, 0.1, 1.3, 0.15, 1),
              widths = c(1, 1, 1, 1, 1, 1, 1, 1)) + 
  patchwork::plot_annotation(tag_levels = "A") &  
  ggplot2::theme(plot.tag = element_text(size = 12, 
                                         face='plain', color='black'))


ggsave("manuscript/figures/Figure2.png", limitsize = FALSE, 
       units = c("px"), Figure2, width = 3000, height = 1700, bg='white')  


ggsave("manuscript/figures/Figure2.pdf", limitsize = FALSE, 
       units = c("px"), Figure2, width = 3000, height = 1600) 







