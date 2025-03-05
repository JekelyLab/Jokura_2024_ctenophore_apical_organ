source("analysis/scripts/packages_and_functions.R")

# calculate average number of post-synaptic sites per synapse  -----

pre_connectors <- catmaid_fetch(
  path = "35/connectors/",
  body = list(
    relation_type = "presynaptic_to",
    with_partners = "false"
  )
)

n_pre_post <- function(connector_id) {
  connector_info <- catmaid_fetch(path = paste("35/connectors/", connector_id, sep = ""))
  # pre should always be 1, but maybe it's good to check if there is any weird stuff going on
  pre <- connector_info$partners %like% "presynaptic_to"
  n_pre <- table(pre)["TRUE"][[1]]
  post <- connector_info$partners %like% "postsynaptic_to"
  n_post <- table(post)["TRUE"][[1]]
  return(n_post)
}

pre_connector_IDs <- lapply(pre_connectors$connectors, "[[", 1)

n_post <- lapply(pre_connector_IDs, n_pre_post)

summary(as.numeric(n_post))

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

#for (mito_done_neu in mito_done) {
#  print(mito_done_neu[1])
#  get_mito_pos(mito_done_neu)
#}

write.csv(mito_vesicle_info, "analysis/data/mito_vesicle_info.csv")
mito_vesicle_info <- read.csv("analysis/data/mito_vesicle_info.csv")

# create tibble for graph ---------------------------------

mito_stats <- mito_vesicle_info %>%
  group_by(skid, mito_type) %>%
  count() %>%
  pivot_wider(names_from = mito_type, values_from = n, values_fill = 0) %>%
  left_join(select(mito_vesicle_info, skid, celltype), by = "skid") %>%
  unique() %>%
  select(celltype, skid, vesicles_syn, vesicles_no_syn, vesicles_unc, vesicles_none) %>%
  arrange(desc(vesicles_syn), desc(vesicles_no_syn), desc(vesicles_none)) %>%
  mutate(vesicles_total = rowSums(across(where(is.numeric))))

ves_syn_tbl <- mito_stats %>%
  group_by(celltype) %>%
  summarize(mean_vesicles_syn = mean(vesicles_syn)) %>%
  arrange(celltype) %>%
  select(mean_vesicles_syn)
ves_no_syn_tbl <- mito_stats %>%
  group_by(celltype) %>%
  summarize(mean_vesicles_no_syn = mean(vesicles_no_syn)) %>%
  arrange(celltype) %>%
  select(mean_vesicles_no_syn)
ves_unc_tbl <- mito_stats %>%
  group_by(celltype) %>%
  summarize(mean_vesicles_unclear = mean(vesicles_unc)) %>%
  arrange(celltype) %>%
  select(mean_vesicles_unclear)
ves_none_tbl <- mito_stats %>%
  group_by(celltype) %>%
  summarize(mean_vesicles_none = mean(vesicles_none)) %>%
  arrange(celltype) %>%
  select(mean_vesicles_none)
celltype_numbers <- mito_stats %>%
  group_by(celltype) %>%
  tally() %>%
  arrange(celltype)
mito_means <- bind_cols(
  celltype_numbers,
  ves_syn_tbl, ves_no_syn_tbl,
  ves_unc_tbl, ves_none_tbl
)


mito_means_tidy <- mito_means %>%
  # select(total) %>%
  pivot_longer(
    cols = c("mean_vesicles_syn", "mean_vesicles_no_syn", "mean_vesicles_unclear", "mean_vesicles_none"),
    names_to = "characteristic",
    values_to = "value"
  )

write.csv(mito_means_tidy, "analysis/data/mito_means_tidy.csv")
mito_means_tidy <- read.csv("analysis/data/mito_means_tidy.csv")

# plot mito vesicles by celltype -----------------------------------------------

plot_mito_stats <- ggplot(mito_means_tidy, aes(
  fill = factor(characteristic,
    levels = c("mean_vesicles_syn", "mean_vesicles_no_syn", "mean_vesicles_unclear", "mean_vesicles_none")
  ),
  y = value, x = interaction(celltype, n)
)) +
  geom_bar(position = "stack", stat = "identity") +
  theme_bw() +
  scale_fill_manual(values = c("#EE6677", "#DDAA33", "#4477AA", "darkgrey")) +
  theme(
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.title = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.17, .85),
    text = element_text(family = "sans", size = 12)
  )

plot_mito_stats

# plot mitochondria positions in SSN and bridge --------------------------------
# I'm not sure it makes sense to do it for other celltypes because not all cells
# for them had annotated mito
for (cell_type in c("SSN", "bridge")) {
  nopen3d()
  # plot3d(outline,
  #  WithNodes = F, add = T, alpha = 0.07, col = "#E2E2E2"
  # )
  skids <- mito_stats %>%
    filter(celltype == cell_type) %>%
    select(skid) %>%
    pull()

  for (skid in skids) {
    cell <- mito_done[[as.character(skid)]] %>% smooth_neuron(sigma = 1000)
    plot_multinucleated_cell(cell, color = "grey", lwd = 2, alpha = 0.6)
  }

  pos_ves_syn <- mito_vesicle_info |>
    filter(celltype == cell_type) |>
    filter(mito_type == "vesicles_syn") |>
    select(x, y, z)
  plot3d(pos_ves_syn,
    add = TRUE,
    col = "#EE6677",
    size = 5,
    alpha = 1
  )

  pos_ves_no_syn <- mito_vesicle_info |>
    filter(celltype == cell_type) |>
    filter(mito_type == "vesicles_no_syn") |>
    select(x, y, z)
  plot3d(pos_ves_no_syn,
    add = TRUE,
    col = "#DDAA33",
    size = 5,
    alpha = 1
  )

  pos_ves_unc <- mito_vesicle_info |>
    filter(celltype == cell_type) |>
    filter(mito_type == "vesicles_unc") |>
    select(x, y, z)
  plot3d(pos_ves_unc,
    add = TRUE,
    col = "#4477AA",
    size = 5,
    alpha = 1
  )

  pos_ves_no <- mito_vesicle_info |>
    filter(celltype == cell_type) |>
    filter(mito_type == "vesicles_none") |>
    select(x, y, z)
  plot3d(pos_ves_no,
    add = TRUE,
    col = "black",
    size = 5,
    alpha = 1
  )
  par3d(zoom=0.50)
  # anterior
  nview3d("anterior", extramat = rotationMatrix(1.05, 250, -200, 1000))
  par3d(windowRect = c(0, 0, 400, 800))
  paste("manuscript/pictures/mito_pos_", cell_type, "_aboral.png", sep = "") %>% rgl.snapshot()
  # sagittal
  nview3d("left", extramat = rotationMatrix(300, 4200, 1800, 800))
  par3d(windowRect = c(0, 0, 450, 450))
  paste("manuscript/pictures/mito_pos_", cell_type, "_sagittal.png", sep = "") %>% rgl.snapshot()
  # tentacular
  nview3d("left", extramat = rotationMatrix(-1.7, 190, -120, -140))
  par3d(windowRect = c(0, 0, 800, 450))
  paste("manuscript/pictures/mito_pos_", cell_type, "_tentacular.png", sep = "") %>% rgl.snapshot()
  close3d()
}



# TODO: automate renaming. Mybe make the crop funtion return a filename?
if (!file.exists("manuscript/pictures/Figure_mito_syn_ves_syn.tiff")) {
  crop_substack("Figure_mito_syn_ves_syn", 700, 700, 0, 0, "manuscript/pictures", 35, 28)
}
if (!file.exists("manuscript/pictures/Figure_mito_syn_ves_no_syn.tiff")) {
  crop_substack("Figure_mito_syn_ves_no_syn", 700, 700, 0, 0, "manuscript/pictures", 35, 28)
}
if (!file.exists("manuscript/pictures/Figure_mito_syn_ves_unc.tiff")) {
  crop_substack("Figure_mito_syn_ves_unc", 700, 700, 0, 0, "manuscript/pictures", 35, 28)
}
if (!file.exists("manuscript/pictures/Figure_mito_syn_ves_none.tiff")) {
  crop_substack("Figure_mito_syn_ves_none", 700, 700, 0, 0, "manuscript/pictures", 35, 28)
}


EM_ves_syn <- magick::image_read("manuscript/pictures/Figure_mito_syn_ves_syn.tiff")
EM_ves_no_syn <- magick::image_read("manuscript/pictures/Figure_mito_syn_ves_no_syn.tiff")
EM_ves_unc <- magick::image_read("manuscript/pictures/Figure_mito_syn_ves_unc.tiff")
EM_ves_none <- magick::image_read("manuscript/pictures/Figure_mito_syn_ves_none.tiff")

panel_EM_ves_syn <- ggdraw() + draw_image(EM_ves_syn)
panel_EM_ves_no_syn <- ggdraw() + draw_image(EM_ves_no_syn)
panel_EM_ves_unc <- ggdraw() + draw_image(EM_ves_unc)
panel_EM_ves_none <- ggdraw() + draw_image(EM_ves_none)

mito_pos_SSN_anterior <- magick::image_read("manuscript/pictures/mito_pos_SSN_anterior.png")
mito_pos_SSN_left <- magick::image_read("manuscript/pictures/mito_pos_SSN_left.png")
mito_pos_bridge_anterior <- magick::image_read("manuscript/pictures/mito_pos_bridge_anterior.png")
mito_pos_bridge_left <- magick::image_read("manuscript/pictures/mito_pos_bridge_left.png")

panel_mito_pos_SSN_anterior <- ggdraw() + draw_image(mito_pos_SSN_anterior)
panel_mito_pos_SSN_left <- ggdraw() + draw_image(mito_pos_SSN_left)
panel_mito_pos_bridge_anterior <- ggdraw() + draw_image(mito_pos_bridge_anterior)
panel_mito_pos_bridge_left <- ggdraw() + draw_image(mito_pos_bridge_left)


layout <- "AABBCC
           AADDEE
           AAFGHI"

figure_mito_in_syn <- plot_mito_stats +
  panel_mito_pos_SSN_anterior + panel_mito_pos_SSN_left + panel_mito_pos_bridge_anterior +  panel_mito_pos_bridge_left +
  panel_EM_ves_syn + panel_EM_ves_no_syn + panel_EM_ves_unc + panel_EM_ves_none +
  plot_layout(design = layout) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 12, face = "plain"))
figure_mito_in_syn
