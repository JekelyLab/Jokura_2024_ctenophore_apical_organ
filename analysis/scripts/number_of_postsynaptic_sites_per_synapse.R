library(catmaid)
library(data.table)
library(tidyverse)

conn <- source("~/R/conn.R")

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


# create tibble for graph ---------------------------------

mito_stats <- mito_vesicle_info %>%
  group_by(skid, mito_type) %>%
  count() %>%
  pivot_wider(names_from = mito_type, values_from = n, values_fill = 0) %>%
  left_join(select(mito_stats2, skid, celltype), by = "skid") %>%
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
  add_count(celltype) %>%
  ungroup() %>%
  select(celltype, n) %>%
  unique() %>%
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


# plot mito vesicles by celltype -----------------------------------------------

ggplot(mito_means_tidy, aes(
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
    legend.position=c(.17,.85),
    text = element_text(family = "sans", size = 12)
  )



# plot mitochondria positions --------------------------------------------------

plot_background()
plot3d(neurons, color = "lightgrey", lwd = 2, soma = TRUE, alpha = 0.2)

pos_ves_syn <- mito_positions |>
  filter(mito_type == "ves_syn") |>
  select(x, y, z)
plot3d(pos_ves_syn,
  add = TRUE,
  col = "#EE6677",
  size = 5,
  alpha = 1
)

pos_ves_no_syn <- mito_positions |>
  filter(mito_type == "ves_no_syn") |>
  select(x, y, z)
plot3d(pos_ves_no_syn,
  add = TRUE,
  col = "#DDAA33",
  size = 5,
  alpha = 1
)

pos_ves_unc <- mito_positions |>
  filter(mito_type == "ves_unc") |>
  select(x, y, z)
plot3d(pos_ves_unc,
  add = TRUE,
  col = "#4477AA",
  size = 5,
  alpha = 1
)

pos_ves_no <- mito_positions |>
  filter(mito_type == "ves_no") |>
  select(x, y, z)
plot3d(pos_ves_no,
  add = TRUE,
  col = "black",
  size = 5,
  alpha = 1
)
