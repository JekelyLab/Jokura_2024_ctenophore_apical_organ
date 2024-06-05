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
  pre <- connector_info$partners %like% 'presynaptic_to'
  n_pre <- table(pre)["TRUE"][[1]]
  post <- connector_info$partners %like% 'postsynaptic_to'
  n_post <- table(post)["TRUE"][[1]]
  return(n_post)
}

pre_connector_IDs <- lapply(pre_connectors$connectors,'[[',1)

n_post <- lapply(pre_connector_IDs, n_pre_post)

summary(as.numeric(n_post))

# calculate percentage of mitochondria which have vesicles or synapses associated with them -----

get_mito_stats <- function(neur) {
  # accepts neuron as input
  n_mito_vesicles_syn <- 0
  n_mito_vesicles_no_syn <- 0
  mito_vesicles <- neur$tags$`mitochondrion vesicles`
  if (length(neur$connectors) > 0) {
    syn <- neur$connectors |> filter(prepost == 0) |> select(treenode_id) |> pull()
    mito_syn <- list()
    for (mito_vesicle in mito_vesicles) {
      if (mito_vesicle %in% syn) {
        n_mito_vesicles_syn <- n_mito_vesicles_syn + 1
      }
      else {
        n_mito_vesicles_no_syn <- n_mito_vesicles_no_syn + 1
      }
    } 
  } else {
    n_mito_vesicles_no_syn <- length(mito_vesicles)
  }
  mito_unclear_vesicles <- neur$tags$`mitochondrion unclear vesicles`
  n_mito_unclear_vesicles <- length(mito_unclear_vesicles)
  mito_no_vesicles <- neur$tags$`mitochondrion no vesicles`
  n_mito_no_vesicles <- length(mito_no_vesicles) 
  
  # for cells without synapses and mitochondria with vesicles, all mitochondria
  # only have generic "mitochondrion" tag
  if (n_mito_no_vesicles == 0) {
    mito_no_vesicles <- neur$tags$`mitochondrion`
    n_mito_no_vesicles <- length(mito_no_vesicles) 
  }
  
  sskid=neur$NeuronName
  celltype <- catmaid_get_annotations_for_skeletons(sskid, pid=35) |> 
    select(annotation) |> 
    filter(grepl("celltype:", annotation)) |> 
    pull()
  celltype <- gsub(".*:","", celltype)
  if (length(celltype)==0) {
    celltype <- "NA"
  }
  
  n_mito_total <- n_mito_vesicles_syn + n_mito_vesicles_no_syn + n_mito_unclear_vesicles + n_mito_no_vesicles
  list(celltype=celltype,
       skid=sskid,
       vesicles_syn=n_mito_vesicles_syn,
       vesicles_no_syn=n_mito_vesicles_no_syn,
       vesicles_unclear=n_mito_unclear_vesicles,
       vesicles_none=n_mito_no_vesicles,
       total=n_mito_total)
}


mito_done <- read.neurons.catmaid("mitochondria done", pid = 35)
mito_stats <- lapply(mito_done, get_mito_stats) |> bind_rows() |>
  arrange(desc(total))


#---------------------------------


# TODO - plot averages instead, and add number of cells per category


ves_syn_tbl <- mito_stats %>% group_by(celltype) %>%
  summarize(mean_vesicles_syn = mean(vesicles_syn)) %>% 
  arrange(celltype) %>% select(mean_vesicles_syn)
ves_no_syn_tbl <- mito_stats %>% group_by(celltype) %>%
  summarize(mean_vesicles_no_syn = mean(vesicles_no_syn)) %>%
  arrange(celltype) %>% select(mean_vesicles_no_syn)
ves_unclear_tbl <- mito_stats %>% group_by(celltype) %>%
  summarize(mean_vesicles_unclear = mean(vesicles_unclear)) %>%
  arrange(celltype) %>% select(mean_vesicles_unclear)
ves_none_tbl <- mito_stats %>% group_by(celltype) %>% 
  summarize(mean_vesicles_none = mean(vesicles_none)) %>%
  arrange(celltype) %>% select(mean_vesicles_none)
celltype_numbers <- mito_stats %>% add_count(celltype) %>%
  select(celltype, n) %>% unique() %>% arrange(celltype)
mito_means <- bind_cols(celltype_numbers,
                        ves_syn_tbl, ves_no_syn_tbl,
                        ves_unclear_tbl, ves_none_tbl)


mito_means_tidy <-  mito_means %>%
  #select(total) %>%
  pivot_longer(
    cols = c("mean_vesicles_syn", "mean_vesicles_no_syn", "mean_vesicles_unclear", "mean_vesicles_none"), 
    names_to = "characteristic", 
    values_to = "value")


mito_in_syn_graph <- mito_in_syn_tidy %>%
  ggplot(aes(as.character(celltype),
             count,
             fill=factor(characteristic,
                         levels=c("vesicles_syn", "vesicles_no_syn", "vesicles_unclear", "vesicles_none")))) +
  geom_bar(position="fill",
           stat = "identity") +
  #scale_x_discrete(limits = as.character(mito_stats$celltype)) +
  #scale_y_reverse() +
  #geom_text(aes(label = count), size = 3, hjust = 0.5, vjust = 3, position = "fill") +
  # geom_text(aes(label = count),
  #           position = "fill",
  #           hjust = 1,
  #           size = 2,
  #           family="sans") +
  #coord_flip() + 
  scale_fill_manual(values = c("#EE6677", "#AA3377", "#4477AA", "lightgrey")) + # Tol
  theme_bw() +
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_blank(),
        text=element_text(family="sans", size = 12))
mito_in_syn_graph



ggplot(mito_means_tidy, aes(fill=factor(characteristic,
                                        levels=c("mean_vesicles_syn", "mean_vesicles_no_syn", "mean_vesicles_unclear", "mean_vesicles_none")),
                            y=value, x=interaction(celltype,n))) + 
  geom_bar(position="stack", stat="identity") +
  theme_bw() +
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_blank(),
        text=element_text(family="sans", size = 12))


# get positions of mitochondria ------------------------------------------------

get_pos_of_tags_in_neuron <- function(neur, tagname) {
  # argument is neuron (data type), not skid
  tag_positions <- tibble(treenodeid=character(),
                          x=double(),
                          y=double(),
                          z=double())
  tags <- neur$tags[[tagname]]
  for (tag in tags) {
    pos <- neur$d |> filter(PointNo == tag) |> select(X,Y,Z)
    tag <- as.character(tag)
    tag_positions <- tag_positions |> bind_rows(list(treenodeid=tag,
                                                     x=pos$X, y=pos$Y, z=pos$Z))
  }
  return(tag_positions)
}

mito_positions <- tibble(celltype=character(),
                         skid=character(),
                         mito_type=character(),
                         treenodeid=character(),
                         x=double(),
                         y=double(),
                         z=double()
)

for (neur in neurons) {
  sskid=neur$NeuronName
  celltype <- catmaid_get_annotations_for_skeletons(sskid, pid=35) |> 
    select(annotation) |> 
    filter(grepl("celltype:", annotation)) |> 
    pull()
  celltype <- gsub(".*:","", celltype)
  if (length(celltype)==0) {
    celltype <- "NA"
  }
  
  ves_no <- get_pos_of_tags_in_neuron(neur, "mitochondrion no vesicles")
  ves_no <- bind_cols(mito_type="ves_no", ves_no)
  
  ves_unc <- get_pos_of_tags_in_neuron(neur, "mitochondrion unclear vesicles")
  ves_unc <- bind_cols(mito_type="ves_unc", ves_unc)
  
  ves_yes <- get_pos_of_tags_in_neuron(neur, "mitochondrion vesicles")
  ves_yes <- bind_cols(mito_type="ves", ves_yes)
  # treenodes with outgoing synapses
  syn_out_treenodes <- neur$connectors |> 
    filter(prepost == 0) |> select(treenode_id) |> pull()
  # check if the mitochondrion is associated with a synapse
  for (treenode in ves_yes$treenodeid) {
    if (treenode %in% syn_out_treenodes) {
      mitotype <- "ves_syn"
    }
    else {
      mitotype <- "ves_no_syn"
    }
    ves_yes <- ves_yes |>
      mutate(mito_type = replace(mito_type, treenodeid == treenode, mitotype))
  }
  mito_pos <- bind_rows(ves_no, ves_unc, ves_yes)
  mito_pos <- bind_cols(celltype=celltype,
                        skid=sskid,
                        mito_pos)
  mito_positions <- bind_rows(mito_positions, mito_pos)
}


# plot mitochondria positions --------------------------------------------------

plot_background()
plot3d(neurons, color = "lightgrey", lwd = 2, soma = TRUE, alpha = 0.2)

pos_ves_syn <- mito_positions |> filter(mito_type=="ves_syn") |> select(x,y,z)
plot3d(pos_ves_syn,
       add = TRUE, 
       col="#EE6677",
       size=5,
       alpha=1)

pos_ves_no_syn <- mito_positions |> filter(mito_type=="ves_no_syn") |> select(x,y,z)
plot3d(pos_ves_no_syn,
       add = TRUE, 
       col="#AA3377", 
       size=5,
       alpha=1)

pos_ves_unc <- mito_positions |> filter(mito_type=="ves_unc") |> select(x,y,z)
plot3d(pos_ves_unc,
       add = TRUE, 
       col="#4477AA", 
       size=5,
       alpha=1)




pos_ves_no <- mito_positions |> filter(mito_type=="ves_no") |> select(x,y,z)
plot3d(pos_ves_no,
       add = TRUE, 
       col="black", 
       size=5,
       alpha=1)

