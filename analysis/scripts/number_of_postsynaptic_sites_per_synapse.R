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

synaptic_neurons <- read.neurons.catmaid("synaptic neuron", pid=35)

mito_in_syn <- data.frame()
for (i in seq_along(synaptic_neurons)) {
  n_mito_vesicles_syn <- 0
  n_mito_vesicles_no_syn <- 0
  mito_vesicles <- synaptic_neurons[[i]]$tags$`mitochondrion vesicles`
  # this doesn't check which type of connector is associated with the node,
  # but in 99.5% cases it will be presynaptic, so it's probably not necessary to do extra checks
  syn <- synaptic_neurons[[i]]$connectors$treenode_id
  mito_syn <- list()
  for (mito_vesicle in mito_vesicles) {
    if (mito_vesicle %in% syn) {
      n_mito_vesicles_syn <- n_mito_vesicles_syn + 1
    }
    else {
      n_mito_vesicles_no_syn <- n_mito_vesicles_no_syn + 1
    }
  }
  mito_unclear_vesicles <- synaptic_neurons[[i]]$tags$`mitochondrion unclear vesicles`
  n_mito_unclear_vesicles <- length(mito_unclear_vesicles)
  mito_no_vesicles <- synaptic_neurons[[i]]$tags$`mitochondrion no vesicles`
  n_mito_no_vesicles <- length(mito_no_vesicles)
  n_mito_total <- n_mito_vesicles_syn + n_mito_vesicles_no_syn + n_mito_unclear_vesicles + n_mito_no_vesicles
  mito_count_neuron <- data.frame(skid=synaptic_neurons[[i]]$NeuronName,
                                  vesicles_syn=n_mito_vesicles_syn,
                                  vesicles_no_syn=n_mito_vesicles_no_syn,
                                  vesicles_unclear=n_mito_unclear_vesicles,
                                  vesicles_no=n_mito_no_vesicles,
                                  total=n_mito_total)
  mito_in_syn <- rbind(mito_in_syn, mito_count_neuron)
}

mito_in_syn_arranged <- arrange(mito_in_syn, desc(total))

mito_in_syn_tidy <-  mito_in_syn_arranged %>%
  #select(total) %>%
  pivot_longer(
    cols = c("vesicles_syn", "vesicles_no_syn", "vesicles_unclear", "vesicles_no"), 
    names_to = "characteristic", 
    values_to = "count")

mito_in_syn_graph <- mito_in_syn_tidy %>%
  ggplot(aes(as.character(skid),
             count,
             fill=factor(characteristic,
                         levels=c("vesicles_syn", "vesicles_no_syn", "vesicles_unclear", "vesicles_no")))) +
  geom_bar(position="fill",
           stat = "identity") +
  scale_x_discrete(limits = as.character(mito_in_syn$skid)) +
  scale_y_reverse() +
  geom_text(aes(label = count),
            position = "fill",
            hjust = 1,
            size = 2,
            family="sans") +
  coord_flip() + 
  scale_fill_manual(values = c("#4477AA","#AA3377","#EE6677","lightgrey")) + # Tol
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

#ggsave("pictures/mito_vesicles_graph.png", limitsize = FALSE, 
#       units = c("px"), desmo_tono_graph, width = 2400, height = 1600, bg = 'white')

