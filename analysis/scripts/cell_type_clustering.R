library(catmaid)
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization

conn <- source("~/R/conn.R")

skids <- unlist(
  catmaid_fetch(path = "/35/skeletons/"))

characters <- list("soma",
                   "centriole",
                   "basal body",
                   "rootlets",
                   "golgi",
                   "lysosome?",
                   "multivesicular body",
                   "SER",
                   "RER",
                   "plume",
                   "large clear vesicles",
                   "small clear vesicles",
                   "vesicles",
                   "large electrondense vesicles",
                   "small dark dots",
                   "cilium tip intracellular",
                   "cilium tip inside different cell")

tag_stats_per_skid <- data.frame(matrix(ncol = length(characters)+2, nrow = 0))
colnames(tag_stats_per_skid) <- c("skid", "neuron_name", characters)

tag_count <- function(tagl, character) {
  return(length(tagl[[character]]))
}

tag_stats <- function(skid) {
  # if skeleton has only one node read.neuron.catmaid gives an error
  try(neuron_info <- read.neuron.catmaid(skid, pid = 35), silent=TRUE)
  if("try-error" %in% class(t)) return(NULL)
  tags <- neuron_info$tags
  
  counts <- mapply(tag_count,
                   list(tags), 
                   characters)

  # ignore cells without a soma
  if (counts[[1]] == 0) {
    return(NULL)
  }
  else {
    return(counts)
  }
}

for (skid in skids) {
  neuron_name <- catmaid_get_neuronnames(skid, pid = 35)
  tag_counts <- tag_stats(skid)
  if (!is.null(tag_counts)) {
    tag_stats_per_skid[nrow(tag_stats_per_skid) + 1,] = c(skid, neuron_name, tag_counts)
  }
}

tag_stats_noskid <- tag_stats_per_skid[,3:ncol(tag_stats_per_skid)]

km <- kmeans(tag_stats_noskid, centers = 20, nstart = 25)
fviz_cluster(km, data = tag_stats_per_skid)


fviz_cluster(km,# clustering result 
             data = tag_stats_noskid, # data 
             ellipse.type = "convex", 
             star.plot = TRUE, 
             repel = TRUE, 
             ggtheme = theme_minimal()
) 


tag_stats_per_skid %>%
  as_tibble() %>%
  mutate(cluster = km$cluster,
         state = row.names(skid)) %>%
  ggplot(aes(SER, centriole, color = factor(cluster), label = skid)) +
  geom_text()

celltypes <- data.frame()
for (i in seq(1:length(tag_stats_per_skid$skid))) {
  skid <- tag_stats_per_skid$skid[[i]]
  nname <- unlist(
    catmaid_fetch(path = "/35/skeleton/neuronnames",
                  body = list(skids = skid)))
  new_row <- c(
    skid,
    unname(nname),
    km$cluster[[i]]
  )
  celltypes <- rbind(celltypes, new_row)
}

sorted <- celltypes[order(celltypes[[3]]),]

