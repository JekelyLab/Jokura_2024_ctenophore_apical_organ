library(catmaid)
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization

conn <- source("~/R/conn.R")

skids <- unlist(
  catmaid_fetch(path = "/35/skeletons/"))

characters <- list("soma",
                   "mitochondrion",
                   "centriole",
                   "basal body",
                   "9+0",
                   "9+2",
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

write.csv(tag_stats_per_skid, "analysis/data/organelle_stats.csv")

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


# visualize different cells ----------------------------------------------------

# all columns are characters, which means that many math functions don't work
# convert comlumns with numers to int
tag_stats_per_skid <-  tag_stats_per_skid |> mutate_at(c(3:21), as.integer)

centriolerich_skids <- tag_stats_per_skid |> filter(centriole > 9) |> 
  select(skid) |> pull()

centriolerich <- nlapply(read.neurons.catmaid(centriolerich_skids, pid = 35),
          function(x) smooth_neuron(x, sigma = 1000))

balancer <- read_smooth_neuron("celltype:balancer")
bridge <- read_smooth_neuron("celltype:bridge")

plot_background()
plot3d(
  centriolerich, soma = TRUE, color = "blue", 
  alpha = 0.3, lwd = c(4,3,3)
)
plot3d(
  balancer, soma = TRUE, color = Okabe_Ito[8], 
  alpha = 0.1, lwd = 3
)

# get cells where all centrioles are annotated and color them by number of centrioles

centriole_skels <- read.neurons.catmaid("centrioles done", pid = 35)
Reds <- brewer.pal(9, 'Reds')
for (centriole_skel in centriole_skels) {
  skel <- smooth_neuron(centriole_skel, sigma=6000)
  ccount <- length(centriole_skel$tags$centriole)
  # pick discrete colors based on normalised degree
  shade <- round(ccount/max(tag_stats_per_skid$centriole)*10)
  # if there are no centrioles, skip it
  if (shade == 0) {
    next
  }
  plot3d(skel, soma = TRUE, lwd = 1, add=T, 
         alpha = 0.7, 
         col=Reds[shade])
}



# mitochondria stats -----------------------------------------------------------

mito_done_cells <- read.neurons.catmaid("mitochondria done", pid=35, fetch.annotations = TRUE)
annotations <- attr(mito_done_cells, 'anndf') |> as.tibble()

mito_counts_tb <- tibble(skid=integer(),
                         neuron_name=character(),
                         celltype=character(),
                         n_mito=integer())
for (mito_done_cell in mito_done_cells) {
  # var named sskid because othrwise filter doesn't work
  sskid <- mito_done_cell$skid
  n_mito <- length(mito_done_cell$tags$mitochondrion)
  celltype <- annotations |> 
    filter(skid==sskid) |> 
    filter(grepl("celltype:", annotation)) |> 
    select(annotation) |> 
    pull()
  celltype <- gsub(".*:","", celltype)
  if (length(celltype)==0) {
    celltype <- "NA"
  }
  neuron_name <- catmaid_get_neuronnames(sskid, pid = 35)
  print(paste(sskid, neuron_name, celltype, n_mito))
  mito_counts_tb <- mito_counts_tb |> add_row(skid=sskid, neuron_name=neuron_name, celltype=celltype, n_mito=n_mito)
}
write.csv(mito_counts_tb, "analysis/data/mito_per_celltype.csv")



# bridge cells --------------------------------------------------------
bridge_cells <- read.neurons.catmaid("celltype:bridge", pid=35, fetch.annotations = TRUE)
annotations <- attr(bridge_cells, 'anndf') |> as.tibble()

bridge_stats_tb <- tibble(skid=integer(),
                   neuron_name=character(),
                   n_centrioles=integer(),
                   centrioles_done=logical(),
                   n_mito=integer(),
                   mito_done=logical(),
                   syn_in=integer(),
                   syn_out=integer())
for (bridge_cell in bridge_cells) {
  # var named sskid because otherwise filter doesn't work
  sskid <- bridge_cell$skid
  neuron_name <- catmaid_get_neuronnames(sskid, pid = 35)
  n_centrioles <- length(bridge_cell$tags$centriole)
  n_mito <- length(bridge_cell$tags$mitochondrion)
  centrioles_done <- annotations |> 
    filter(skid==sskid) |> 
    filter(annotation=="centrioles done") |> 
    select(annotation) |> 
    pull() |>
    length() |>
    as.logical()
  mito_done <- annotations |> 
    filter(skid==sskid) |> 
    filter(annotation=="mitochondria done") |> 
    select(annotation) |> 
    pull() |>
    length() |>
    as.logical()
  # pre is 0 and post is 1
  prepost <- bridge_cell$connectors$prepost
  n_pre <- sum(prepost == 0)
  n_post <- sum(prepost == 1)
  bridge_stats_tb <- bridge_stats_tb |> add_row(skid=sskid, 
                                  neuron_name=neuron_name, 
                                  n_centrioles=n_centrioles,
                                  centrioles_done=centrioles_done,
                                  n_mito=n_mito,
                                  mito_done=mito_done,
                                  syn_in=n_pre,
                                  syn_out=n_post)
}
write.csv(bridge_stats_tb, "analysis/data/bridge_stats.csv")
