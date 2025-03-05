source("analysis/scripts/packages_and_functions.R")

skid_Q1234 <- 2496955
skid_Q12 <- 2436172
skid_Q34 <- 2436531

neuron_Q1234 <- read.neuron.catmaid(skid_Q1234, pid = 35)
neuron_Q12 <- read.neuron.catmaid(skid_Q12, pid = 35)
neuron_Q34 <- read.neuron.catmaid(skid_Q34, pid = 35)

cablelength_Q1234 <- summary(neuron_Q1234) %>%
  select(cable.length) %>%
  pull()

cablelength_Q12 <- summary(neuron_Q12) %>%
  select(cable.length) %>%
  pull()

cablelength_Q34 <- summary(neuron_Q34) %>%
  select(cable.length) %>%
  pull()

# compare number of mitochondria per cable length-----------------------------
mito_counts <- read.csv("analysis/data/mito_per_celltype.csv")

mito_Q1234 <- mito_counts %>%
  filter(neuron_name == "SSN_Q1-4") %>%
  select(n_mito) %>%
  pull()

mito_Q12 <- mito_counts %>%
  filter(neuron_name == "SSN_Q1-2") %>%
  select(n_mito) %>%
  pull()

mito_Q34 <- mito_counts %>%
  filter(neuron_name == "SSN_Q3-4") %>%
  select(n_mito) %>%
  pull()

mito_ratio_Q1234 <- mito_Q1234/cablelength_Q1234
mito_ratio_Q12 <- mito_Q12/cablelength_Q12
mito_ratio_Q34 <- mito_Q34/cablelength_Q34

paste(mito_ratio_Q1234, mito_ratio_Q12, mito_ratio_Q34) %>% print()

# compare number of outgoing synapses per cable length--------------------------

pre_Q1234 <- sum(neuron_Q1234$connectors$prepost == 0)
pre_Q12 <- sum(neuron_Q12$connectors$prepost == 0)
pre_Q34 <- sum(neuron_Q34$connectors$prepost == 0)

ratio_pre_Q1234 <- pre_Q1234/cablelength_Q1234
ratio_pre_Q12 <- pre_Q12/cablelength_Q12
ratio_pre_Q34 <- pre_Q12/cablelength_Q34

paste(ratio_pre_Q1234, ratio_pre_Q12, ratio_pre_Q34) %>% print()

# compare number of incoming synapses per cable length

post_Q1234 <- sum(neuron_Q1234$connectors$prepost == 1)
post_Q12 <- sum(neuron_Q12$connectors$prepost == 1)
post_Q34 <- sum(neuron_Q34$connectors$prepost == 1)

ratio_post_Q1234 <- post_Q1234/cablelength_Q1234
ratio_post_Q12 <- post_Q12/cablelength_Q12
ratio_post_Q34 <- post_Q12/cablelength_Q34

paste(ratio_post_Q1234, ratio_post_Q12, ratio_post_Q34) %>% print()

# compare appearance of synapses------------------------------------------------
dir.create("analysis/data/synapse_crop")

for (neu in neuronlist(neuron_Q1234, neuron_Q12, neuron_Q34)) {
  half_bb_size_xy=800
  half_bb_size_z=0
  print(neu)
  connectors_info <- neu$connectors
  for (cid in connectors_info$connector_id) {
    # xzy positions in the connectors table are of the connector,
    # not of the treenode
    x <- connectors_info %>% filter(connector_id == cid) %>% select(x) %>% pull()
    y <- connectors_info %>% filter(connector_id == cid) %>% select(y) %>% pull()
    z <- connectors_info %>% filter(connector_id == cid) %>% select(z) %>% pull()
    catmaid_fetch(
      path = "35/crop/",
      body = list(
        stack_ids=28,
        min_x=x-half_bb_size_xy,
        min_y=y-half_bb_size_xy,
        min_z=z-half_bb_size_z,
        max_x=x+half_bb_size_xy,
        max_y=y+half_bb_size_xy,
        max_z=z+half_bb_size_z,
        zoom_level=0,
        single_channel='true',
        rotationcw=0
      )
    )
  }
  
  print(paste("Found", length(connectors_info$connector_id), "connectors"))
  print("Downloading the following ids:")
  print(connectors_info$connector_id)
  print("This will take some time")
  
  # server needs time to process this
  Sys.sleep(70)
  
  # download stacks ---------------------------------------------------------
  # the exact file names are needed for download,
  # to get them we need to check messages from the server
  # originally I was checking which messages are new, and downloading those
  # but this could lead to false positives
  # so I decided that's it's probably safe to download the last n messages
  # where n is the number of jobs
  
  # It would be cool if I can put the synapse ID as file name,
  # but it's not guaranteed that jobs will return in the same order
  
  cat_messages <- catmaid_fetch(
    path = "messages/list"
  )
  
  for (i in seq_along(connectors_info$connector_id)) {
    cat_message <- cat_messages[[i]]$text
    link <- regmatches(cat_message,
                       regexpr("/crop/download/crop_.{6}.tiff", cat_message))
    filename <- regmatches(link,
                           regexpr("crop_.{6}.tiff", link))
    full_link <- paste("https:/catmaid.ex.ac.uk", link, sep = "")
    file_path <- paste("analysis/data/synapse_crop/", filename, sep = "")
    download.file(full_link, destfile = file_path)
    file_path_new <- paste("analysis/data/synapse_crop/", neu$skid, "-", filename, sep = "")
    file.rename(file_path, file_path_new)
  }
}

syn_img_Q1234 <- list.files(path = "analysis/data/synapse_crop", pattern = as.character(skid_Q1234))
pattern <- paste(skid_Q12, "|", skid_Q34, sep = "")
# sample() to randomize synapses from the 2 small neurons
syn_img_small <- list.files(path = "analysis/data/synapse_crop", pattern = pattern) %>% sample()


dir.create("analysis/data/syn_comparison_images")
library(foreach)
library(magick)
i=1
foreach(filename_Q1234 = syn_img_Q1234, filename_small = syn_img_small) %do% {
  path_Q1234 <- paste("analysis/data/synapse_crop/", filename_Q1234, sep = "")
  path_small <- paste("analysis/data/synapse_crop/", filename_small, sep = "")
  print(paste(path_Q1234, path_small))
  panel1 <- ggdraw() + draw_image(image_read(path_Q1234))
  panel2 <- ggdraw() + draw_image(image_read(path_small))
  layout <- "A#B"
  fig <- panel1 + panel2 + plot_layout(design = layout, heights = 800, widths = c(1, 0.05, 1))
  save_plot(paste("analysis/data/syn_comparison_images/", i, ".jpg", sep = ""), fig)
  i=i+1
}
