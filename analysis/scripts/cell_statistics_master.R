source("analysis/scripts/packages_and_functions.R")

# statistics (one number per cell): ---------------------------------------------
# - number of cilia
# - avg cilium length
# - avg ciliary pocket depth
# - number of centrioles
# - number of mitochondria - NA if not done
# - number of pre synapses
# - number of post synapses
# - synapse polyadicity
# - pre partner celltypes
# - post partner celltypes
# - celltype

stats_cilia <- read.csv("analysis/data/cilium_lengths.csv")
stats_organelles <- read.csv("analysis/data/organelle_stats.csv")

# process only cells with somas:
tags <- catmaid_get_label_stats(pid = 35)
skids <- tags %>% filter(labelName == "soma") %>% 
  select(skeletonID) %>% 
  unique %>% 
  pull


# get synapse info, but only process pre connectors so we don't duplicate things
get_syn_stats <- function(skid) {
  nneuron <- read.neuron.catmaid(skid, pid = 35)
  connectors_info <- nneuron$connectors
  if (length(connectors_info) == 0) {
    connectors_info <- NULL
  } else {
  connectors_info$skid <- nneuron$skid
  }
  return(connectors_info)
}

stats_synapse <- data.frame()
for (skid in skids) {
  stats_synapse <- rbind(stats_synapse, get_syn_stats(skid))
}

stats_synapse <- stats_synapse %>% group_by(connector_id, skid) %>% slice(1)

write.csv(stats_synapse, "analysis/data/stats_synapse.csv")
stats_synapse <- read.csv("analysis/data/stats_synapse.csv")

# compile stats ------------------------------
get_celltype_annot_for_skid <- function(skid) {
  annot <- catmaid_fetch(path = "/35/annotations/forskeletons",
                         body = list(skeleton_ids=skid))
  result <- annot$annotations %>%
    map(~ str_extract(.x, "(?<=celltype:).*")) %>%
    unlist() %>%
    `[`(., . != "") %>%
    na.omit()
  if (length(result) > 0) {
    celltype <- as.character(result)
  } else {
    celltype <- NA
  }
  return(celltype)
}

get_stats <- function(skid) {
  sskid=skid
  annot <- catmaid_fetch(path = "/35/annotations/forskeletons",
                         body = list(skeleton_ids=sskid))
  # get celltype from server because none of the current tables have annotations for all cells
  celltype <- get_celltype_annot_for_skid(sskid)
  basal_bodies <- stats_organelles %>% filter(skid==sskid) %>% select(basal.body) %>% pull() %>% sum()
  avg_cil_length <- stats_cilia %>% 
    filter(skid==sskid) %>% 
    select(cil_length) %>% 
    pull() %>% 
    mean(na.rm = TRUE) %>% 
    ifelse(is.nan(.), NA, .)
  avg_cil_pocket_depth <- stats_cilia %>%
    filter(skid==sskid) %>% 
    select(pocket_length) %>% 
    pull() %>% 
    mean(na.rm = TRUE) %>% 
    ifelse(is.nan(.), NA, .)
  # get centrioles only if all annotated
  if (any(annot$annotations == "centrioles done")) {
    centrioles <- stats_organelles %>% filter(skid==sskid) %>% select(centriole) %>% pull() %>% sum()
  } else {
    centrioles <- NA
  }
  # only count mitochondria if mito done
  if (any(annot$annotations == "mitochondria done")) {
    mito <- stats_organelles %>% filter(skid==sskid) %>% select(mitochondrion) %>% pull() %>% sum()
  } else {
    mito <- NA
  }
  pre_syn_count <- stats_synapse %>% filter(skid == sskid) %>% filter(prepost == 0) %>% count() %>% as.numeric()
  post_syn_count <- stats_synapse %>% filter(skid == sskid) %>% filter(prepost == 1) %>% count() %>% as.numeric()
  # careful, pre-synaptic connectors from our skid are connecting to post partners
  # presynaptic partners
  post_connectors <- stats_synapse %>% 
    filter(skid == sskid) %>% 
    filter(prepost == 1) %>% 
    select(connector_id) %>% 
    pull() %>% 
    unique()
  pre_partner_skids <- stats_synapse %>% 
    filter(connector_id %in% post_connectors) %>%
    filter(prepost == 0) %>%
    select(skid) %>%
    pull() %>%
    unique()
  pre_partner_celltypes <- list()
  for (pre_partner_skid in pre_partner_skids) {
    pre_partner_celltypes <- c(pre_partner_celltypes, get_celltype_annot_for_skid(pre_partner_skid))
  }
  pre_partner_celltypes <- pre_partner_celltypes %>% unlist() %>% sort(decreasing = T) %>% unique() %>% list()
  # postsynaptic partners
  pre_connectors <- stats_synapse %>% 
    filter(skid == sskid) %>% 
    filter(prepost == 0) %>% 
    select(connector_id) %>% 
    pull() %>% 
    unique()
  post_partner_skids <- stats_synapse %>% 
    filter(connector_id %in% pre_connectors) %>%
    filter(prepost == 1) %>%
    select(skid) %>%
    pull() %>%
    unique()
  post_partner_celltypes <- list()
  for (post_partner_skid in post_partner_skids) {
    post_partner_celltypes <- c(post_partner_celltypes, get_celltype_annot_for_skid(post_partner_skid))
  }
  post_partner_celltypes <- post_partner_celltypes %>% unlist() %>% sort(decreasing = T) %>% unique() %>% list()
  # get polyadicity
  post_targets_count <- stats_synapse %>% 
    filter(connector_id %in% pre_connectors) %>%
    filter(prepost == 1) %>%
    nrow()
  polyadicity <- post_targets_count/length(pre_connectors)
  # use incoming/outgoing and upstream/downstream instead of pre/post to be more clear
  # because pre_syn and pre_partners are opposite and you have to think hard to make sense of it
  stats <- tibble(skid=sskid,
                  celltype=celltype,
                  cilia=basal_bodies,
                  avg_cil_length=avg_cil_length,
                  avg_cil_pocket_depth=avg_cil_pocket_depth,
                  centrioles=centrioles,
                  mitochondria=mito,
                  incoming_syn=post_syn_count,
                  outgoing_syn=pre_syn_count,
                  polyadicity=polyadicity,
                  upstream_celltypes=c(pre_partner_celltypes),
                  downstream_celltypes=c(post_partner_celltypes),
  )

  stats <- stats %>% mutate(
    upstream_celltypes = map(upstream_celltypes, ~ paste(., collapse = ", ")),
    downstream_celltypes = map(downstream_celltypes, ~ paste(., collapse = ", "))
  )
  return(stats)
}

stats_all <- bind_rows(lapply(skids, get_stats))

stats_all <- stats_all %>% mutate_if(is.numeric, ~ ifelse(is.nan(.), NA, .))
stats_all <- stats_all %>% mutate_if(is.character, ~ ifelse(is.null(.), NA, .))
# when calculating polyadicity it could be possible to divide by 0
# this issue could only arise as a tracing mistake, but anyway, let's convert any such cases to NA
stats <- stats_all %>% mutate(polyadicity = ifelse(abs(polyadicity) == Inf, NA, polyadicity))
# keep in mind that cilium length can be NA even if there is a cilium, 
# because some cilia extend outside of the scanned volume, so we cannot determine their length
stats_all$upstream_celltypes <- sapply(stats_all$upstream_celltypes, paste, collapse = ", ")
stats_all$downstream_celltypes <- sapply(stats_all$downstream_celltypes, paste, collapse = ", ")
write.csv(stats_all, "analysis/data/stats_master.csv")

