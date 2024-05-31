source("analysis/scripts/packages_and_functions.R")

# find skeletons without annotations -------------------------------------------
skids <- unlist(
  catmaid_fetch(path = "/35/skeletons/"))

annot <- catmaid_get_annotations_for_skeletons(skids, pid = 35)

# no annotations at all
skids_with_annot <- annot$skid |> unique()
# skids_without_annot <- setdiff(skids, skids_with_annot)

# without specific annotations
no_Q_annot <- annot |> filter(annotation != "Q1" & annotation != "Q2" & annotation != "Q3" & annotation != "Q4" ) |> 
  select(skid) |> pull() |> unique()


# centrioles -------------------------------------------------------------------
read_smooth_neuron <- function(annotation){
  nlapply(read.neurons.catmaid(annotation, pid = 35),
          function(x)
            smooth_neuron(x, sigma = 1000))
}

balancer <- read_smooth_neuron("celltype:balancer")
bridge <- read_smooth_neuron("celltype:bridge")

plot_background()
plot3d(
  bridge, soma = TRUE, color = Okabe_Ito[c(1,5,8)], 
  alpha = 0.3, lwd = c(4,3,3)
)
plot3d(
  balancer, soma = TRUE, color = Okabe_Ito[8], 
  alpha = 0.1, lwd = 3
)


get_treenode_pos <- function(treenodeid, pid) {
  path <- paste("/", pid, "/node/get_location", sep = "")
  pos <- catmaid_fetch(path = path,
                       body = list(tnid=treenodeid))
  xyz <- data.frame(x = pos[[2]],
                    y = pos[[3]],
                    z = pos[[4]])
}


plot_tag <- function(tagname, pid, color) {
  tags_all <- catmaid_get_label_stats(pid = pid)
  treenodes_with_tags_all <- tags_all |> 
    filter(labelName==tagname) |> 
    select(treenodeID) |>
    unlist() |>
    unique()
  
  pos <- lapply(treenodes_with_tags_all, get_treenode_pos, 35) |>
    bind_rows()
  
  plot3d(pos,
         add = TRUE, 
         col=color, 
         size=5,
         alpha=0.9)
}

plot_tag("centriole", 35, "magenta")
plot_tag("basal body", 35, "cyan")




# find cilium tips which are not tagged ----------------------------------------
tags <- catmaid_get_label_stats(pid = 35)

skids <- tags[tags$labelName == "basal body", "skeletonID"] |> unique()

cilium_annot_mistake <- function(skid) {
  neurons_info <- read.neurons.catmaid(skid, pid = 35)
  neuron_info <- neurons_info[[1]]
  annotations <- neurons_info %>% attr("anndf")
  if (length(annotations) == 0) {
    c_pocket <- 0
  } else {
    c_pocket <- annotations %>% select(annotation) %>% 
      filter(annotation == "ciliary_pocket_yes") %>% pull() %>% length()
  }
  neuron_name <- catmaid_fetch(path = "/35/skeleton/neuronnames",
                               body = list(skids=skid))[[1]]
  n_bb <- neuron_info$tags$`basal body` |> length()
  n_ctip <- neuron_info$tags$`cilium tip` |> length()
  n_ccut <- neuron_info$tags$`cilium cut` |> length()
  n_ct_extra <- neuron_info$tags$`cilium tip extracellular` |> length()
  n_c_exit <- neuron_info$tags$`exit_ciliary_pocket` |> length()
  if (n_ctip + n_ccut != n_bb || (c_pocket == 1 && n_ct_extra != n_c_exit)) {
    return(list(skid=skid, 
               neuron_name=neuron_name, 
               n_ctip=n_ctip, 
               n_ccut=n_ccut, 
               n_bb=n_bb,
               n_ct_extra=n_ct_extra,
               n_c_exit=n_c_exit))
  }
}

c_annot_mistakes <- lapply(skids, cilium_annot_mistake) %>% bind_rows()
