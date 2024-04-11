source("analysis/scripts/packages_and_functions.R")

# find skeletons without annotations -------------------------------------------
skids <- unlist(
  catmaid_fetch(path = "/35/skeletons/"))

annot <- catmaid_get_annotations_for_skeletons(skids, pid = 35)

# no annotations at all
skids_with_annot <- annot$skid |> unique()
skids_without_annot <- setdiff(skids, skids_with_annot)

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
  path <- paste(pid, "/treenodes/compact-detail", sep = "")
  treenode_detail <- catmaid_fetch(path = path, 
                                   body = list(treenode_ids=treenodeid))
  xyz <- data.frame(x = treenode_detail[[1]][[3]],
                    y = treenode_detail[[1]][[4]],
                    z = treenode_detail[[1]][[5]])
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
