# code to measure length of cilia
# WARNING: It ia absolutely essential that the root node is proximal to the basal body 
# i.e., the order of nodes needs to be "root" -> "basal body" -> "cilium tip"
# it should work fine if skeleton has multiple cilia

library(catmaid)

conn <- source("~/R/conn.R")

tags <- catmaid_get_label_stats(pid = 35)

# using "cilium tip" tag because that shows that cilium is complete, cut off cilia don't have that tag
# measuring from tips because it's easy to get parent node, but harder to get child
ctips <- tags[tags$labelName == "cilium tip", "treenodeID"]

cilium_length <- function(ctip_id) {
  # this will only work if root node is proximal to basal body tag
  # treenode compact detail response:
  # [ID, parent ID, x, y, z, confidence, radius, skeleton_id, edition_time, user_id]
  path <- paste("/35/treenodes/", ctip_id, "/compact-detail", sep = "")
  ctip_treenode <- catmaid_fetch(path = path)
  skid <- ctip_treenode[[8]]
  neuron <- read.neuron.catmaid(skid, pid = 35)
  # check which treenodes have basal body annotation, because this is the stop point
  basalb_ids <- neuron$tags$`basal body`
  treenode_id <- ctip_id
  distance <- 0
  while (!(treenode_id %in% basalb_ids)) {
    index <- match(treenode_id, neuron$d$PointNo)
    parent_id <- neuron$d$Parent[[index]]
    if (parent_id == -1) {
      print("root reached before basal body tag, from ciium tip")
      print(ctip_id)
      stop()
    }
    index_parent <- match(parent_id, neuron$d$PointNo)
 
    distance_2nodes <- dist(rbind(c(neuron$d$X[[index]],
                                neuron$d$Y[[index]],
                                neuron$d$Z[[index]]),
                              c(neuron$d$X[[index_parent]],
                                neuron$d$Y[[index_parent]],
                                neuron$d$Z[[index_parent]])
                                ))
    distance <- distance + distance_2nodes
    #print(distance)
    treenode_id <- parent_id
  }
  return(distance)
}

cilium_lenghts <- list()
for (ctip_id in ctips) {
  c_length <- cilium_length(ctip_id)
  cilium_lenghts <- c(cilium_lenghts, ctip_id = c_length)
}

hist(as.numeric(cilium_lenghts))
