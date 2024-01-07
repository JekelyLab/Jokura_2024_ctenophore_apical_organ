
library(catmaid)

conn <- source("~/R/conn.R")


pre_connectors <- catmaid_fetch(
  path = "35/connectors/links",
  body = list(
    relation_type = "presynaptic_to",
    with_partners = "false"
  )
)

n_pre_post <- function(connector_id) {
  connector_info <- catmaid_fetch(path = paste("35/connectors/", connector_id, sep = ""))
  # pre should always be 1, but maybe it's good to check if there is any weird stuff going on
  pre <- connector_info$partners %like% 'presynaptic_to'
  n_pre <- table(pre)["TRUE"]
  post <- connector_info$partners %like% 'postsynaptic_to'
  n_post <- table(post)["TRUE"]
  return(n_post)
}
