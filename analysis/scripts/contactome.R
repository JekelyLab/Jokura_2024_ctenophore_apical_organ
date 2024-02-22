source("analysis/scripts/packages_and_functions.R")

# get all abutting contacts from CATMAID ------------

contacts <- catmaid_fetch(
  path = "35/connectors/",
  body = list(
    relation_type = "abutting",
    with_partners = "true"
  )
)
length(contacts$connectors)

skid1 <- contacts$partners[[1]][[1]][[3]]

skid2 <- contacts$partners[[1]][[2]][[3]]

