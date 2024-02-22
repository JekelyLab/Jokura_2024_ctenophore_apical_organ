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

# get connector ids
connector_ids <- c()
for(i in 1:length(contacts$connectors)){
  connector_ids[i] <- contacts$connectors[[i]][[1]]
  }
connector_ids

partner1_skid <- c()
for(i in 1:length(contacts$connectors)){
  partner1_skid[i] <- contacts$partners[[i]][[1]][[3]]
}
partner1_skid

partner2_skid <- c()
for(i in 1:length(contacts$connectors)){
  partner2_skid[i] <- contacts$partners[[i]][[2]][[3]]
}
partner2_skid

node_id_of_contact_partner1 <- c()
for(i in 1:length(contacts$connectors)){
  node_id_of_contact_partner1[i] <- contacts$partners[[i]][[1]][[2]]
}
node_id_of_contact_partner1

node_id_of_contact_partner2 <- c()
for(i in 1:length(contacts$connectors)){
  node_id_of_contact_partner2[i] <- contacts$partners[[i]][[2]][[2]]
}
node_id_of_contact_partner2


#check for tags enter_cil and enter_nocil
partner1_pre_cil <- c()
partner1_pre_nocil <- c()
for(i in 1:length(contacts$connectors)){
  neuron <- read.neuron.catmaid(partner1_skid[i], pid = 35)
  partner1_pre_cil[i] <- node_id_of_contact_partner1[i] %in% neuron$tags$enter_cil
  partner1_pre_nocil[i] <- node_id_of_contact_partner1[i] %in% neuron$tags$enter_nocil
}
partner1_pre_cil
partner1_pre_nocil

partner2_pre_cil <- c()
partner2_pre_nocil <- c()
for(i in 1:length(contacts$connectors)){
  neuron <- read.neuron.catmaid(partner2_skid[i], pid = 35)
  partner2_pre_cil[i] <- node_id_of_contact_partner2[i] %in% neuron$tags$enter_cil
  partner2_pre_nocil[i] <- node_id_of_contact_partner2[i] %in% neuron$tags$enter_nocil
}
partner2_pre_cil
partner2_pre_nocil

# define connectivity tibble

connector_tb <- tibble(
  connectors = connector_ids,
  skid1 = partner1_skid,
  skid2 = partner2_skid,
  partner1_pre_cil = partner1_pre_cil,
  partner1_pre_nocil = partner1_pre_nocil,
  partner2_pre_cil = partner2_pre_cil,
  partner2_pre_nocil = partner2_pre_nocil
)
connector_tb

connector_tb |>
  select(skid1, skid2) |>
  pull()

contactome_tbl_graph 

