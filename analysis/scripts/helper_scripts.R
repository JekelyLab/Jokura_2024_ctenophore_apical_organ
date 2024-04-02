source("analysis/scripts/packages_and_functions.R")

# find skeletons without annotations ----------------------------------------
skids <- unlist(
  catmaid_fetch(path = "/35/skeletons/"))

annot <- catmaid_get_annotations_for_skeletons(skids, pid = 35)

# no annotations at all
skids_with_annot <- annot$skid |> unique()
skids_without_annot <- setdiff(skids, skids_with_annot)

# without specific annotations
no_Q_annot <- annot |> filter(annotation != "Q1" & annotation != "Q2" & annotation != "Q3" & annotation != "Q4" ) |> 
  select(skid) |> pull() |> unique()
