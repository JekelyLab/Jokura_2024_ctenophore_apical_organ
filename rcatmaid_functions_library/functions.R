library(catmaid)
library(data.table)
library(tidyverse)

# get list of cell types in project --------------------------------------------
# assumes they have an annotation starting with "celltype:"
get_celltypes <- function(pid) {
  annotations <- catmaid_get_annotationlist(pid = pid)
  celltypes <- annotations$annotations |> filter(grepl("^celltype:", name)) |>
    select(name) |> pull() |> sub(pattern="^celltype:", replacement="")
  return(celltypes)
}

# get skids with multiple annotations ------------------------------------------
# note that second argument needs a vector if you are using more than 1 annotation!
get_skids_with_annot <- function(pid, annotations) {
  skids_list <- lapply(annotations, catmaid_skids, pid=pid)
  intersection <- Reduce(intersect,skids_list)
  return(intersection)
}

# find skeletons without annotation --------------------------------------------
# if used without annotation name, returns skids without any annotation
get_skels_without_annot <- function(pid, annotationname = "") {
  skids <- unlist(catmaid_fetch(path = paste("/", pid, "/skeletons/", sep="")))
  # I'm using the above method of getting all skids instead of catmaid_skids(".*", pid)
  # because the above method is 100x faster (~0.06s vs ~6.05s for 988 skeletons)
  # it is not posible to omit the above line because 
  # catmaid_get_annotations_for_skeletons(".*", pid) doesn't return skeletons with no annotations at all
  skids_with_annot <- catmaid_get_annotations_for_skeletons(skids, pid = pid) |>
    filter(grepl(annotationname, annotation)) |> select(skid) |> 
    pull() |> unique()
  skids_without_annot <- setdiff(skids, skids_with_annot)
  return(skids_without_annot)
}


# crop substack ----------------------------------------------------------------
crop_substack <- function(tagname,
                         half_bb_size_x, half_bb_size_y, half_bb_size_z,
                         zoomlevel,
                         dest_dir,
                         pid, stackid) {
  # tag search only returns treenodes with tag, not connectors
  # however, we want to get either treenodes or connectors, so we have to use generic search
  tagname_edited <- gsub(" ", "%20", tagname)
  search_results <- catmaid_fetch(path = paste(pid, "/search?substring=", tagname_edited, sep = ""))
  # search doesn't accept regex, so we have to filter results
  ids=c() # keep track of the number of substacks we have to download from the server later
  # use ids rather than simple counter because it helps with debugging
  for (i in seq_along(search_results)) {
    if (search_results[[i]]$name == tagname && search_results[[i]]$class_name == "label") {
      for (j in seq_along(search_results[[i]]$nodes)) {
        ids <- c(ids, search_results[[i]]$nodes[[j]]$id)
        x <- search_results[[i]]$nodes[[j]]$x
        y <- search_results[[i]]$nodes[[j]]$y
        z <- search_results[[i]]$nodes[[j]]$z
        catmaid_fetch(
          path = paste(pid, "/crop/", sep = ""),
          body = list(
            stack_ids=stackid,
            min_x=x-half_bb_size_x,
            min_y=y-half_bb_size_y,
            min_z=z-half_bb_size_z,
            max_x=x+half_bb_size_x,
            max_y=y+half_bb_size_y,
            max_z=z+half_bb_size_z,
            zoom_level=zoomlevel,
            single_channel='true',
            rotationcw=0
          )
        )
      }
      for (j in seq_along(search_results[[i]]$connectors)) {
        ids <- c(ids, search_results[[i]]$connectors[[j]]$id)
        x <- search_results[[i]]$connectors[[j]]$x
        y <- search_results[[i]]$connectors[[j]]$y
        z <- search_results[[i]]$connectors[[j]]$z
        catmaid_fetch(
          path = paste(pid, "/crop/", sep = ""),
          body = list(
            stack_ids=stackid,
            min_x=x-half_bb_size_x,
            min_y=y-half_bb_size_y,
            min_z=z-half_bb_size_z,
            max_x=x+half_bb_size_x,
            max_y=y+half_bb_size_y,
            max_z=z+half_bb_size_z,
            zoom_level=zoomlevel,
            single_channel='true',
            rotationcw=0
          )
        )
      }
    }
  }
  
  if (length(ids) == 0) {
    print(paste("No treenodes or connectors tagged with", tagname, "were found"))
    return()
  }
  
  print(paste("Found", length(ids), "treenodes or connectors tagged with", tagname))
  print("Downloading the following ids:")
  print(ids)
  print("This will take some time")
  
  # server needs time to process this
  Sys.sleep(90)
  
  # download stacks ---------------------------------------------------------
  # the exact file names are needed for download,
  # to get them we need to check messages from the server
  
  # It would be cool if I can put the synapse ID as file name,
  # but it's not guaranteed that jobs will return in the same order
  
  cat_messages <- catmaid_fetch(
    path = "messages/list"
  )
  
  for (k in seq_along(ids)) {
    cat_message <- cat_messages[[k]]$text
    link <- regmatches(cat_message, 
                       regexpr("/crop/download/crop_.{6}.tiff", cat_message))
    filename <- regmatches(link, 
                           regexpr("crop_.{6}.tiff", link))
    full_link <- paste("https:/catmaid.ex.ac.uk", link, sep = "")
    destpath <- paste(dest_dir, "/", filename, sep = "")
    download.file(full_link, destfile = destpath)
  }
}


# number of postsynaptic sites of a connector ----------------------------------
n_pre_post <- function(connector_id, pid) {
  connector_info <- catmaid_fetch(path = paste(pid, "/connectors/", connector_id, sep = ""))
  # pre should always be 1, but maybe it's good to check if there is any weird stuff going on
  pre <- connector_info$partners %like% 'presynaptic_to'
  n_pre <- table(pre)["TRUE"][[1]]
  post <- connector_info$partners %like% 'postsynaptic_to'
  n_post <- table(post)["TRUE"][[1]]
  return(n_post)
}

# synapse polyadicity statistics for a project ---------------------------------
get_syn_polyadicity <- function(pid) {
  pre_connectors <- catmaid_fetch(
    path = paste(pid, "/connectors/", sep = ""),
    body = list(
      relation_type = "presynaptic_to",
      with_partners = "false"
    )
  )

  pre_connector_IDs <- lapply(pre_connectors$connectors,'[[',1)
  
  n_post <- lapply(pre_connector_IDs, n_pre_post, pid)
  
  n_post |> as.numeric() |> summary()
}


# plotting multinucleated cells ------------------------------------------------
soma_positions <- function(nneuron) {
  soma_treenodes <- nneuron$tags$soma
  if (is.null(soma_treenodes)) {
    print("There are no treenodes with soma tag")
    return()
  }
  soma_info <- function(soma_treenode) {
    index <- which(nneuron$d$PointNo == soma_treenode)
    x <- nneuron$d$X[[index]]
    y <- nneuron$d$Y[[index]]
    z <- nneuron$d$Z[[index]]
    radius <- nneuron$d$W[[index]]
    return(data.frame(x=x, y=y, z=z, radius=radius))
  }
  somapos <- ldply(soma_treenodes, soma_info)
  return(somapos)
}

plot_multinucleated_cell <- function(nneuron,
                                     sigma = 1, 
                                     color = "gray",
                                     lwd = 1,
                                     alpha = 1) {
  if ("neuron" %in% class(nneuron)) {
    smoothed_neuron <- smooth_neuron(nneuron, sigma = sigma)
    # plot directly if there is not more than 1 soma
    # not elegant at all, but this way saves a lot of processing time
    if (length(nneuron$tags$soma) < 2) {
      plot3d(smoothed_neuron, lwd = lwd, color = color, alpha = alpha, soma = TRUE)
    } else {
      plot3d(smoothed_neuron, lwd = lwd, color = color, alpha = alpha)
      somapos <- soma_positions(nneuron)
      for (i in seq(nrow(somapos))) {
        spheres3d(x = somapos$x[[i]], 
                  y = somapos$y[[i]], 
                  z = somapos$z[[i]],
                  radius = somapos$radius[[i]]/2,
                  color = color,
                  alpha = alpha)
      }
    }
  } else if ("neuronlist" %in% class(nneuron)) {
    smoothed_neuron <- nlapply(nneuron, smooth_neuron, sigma = 1000)
    plot3d(smoothed_neuron, lwd = lwd, color = color, alpha = alpha)
    for (neu in nneuron) {
      somapos <- soma_positions(neu)
      for (i in seq(nrow(somapos))) {
        spheres3d(x = somapos$x[[i]], 
                  y = somapos$y[[i]], 
                  z = somapos$z[[i]],
                  radius = somapos$radius[[i]]/2,
                  color = color,
                  alpha = alpha)
      }
    }
  }
}


segments_between_tags <- function(nneuron, tag1, tag2) {
  # this assumes we only want the segments between each tag1 to its closest downstream tag2
  segments_neuronlist <- neuronlist()
  ng <- as.ngraph(nneuron)
  tag1_ids <- nneuron$tags[[tag1]]
  tag2_ids <- nneuron$tags[[tag2]]
  tag1_ranks <- match(tag1_ids, nneuron$d$PointNo)
  tag2_ranks <- match(tag2_ids, nneuron$d$PointNo)
  for (tag1_rank in tag1_ranks) {
    proximal_points <- igraph::graph.dfs(ng,
                                         root=tag1_rank, 
                                         unreachable=FALSE,
                                         mode='in')$order
    tag2_matches_indices <- match(tag2_ranks, proximal_points)
    tag2_matches_indices <- tag2_matches_indices[!is.na(tag2_matches_indices)]
    if (length(tag2_matches_indices) == 0) {
      segment_neuronlist <- neuronlist()
    }
    else {
      # get only first match
      tag2_match_index <- sort(tag2_matches_indices)[[1]]
      segment_points <- proximal_points[1:tag2_match_index]
      segment_neuronlist <- subset(nneuron, segment_points)
      #print(xyzmatrix(segment_tree))
      # add information from original neuron
      segment_neuronlist <- c(segment_neuronlist, skid=nneuron$skid)
      # add back tags
      tags_to_add <- list()
      for (i in seq_along(nneuron$tags)) {
        tagname <- names(nneuron$tags)[[i]]
        tag_trids <- nneuron$tags[[i]]
        tag_trids_subset <- intersect(tag_trids, segment_neuronlist$d$PointNo)
        if (length(tag_trids_subset) == 0) {
          next
        }
        tags_to_add <- append(tags_to_add, tag_trids_subset)
        names(tags_to_add)[[length(tags_to_add)]] <- tagname
      }
      if (length(tags_to_add) != 0) {
        segment_neuronlist[[length(segment_neuronlist)+1]] <- as.list(tags_to_add)
        names(segment_neuronlist)[[length(segment_neuronlist)]] <- "tags"
      }
    }

    segments_neuronlist <- c(segments_neuronlist, as.neuronlist(segment_neuronlist))
  }
  return(as.neuronlist(segments_neuronlist))
}

