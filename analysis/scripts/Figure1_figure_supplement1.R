# Code to generate Figure 1 Supplements of the Jokura et al 2024 Ctenophore apical organ connectome paper

# source packages and functions ------------------------------------------------
source("analysis/scripts/packages_and_functions.R")


# show just a subset of celltypes, just the interesting ones

celltype_map <- c(
  "balancer" = "balancer cells",
  "bridge" = "bridge cells",
  "large_granular_cell" = "large granular cells",
  "Cgroove" = "ciliated groove cells",
  "dense_vesicle" = "dense vesicle cells",
  "dome" = "dome cells",
  "intra-multi-ciliated" = "intracellular multiciliated cells",
  "lamellate" = "lamellate bodies",
  "lithocyte" = "lithocytes",
  "plumose" = "plumose cells",
  "SSN" = "aboral nerve net neurons",
  "epithelial_floor" = "epithelial floor cells"
)

# celltypes short, because these are not all celltypes in the animal
celltypes_short <- names(celltype_map)


skids_with_soma <- catmaid_get_label_stats(pid=35) %>%
  filter(labelName == "soma") %>%
  select(skeletonID) %>%
  pull() %>%
  unique()

skids_outside <- catmaid_skids("outside", pid=35)

skids_soma_AO <- setdiff(skids_with_soma, skids_outside)

with_soma <- read_smooth_neuron(skids_soma_AO)

celltype_cells <- list()
for (celltype in celltypes_short) {
  nneuron <- read_smooth_neuron(paste("celltype:", celltype, sep=""))
  #assign(celltype, nneuron, envir = .GlobalEnv)
  celltype_cells[[celltype]]=nneuron
}


# plot cells -----
for (i in seq_along(celltype_cells)) {
  celltype <- names(celltype_cells)[i]
  cells <- celltype_cells[[i]]
  close3d()
  nopen3d()
  mfrow3d(1, 3)
  par3d(windowRect = c(0, 0, 2400, 700))
  
  for (view in c("aboral", "sagittal", "tentacular")) {
    # plot3d doesn't plot all nuclei in multinucleated cells
    # we have to use plot_multinucleated cell instead
    # but plot_multinucleated_cell is much slower, so we want to use it only on cells which are actually multinucleated
    if (celltype == "SSN") {
      plot_multinucleated_cell(cells, lwd = 1, alpha = 0.8, col = Okabe_Ito[5])
    } else {
      plot3d(cells, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.8, col = Okabe_Ito[5], WithConnectors = FALSE, WithNodes = FALSE)
    }
    plot3d(with_soma, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.05, col = Okabe_Ito[8])
    plot3d(outline, add = TRUE, alpha = 0.025, col = Okabe_Ito[8])
    
    switch(view,
           "aboral" = aboral(),
           "sagittal" = sagittal(),
           "tentacular" = tentacular())
    
    par3d(zoom = 0.61)
    if (view != "tentacular") next3d(clear = FALSE)
  }
  
  rgl.snapshot(paste0("manuscript/pictures/3d_plot/plot_", celltype, ".png"))
  close3d()
}


# Display all cell types together in one 3d plot--------------------------------
# get all cell types 

celltypes <- get_celltypes(35)
celltypes_list <- list()

for (celltype in celltypes) {
  celltypes_list[[celltype]] <- read.neurons.catmaid(paste("celltype:", celltype, sep = ""), pid = 35)
}

print(celltypes)

# Define colors

display.brewer.pal(5, "Blues")
blues <- brewer.pal(5, "Blues")

display.brewer.pal(5, "BuPu")
bluepurple <- brewer.pal(5, "BuPu")

colour_palettes <- c(Okabe_Ito[-8], bluepurple, blues)

length(unique(colour_palettes))


# 3d plots of three planes for all cell types

plot_celltypes_views <- function(celltypes_list, colour_palettes, output_file) {
  nopen3d()
  mfrow3d(1, 3)
  par3d(windowRect = c(0, 0, 2400, 700))
  
  views <- c("aboral", "sagittal", "tentacular")
  
  for (view in views) {
    start.time <- Sys.time()
    
    for (i in seq_along(celltypes_list)) {
      if (names(celltypes_list)[[i]] %in% c("SSN", "SNN")) {
        plot_multinucleated_cell(celltypes_list[[i]], alpha = 0.6, color = colour_palettes[i])
      } else {
        plot3d(celltypes_list[[i]], alpha = 0.6, color = colour_palettes[i], soma = TRUE)
      }
    }
    
    end.time <- Sys.time()
    print(end.time - start.time)
    
    switch(view,
           "aboral" = aboral(),
           "sagittal" = sagittal(),
           "tentacular" = tentacular())
    
    par3d(zoom = 0.61)
    
    if (view != "tentacular") next3d(clear = FALSE)
  }
  
  rgl.snapshot(output_file)
  close3d()
}


plot_celltypes_views(celltypes_list, colour_palettes, "manuscript/pictures/3d_plot/all_cells_3_views_alt.png")


# assemble figure -------------------------------------------------------------


#Figure1 Supplement1
#read pics
panel_balancer <- ggdraw() + draw_image(readPNG("manuscript/pictures/3d_plot/plot_balancer.png")) +
  draw_label("balancer cells", x = 0.5, y = 0.95, size = 8.5, fontface="bold", hjust = 0.5) +
  draw_label("aboral view", x = 0.01, y = 0.86, color="black", size = 6, fontface="plain", hjust = 0) +
  draw_label("sagittal plane", x = 0.33, y = 0.86, color="black", size = 6, fontface="plain", hjust = 0) +
  draw_label("tentacular plane", x = 0.69, y = 0.86, color="black", size = 6, fontface="plain", hjust = 0) +
  draw_line(x = c(0.85, 0.95), y = c(0.1, 0.1), color = "black", size = 0.5) +
  draw_label(expression(paste("25 ", mu, "m")), x = 0.9, y = 0.14, color = "black", size = 7, hjust = 0.5) +
  draw_label("TA", x = 0.325, y = 0.16, size = 6, color = "black", hjust = 0.5) +
  geom_segment(aes(x = 0.25, y = 0.16, xend = 0.31, yend = 0.16), color = "black", 
               arrow = arrow(ends = "both", type = "closed", length = unit(0.1,"cm")),
               lineend = "butt",
               linejoin = "mitre",
               arrow.fill = "black", linewidth = 0.175) +
  draw_label("SA", x = 0.28, y = 0.05, size = 6, color = "black", hjust = 0.5) +
  geom_segment(aes(x = 0.28, y = 0.08, xend = 0.28, yend = 0.24), color = "black", 
               arrow = arrow(ends = "both", type = "closed", length = unit(0.1,"cm")),
               lineend = "butt",
               linejoin = "mitre",
               arrow.fill = "black", linewidth = 0.175) +
  draw_label("A", x = 0.66, y = 0.25, size = 6, color = "black", hjust = 0.5) +
  draw_label("O", x = 0.66, y = 0.05, size = 6, color = "black", hjust = 0.5) +
  geom_segment(aes(x = 0.66, y = 0.09, xend = 0.66, yend = 0.21), color = "black", 
               arrow = arrow(ends = "both", type = "closed", length = unit(0.1,"cm")),
               lineend = "butt",
               linejoin = "mitre",
               arrow.fill = "black", linewidth = 0.175)

panel_celltype <- list()
for (celltype in names(celltype_cells)) {
  name_text <- rename_map[celltype]
  panel_celltype[[celltype]] <- ggdraw() +
    draw_image(readPNG(paste("manuscript/pictures/3d_plot/plot_", celltype, ".png", sep = ""))) +
    draw_label(name_text, x = 0.5, y = 0.95, size = 8, fontface="bold", hjust = 0.5)
}


layout <- "
A#B#C
#####
D#E#F
#####
G#H#I
#####
J#K#L
"
# panel_balancer is different, because we draw coodrinates and scale bar above
Fig1_Sup1 <- panel_balancer + panel_celltype[["bridge"]] + panel_celltype[["large_granular_cell"]] +
  panel_celltype[["Cgroove"]] + panel_celltype[["dense_vesicle"]] + panel_celltype[["dome"]] +
  panel_celltype[["intra-multi-ciliated"]] + panel_celltype[["lamellate"]] + panel_celltype[["lithocyte"]] +
  panel_celltype[["plumose"]] + panel_celltype[["SSN"]] + panel_celltype[["epithelial_floor"]] +
  patchwork::plot_layout(design = layout, 
                         heights = c(1,0.05,1,0.05,1,0.05,1,0.05,1), 
                         widths = c(1,0.05,1,0.05,1)) + 
  patchwork::plot_annotation(tag_levels = "A") &  
  ggplot2::theme(plot.tag = element_text(size = 12, 
                                         face='plain', color='black'))

ggsave("manuscript/figures/Figure1_Supplement1.png", limitsize = FALSE, 
       units = c("px"), Fig1_Sup1, width = 3800, height = 2000, bg='white') 

ggsave("manuscript/figures/Figure1_Supplement1.pdf", limitsize = FALSE, 
      units = c("px"), Fig1_Sup1, width = 3800, height = 2000) 



#Figure1 Supplement2

panel_all_cells <- ggdraw() + draw_image(readPNG("manuscript/pictures/3d_plot/all_cells_3_views_alt.png")) +
  draw_label("all cells", x = 0.5, y = 0.95, size = 8.5, fontface="bold", hjust = 0.5) +
  draw_label("aboral view", x = 0.01, y = 0.86, color="black", size = 6, fontface="plain", hjust = 0) +
  draw_label("sagittal plane", x = 0.33, y = 0.86, color="black", size = 6, fontface="plain", hjust = 0) +
  draw_label("tentacular plane", x = 0.69, y = 0.86, color="black", size = 6, fontface="plain", hjust = 0) +
  draw_line(x = c(0.85, 0.95), y = c(0.1, 0.1), color = "black", size = 0.5) +
  draw_label(expression(paste("25 ", mu, "m")), x = 0.9, y = 0.14, color = "black", size = 7, hjust = 0.5) +
  draw_label("TA", x = 0.325, y = 0.16, size = 6, color = "black", hjust = 0.5) +
  geom_segment(aes(x = 0.25, y = 0.16, xend = 0.31, yend = 0.16), color = "black", 
               arrow = arrow(ends = "both", type = "closed", length = unit(0.1,"cm")),
               lineend = "butt",
               linejoin = "mitre",
               arrow.fill = "black", linewidth = 0.175) +
  draw_label("SA", x = 0.28, y = 0.05, size = 6, color = "black", hjust = 0.5) +
  geom_segment(aes(x = 0.28, y = 0.08, xend = 0.28, yend = 0.24), color = "black", 
               arrow = arrow(ends = "both", type = "closed", length = unit(0.1,"cm")),
               lineend = "butt",
               linejoin = "mitre",
               arrow.fill = "black", linewidth = 0.175) +
  draw_label("A", x = 0.66, y = 0.25, size = 6, color = "black", hjust = 0.5) +
  draw_label("O", x = 0.66, y = 0.05, size = 6, color = "black", hjust = 0.5) +
  geom_segment(aes(x = 0.66, y = 0.09, xend = 0.66, yend = 0.21), color = "black", 
               arrow = arrow(ends = "both", type = "closed", length = unit(0.1,"cm")),
               lineend = "butt",
               linejoin = "mitre",
               arrow.fill = "black", linewidth = 0.175)

layout <- "
A
"

Fig1_Sup2 <- panel_all_cells +
  plot_layout(design = layout,
              heights = c(),
              widths = c()) + 
  patchwork::plot_annotation(tag_levels = "A") &  
  ggplot2::theme(plot.tag = element_text(size = 12, 
                                         face='plain', color='black'))

ggsave("manuscript/figures/Figure1_Supplement2.png", limitsize = FALSE, 
       units = c("px"), Fig1_Sup2, width = 1200, height = 500, bg='white')  


ggsave("manuscript/figures/Figure1_Supplement2.pdf", limitsize = FALSE, 
       units = c("px"), Fig1_Sup2, width = 1200, height = 500) 



