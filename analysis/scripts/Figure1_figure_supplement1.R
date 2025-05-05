# Code to generate Figure 1 Supplements of the Jokura et al 2024 Ctenophore apical organ connectome paper

# source packages and functions ------------------------------------------------
source("analysis/scripts/packages_and_functions.R")


# load all cell types to plot individually -------------------------------------
#twelve celltypes
balancer <- read_smooth_neuron("celltype:balancer")
bridge <- read_smooth_neuron("celltype:bridge")
bristle <- read_smooth_neuron("celltype:bristle")
Cgroove <- read_smooth_neuron("celltype:Cgroove")
dense_vesicle <- read_smooth_neuron("celltype:dense_vesicle")
dome <- read_smooth_neuron("celltype:dome")
intra_multi_ciliated <- read_smooth_neuron("celltype:intra-multi-ciliated")
lamellate <- read_smooth_neuron("celltype:lamellate")
lithocyte <- read_smooth_neuron("celltype:lithocyte")
plumose <- read_smooth_neuron("celltype:plumose")
SSN <- read_smooth_neuron("celltype:SSN")
epithelial_floor <- read_smooth_neuron("celltype:epithelial_floor")

with_soma <- read_smooth_neuron("with_soma")

# function individual plot process --------------------------------------------------------

plot_neuron_views <- function(neuron_name, neuron_data) {
  close3d()
  nopen3d()
  mfrow3d(1, 3)
  par3d(windowRect = c(0, 0, 2400, 700))
  
  for (view in c("aboral", "sagittal", "tentacular")) {
    plot3d(neuron_data, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.8, col = Okabe_Ito[5], WithConnectors = FALSE, WithNodes = FALSE)
    plot3d(with_soma, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.05, col = Okabe_Ito[8])
    plot3d(outline, add = TRUE, alpha = 0.025, col = Okabe_Ito[8])
    
    switch(view,
           "aboral" = aboral(),
           "sagittal" = sagittal(),
           "tentacular" = tentacular())
    
    par3d(zoom = 0.61)
    if (view != "tentacular") next3d(clear = FALSE)
  }
  
  rgl.snapshot(paste0("manuscript/pictures/3d_plot/plot_", neuron_name, ".png"))
  close3d()
}

# individual plot each cells except SSNs---------------------------------------------------

plot_neuron_views("balancer", balancer)
plot_neuron_views("bridge", bridge)
plot_neuron_views("bristle", bristle)
plot_neuron_views("Cgroove", Cgroove)
plot_neuron_views("dense_vesicle", dense_vesicle)
plot_neuron_views("dome", dome)
plot_neuron_views("intra_multi_ciliated", intra_multi_ciliated)
plot_neuron_views("lamellate", lamellate)
plot_neuron_views("lithocyte", lithocyte)
plot_neuron_views("plumose", plumose)
plot_neuron_views("epithelial_floor", epithelial_floor)

# plot SSNs --------------------------------------------------------------------
plot_multinucleated_views <- function(neuron_name, neuron_data) {
  close3d()
  nopen3d()
  mfrow3d(1, 3)
  par3d(windowRect = c(0, 0, 2400, 700))
  
  views <- c("aboral", "sagittal", "tentacular")
  
  for (view in views) {
    plot_multinucleated_cell(neuron_data, lwd = 1, alpha = 0.8, col = Okabe_Ito[5])
    plot3d(with_soma, soma = TRUE, lwd = 1, add = TRUE, alpha = 0.05, col = Okabe_Ito[8])
    plot3d(outline, add = TRUE, alpha = 0.025, col = Okabe_Ito[8])
    
    switch(view,
           "aboral" = aboral(),
           "sagittal" = sagittal(),
           "tentacular" = tentacular())
    
    par3d(zoom = 0.61)
    
    if (view != "tentacular") next3d(clear = FALSE)
  }
  
  rgl.snapshot(paste0("manuscript/pictures/3d_plot/plot_", neuron_name, ".png"))
  close3d()
}

plot_multinucleated_views("SSN", SSN)


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

panel_bridge <- ggdraw() + draw_image(readPNG("manuscript/pictures/3d_plot/plot_bridge.png")) +
  draw_label("bridge cells", x = 0.5, y = 0.95, size = 8, fontface="bold", hjust = 0.5)
panel_bristle <- ggdraw() + draw_image(readPNG("manuscript/pictures/3d_plot/plot_bristle.png")) +
  draw_label("bristle cells", x = 0.5, y = 0.95, size = 8, fontface="bold", hjust = 0.5)
panel_Cgroove <- ggdraw() + draw_image(readPNG("manuscript/pictures/3d_plot/plot_Cgroove.png")) +
  draw_label("ciliated groove cells", x = 0.5, y = 0.95, size = 8, fontface="bold", hjust = 0.5)
panel_dense_vesicle <- ggdraw() + draw_image(readPNG("manuscript/pictures/3d_plot/plot_dense_vesicle.png")) +
  draw_label("dense vesicle cells", x = 0.5, y = 0.95, size = 8, fontface="bold", hjust = 0.5)
panel_dome <- ggdraw() + draw_image(readPNG("manuscript/pictures/3d_plot/plot_dome.png")) +
  draw_label("dome cells", x = 0.5, y = 0.95, size = 8, fontface="bold", hjust = 0.5)
panel_intra_multi_ciliated <- ggdraw() + draw_image(readPNG("manuscript/pictures/3d_plot/plot_intra_multi_ciliated.png")) +
  draw_label("intra-multiciliated cells", x = 0.5, y = 0.95, size = 8, fontface="bold", hjust = 0.5)
panel_lamellate <- ggdraw() + draw_image(readPNG("manuscript/pictures/3d_plot/plot_lamellate.png")) +
  draw_label("lamellate bodies", x = 0.5, y = 0.95, size = 8, fontface="bold", hjust = 0.5)
panel_lithocyte <- ggdraw() + draw_image(readPNG("manuscript/pictures/3d_plot/plot_lithocyte.png")) +
  draw_label("lithocytes", x = 0.5, y = 0.95, size = 8, fontface="bold", hjust = 0.5)
panel_plumose <- ggdraw() + draw_image(readPNG("manuscript/pictures/3d_plot/plot_plumose.png")) +
  draw_label("plumose cells", x = 0.5, y = 0.95, size = 8, fontface="bold", hjust = 0.5)
panel_SSN <- ggdraw() + draw_image(readPNG("manuscript/pictures/3d_plot/plot_SSN.png")) +
  draw_label("aboral nerve net neurons", x = 0.5, y = 0.95, size = 8, fontface="bold", hjust = 0.5)
panel_epithelial_floor <- ggdraw() + draw_image(readPNG("manuscript/pictures/3d_plot/plot_epithelial_floor.png")) +
  draw_label("epithelial floor cells", x = 0.5, y = 0.95, size = 8, fontface="bold", hjust = 0.5)

layout <- "
A#B#C
#####
D#E#F
#####
G#H#I
#####
J#K#L
"

Fig1_Sup1 <- panel_balancer + panel_bridge + panel_bristle +
  panel_Cgroove + panel_dense_vesicle + panel_dome + 
  panel_intra_multi_ciliated + panel_lamellate + panel_lithocyte + 
  panel_plumose + panel_SSN + panel_epithelial_floor +
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



