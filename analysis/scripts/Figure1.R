# Code to generate Figure 1 of the Jokura et al 2024 Ctenophore apical organ connectome paper

# source packages and functions -----------------
source("analysis/scripts/packages_and_functions.R")

# load cell type ----------------------------------

read_smooth_neuron <- function(annotation){
  nlapply(read.neurons.catmaid(annotation, pid = 35),
          function(x)
            smooth_neuron(x, sigma = 1000))
}

balancer <- read_smooth_neuron("celltype:balancer")
bridge <- read_smooth_neuron("celltype:bridge")
bristle <- read_smooth_neuron("celltype:bristle")
dome <- read_smooth_neuron("celltype:dome")
Cgroove_tag <- read_smooth_cell("celltype:Cgroove-tag")
Cgroove_sag <- read_smooth_cell("celltype:Cgroove-sag")
intramulticilia <- read_smooth_neuron("celltype:intra-multi-ciliated")
lamellate <- read_smooth_neuron("celltype:lamellate")
lithocyte <- read_smooth_neuron("celltype:lithocyte")
neuron <- read_smooth_neuron("celltype:neuron")
plumose <- read_smooth_neuron("celltype:plumose")
dense_vesicle <- read_smooth_neuron("celltype:dense_vesicle")
monocilia <- read_smooth_neuron("celltype:monociliated")
bicilia <- read_smooth_neuron("celltype:biciliated")
non_cilia <- read_smooth_neuron("celltype:nonciliated")


Q1 <- nlapply(read.neurons.catmaid("Q1", pid = 35),
              function(x)
                smooth_neuron(x, sigma = 1000))
Q2 <- nlapply(read.neurons.catmaid("Q2", pid = 35),
              function(x)
                smooth_neuron(x, sigma = 1000))
Q3 <- nlapply(read.neurons.catmaid("Q3", pid = 35),
              function(x)
                smooth_neuron(x, sigma = 1000))
Q4 <- nlapply(read.neurons.catmaid("Q4", pid = 35),
              function(x)
                smooth_neuron(x, sigma = 1000))


balancer_Q1 <- nlapply(read.neurons.catmaid(skids_by_2annotations("celltype:balancer", "Q1"),
                                            pid = 35),
                       function(x)
                         smooth_neuron(x, sigma = 1000))
balancer_Q2 <- nlapply(read.neurons.catmaid(skids_by_2annotations("celltype:balancer", "Q2"),
                                            pid = 35),
                       function(x)
                         smooth_neuron(x, sigma = 1000))
balancer_Q3 <- nlapply(read.neurons.catmaid(skids_by_2annotations("celltype:balancer", "Q3"),
                                            pid = 35),
                       function(x)
                         smooth_neuron(x, sigma = 1000))
balancer_Q4 <- nlapply(read.neurons.catmaid(skids_by_2annotations("celltype:balancer", "Q4"),
                                            pid = 35),
                       function(x)
                         smooth_neuron(x, sigma = 1000))


all_celltypes <- list(balancer, bridge, bristle, dome, Cgroove_tag, Cgroove_sag,
                      intramulticilia, lamellate, lithocyte, neuron, 
                      plumose, dense_vesicle, monocilia, bicilia, non_cilia)


# assemble figure -------------------------------------------------------------

#read pics
panel_balancer <- ggdraw() + draw_image(readPNG("manuscript/pictures/balancer.png")) +
  draw_label("balancer cells", x = 0.5, y = 0.95, size = 8.5, fontface="bold", hjust = 0.5) +
  draw_label("aboral view", x = 0.01, y = 0.86, color="black", size = 6, fontface="plain", hjust = 0) +
  draw_label("lateral view of PA plane", x = 0.33, y = 0.86, color="black", size = 6, fontface="plain", hjust = 0) +
  draw_label("lateral view of TA plane", x = 0.69, y = 0.86, color="black", size = 6, fontface="plain", hjust = 0) +
  draw_line(x = c(0.85, 0.95), y = c(0.1, 0.1), color = "black", size = 0.5) +
  draw_label(expression(paste("25 ", mu, "m")), x = 0.9, y = 0.14, color = "black", size = 7, hjust = 0.5) +
  draw_label("TA", x = 0.325, y = 0.16, size = 6, color = "black", hjust = 0.5) +
  geom_segment(aes(x = 0.25, y = 0.16, xend = 0.31, yend = 0.16), data = arrow_TA, color = "black", 
               arrow = arrow(ends = "both", type = "closed", length = unit(0.1,"cm")),
               lineend = "butt",
               linejoin = "mitre",
               arrow.fill = "black", size = 0.175) +
  draw_label("PA", x = 0.28, y = 0.05, size = 6, color = "black", hjust = 0.5) +
  geom_segment(aes(x = 0.28, y = 0.08, xend = 0.28, yend = 0.24), data = arrow_PA, color = "black", 
               arrow = arrow(ends = "both", type = "closed", length = unit(0.1,"cm")),
               lineend = "butt",
               linejoin = "mitre",
               arrow.fill = "black", size = 0.175) +
  draw_label("A", x = 0.66, y = 0.25, size = 6, color = "black", hjust = 0.5) +
  draw_label("O", x = 0.66, y = 0.05, size = 6, color = "black", hjust = 0.5) +
  geom_segment(aes(x = 0.66, y = 0.09, xend = 0.66, yend = 0.21), data = arrow_AO, color = "black", 
               arrow = arrow(ends = "both", type = "closed", length = unit(0.1,"cm")),
               lineend = "butt",
               linejoin = "mitre",
               arrow.fill = "black", size = 0.175)



library(jpeg)
 
panel_larva_5dpf <- ggdraw() + draw_image(readJPEG("manuscript/pictures/Mnemiopsis_larva_5dpf.jpg"))

library(magick)

panel_AO <- ggdraw() + draw_image(image_read("manuscript/pictures/Mnemiopsis_larva_5dpf_AO.jpg"))

panel_AO_schema <- ggdraw() + draw_image(readPNG("manuscript/pictures/Mnemiopsis_larva_5dpf_AO_schematic.png"))
panel_AO_dimention_schema <- ggdraw() + draw_image(readPNG("manuscript/pictures/AO_connectome_paper_Figure1.png"))
panel_catmaid_overview <- ggdraw() + draw_image(readPNG("manuscript/pictures/overview.png"))
panel_catmaid_all_cells <- ggdraw() + draw_image(readPNG("manuscript/pictures/all_cells_3_views.png"))

layout <- "
AABBCCDDD
#########
EEEFFFFFF
"




Figure1 <- panel_larva_5dpf + panel_AO + panel_AO_schema + panel_AO_dimention_schema + 
  panel_catmaid_overview + panel_catmaid_all_cells +
  patchwork::plot_layout(design = layout, 
                         heights = c(1,0.1,1), 
                         widths = c(1,1,1,1,1,1,1)) + 
  patchwork::plot_annotation(tag_levels = "A") &  
  ggplot2::theme(plot.tag = element_text(size = 12, 
                                         face='plain', color='black'))


ggsave("manuscript/figures/Figure1.png", limitsize = FALSE, 
       units = c("px"), Figure1, width = 2400, height = 1000, bg='white')  


ggsave("manuscript/figures/Fig1.pdf", limitsize = FALSE, 
       units = c("px"), Figure1, width = 2400, height = 1000) 






















img_bridge_aboral <- readPNG("manuscript/pictures/bridge_aboral_view.png")
img_bridge_tentacular <- readPNG("manuscript/pictures/bridge_tentacular_plane.png")
img_bridge_sagittal <- readPNG("manuscript/pictures/bridge_sagittal_plane.png")
img_bristle_aboral <- readPNG("manuscript/pictures/bristle_aboral_view.png")
img_bristle_tentacular <- readPNG("manuscript/pictures/bristle_tentacular_plane.png")
img_bristle_sagittal <- readPNG("manuscript/pictures/bristle_sagittal_plane.png")
img_dome_aboral <- readPNG("manuscript/pictures/dome_aboral_view.png")
img_dome_tentacular <- readPNG("manuscript/pictures/dome_tentacular_plane.png")
img_dome_sagittal <- readPNG("manuscript/pictures/dome_sagittal_plane.png")
img_groove_aboral <- readPNG("manuscript/pictures/groove_aboral_view.png")
img_groove_tentacular <- readPNG("manuscript/pictures/groove_tentacular_plane.png")
img_groove_sagittal <- readPNG("manuscript/pictures/groove_sagittal_plane.png")
img_intramulticilia_aboral <- readPNG("manuscript/pictures/intramulticilia_aboral_view.png")
img_intramulticilia_tentacular <- readPNG("manuscript/pictures/intramulticilia_tentacular_plane.png")
img_intramulticilia_sagittal <- readPNG("manuscript/pictures/intramulticilia_sagittal_plane.png")
img_lamellate_aboral <- readPNG("manuscript/pictures/lamellate_aboral_view.png")
img_lamellate_tentacular <- readPNG("manuscript/pictures/lamellate_tentacular_plane.png")
img_lamellate_sagittal <- readPNG("manuscript/pictures/lamellate_sagittal_plane.png")
img_lithocyte_aboral <- readPNG("manuscript/pictures/lithocyte_aboral_view.png")
img_lithocyte_tentacular <- readPNG("manuscript/pictures/lithocyte_tentacular_plane.png")
img_lithocyte_sagittal <- readPNG("manuscript/pictures/lithocyte_sagittal_plane.png")
img_neuron_aboral <- readPNG("manuscript/pictures/neuron_aboral_view.png")
img_neuron_tentacular <- readPNG("manuscript/pictures/neuron_tentacular_plane.png")
img_neuron_sagittal <- readPNG("manuscript/pictures/neuron_sagittal_plane.png")
img_plumose_aboral <- readPNG("manuscript/pictures/plumose_aboral_view.png")
img_plumose_tentacular <- readPNG("manuscript/pictures/plumose_tentacular_plane.png")
img_plumose_sagittal <- readPNG("manuscript/pictures/plumose_sagittal_plane.png")
img_dense_vesicle_aboral <- readPNG("manuscript/pictures/dense_vesicle_aboral_view.png")
img_dense_vesicle_tentacular <- readPNG("manuscript/pictures/dense_vesicle_tentacular_plane.png")
img_dense_vesicle_sagittal <- readPNG("manuscript/pictures/dense_vesicle_sagittal_plane.png")
img_monocilia_aboral <- readPNG("manuscript/pictures/monocilia_aboral_view.png")
img_monocilia_tentacular <- readPNG("manuscript/pictures/monocilia_tentacular_plane.png")
img_monocilia_sagittal <- readPNG("manuscript/pictures/monocilia_sagittal_plane.png")
img_bicilia_aboral <- readPNG("manuscript/pictures/bicilia_aboral_view.png")
img_bicilia_tentacular <- readPNG("manuscript/pictures/bicilia_tentacular_plane.png")
img_bicilia_sagittal <- readPNG("manuscript/pictures/bicilia_sagittal_plane.png")
img_non_cilia_aboral <- readPNG("manuscript/pictures/non_cilia_aboral_view.png")
img_non_cilia_tentacular <- readPNG("manuscript/pictures/non_cilia_tentacular_plane.png")
img_non_cilia_sagittal <- readPNG("manuscript/pictures/non_cilia_sagittal_plane.png")


#make panels

panel_balancer <- ggdraw() + draw_image(img_balancer_aboral) +
  draw_label("aboral view", x = 0.3, y = 0.62, size = 6, fontface = "plain", hjust = 0.5) +
  ggdraw() + draw_image(img_balancer_sagittal) +
  draw_label("saggital plane", x = 0.3, y = 0.62, size = 6, fontface = "plain", hjust = 0.5) + 
  ggdraw() + draw_image(img_balancer_tentacular) +
  draw_label("tentacular plane", x = 0.3, y = 0.62, size = 6, fontface = "plain", hjust = 0.5) +
  draw_line(x = c(0.58, 0.82), y = c(0.395, 0.395), color = "black", size = 0.5) +
  draw_label(expression(paste("20 ", mu, "m")), x = 0.7, y = 0.415, color = "black", size = 6, hjust = 0.5)

panel_bridge <- ggdraw() + draw_image(img_bridge_aboral) +
  ggdraw() + draw_image(img_bridge_sagittal) + 
  ggdraw() + draw_image(img_bridge_tentacular)

panel_bristle <- ggdraw() + draw_image(img_bristle_aboral) +
  ggdraw() + draw_image(img_bristle_sagittal) + 
  ggdraw() + draw_image(img_bristle_tentacular)

panel_dome <- ggdraw() + draw_image(img_dome_aboral) +
  ggdraw() + draw_image(img_dome_sagittal) + 
  ggdraw() + draw_image(img_dome_tentacular)

panel_groove <- ggdraw() + draw_image(img_groove_aboral) +
  ggdraw() + draw_image(img_groove_sagittal) + 
  ggdraw() + draw_image(img_groove_tentacular)

panel_intramulticilia <- ggdraw() + draw_image(img_intramulticilia_aboral) +
  ggdraw() + draw_image(img_intramulticilia_sagittal) + 
  ggdraw() + draw_image(img_intramulticilia_tentacular)

panel_lamellate <- ggdraw() + draw_image(img_lamellate_aboral) +
  ggdraw() + draw_image(img_lamellate_sagittal) + 
  ggdraw() + draw_image(img_lamellate_tentacular)

panel_lithocyte <- ggdraw() + draw_image(img_lithocyte_aboral) +
  ggdraw() + draw_image(img_lithocyte_sagittal) + 
  ggdraw() + draw_image(img_lithocyte_tentacular)

panel_neuron <- ggdraw() + draw_image(img_neuron_aboral) +
  ggdraw() + draw_image(img_neuron_sagittal) + 
  ggdraw() + draw_image(img_neuron_tentacular)

panel_plumose <- ggdraw() + draw_image(img_plumose_aboral) +
  ggdraw() + draw_image(img_plumose_sagittal) + 
  ggdraw() + draw_image(img_plumose_tentacular)

panel_dense_vesicle <- ggdraw() + draw_image(img_dense_vesicle_aboral) +
  ggdraw() + draw_image(img_dense_vesicle_sagittal) + 
  ggdraw() + draw_image(img_dense_vesicle_tentacular)

panel_monocilia <- ggdraw() + draw_image(img_monocilia_aboral) +
  ggdraw() + draw_image(img_monocilia_sagittal) + 
  ggdraw() + draw_image(img_monocilia_tentacular)

panel_bicilia <- ggdraw() + draw_image(img_bicilia_aboral) +
  ggdraw() + draw_image(img_bicilia_sagittal) + 
  ggdraw() + draw_image(img_bicilia_tentacular)

panel_non_cilia <- ggdraw() + draw_image(img_non_cilia_aboral) +
  ggdraw() + draw_image(img_non_cilia_sagittal) + 
  ggdraw() + draw_image(img_non_cilia_tentacular)




# define layout with textual representation for pathchwork assembly of figure
layout <- "
ABC
DEF
GHI
JKL
M##
"

Figure1 <-
  panel_balancer + panel_bridge + panel_bristle + panel_dome + 
  panel_groove + panel_intramulticilia + panel_lamellate + 
  panel_lithocyte + panel_neuron + panel_plumose + 
  panel_dense_vesicle + panel_monocilia + panel_bicilia + panel_non_cilia +
  plot_layout(design = layout) +
  patchwork::plot_layout(design = layout, 
                         heights = c(1,1,1,1,1), 
                         widths = c(1,1,1)) +
  patchwork::plot_annotation(tag_levels = "A") &
  ggplot2::theme(plot.tag = element_text(size = 12, face='plain', color='black'))


ggsave("manuscript/figures/Figure1.png", limitsize = FALSE, 
       units = c("px"), Figure1, width = 3200, height = 3200, bg = "white")


ggsave(
  "manuscript/figures/Figure1.pdf",
  limitsize = FALSE,
  units = c("px"),
  Figure1,
  width = 3300,
  height = 1600
)
