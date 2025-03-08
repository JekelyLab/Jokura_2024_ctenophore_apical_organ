# Code to generate Figure 4 of the Jokura et al 2024 Ctenophore apical organ connectome paper

# source packages and functions ------------------------------------------------
source("analysis/scripts/packages_and_functions.R")


# read data---------------------------------------------------------------------

df_time_dif <- read_csv("analysis/data/arrest_rebeat_time_differences.csv")


# statistics analysis of the L-NAME tracking mean vertical position (side UV after 30sec)--------------------------------------

#Change the name to perform the Dunnett
#df_L_NAME_stat <- df_L_NAME_tracking %>%
#  filter(sec >= 39.92 & sec < 40)
#
##Dunnett's test
#L_NAME = factor(df_L_NAME_stat$'L-NAME')
#y_axis = df_L_NAME_stat$y_axis
#summary(glht(aov(y_axis~L_NAME),linfct=mcp(L_NAME = "Dunnett")))
#
##                   Estimate Std. Error t value Pr(>|t|)    
##0.1 mM - 0 mM == 0   13.805      8.683   1.590 0.223749    
##1 mM - 0 mM == 0     39.897      8.683   4.595 0.000787 ***
#
##‘***’ <0.001 ‘**’ <0.01 ‘*’ <0.05

install.packages("ggbeeswarm")

library("ggbeeswarm")

plot_arrest_rebeat <- 
  ggplot(df_time_dif) +
  aes(x = interaction(plane, response), y = time, fill = response) +
  geom_boxplot(size = 0.5, outlier.shape = NA) +
  geom_beeswarm(aes(color = time),
                size = 2,
                cex = 3.2,
                alpha = 0.3,
                color = "gray5") +
  labs(y = "time difference (sec)") +
  theme_minimal() +
  theme(legend.position = 'none',
        axis.title.x = element_text(size = 9)) +
  scale_fill_manual(
    values = c("re-beat" = "deepskyblue2", 
               "arrest" = "magenta"))


plot_arrest_rebeat

#  scale_x_discrete(labels = c("0 mM" = "0 mM", 
#                              "0.1 mM" = "0.1 mM",
#                              "1 mM" = "1.0 mM")) +
#  geom_signif(comparisons = list(c("0 mM", "0.1 mM")),
#              annotations = "0.22",
#              y_position = 30) +
#  geom_signif(comparisons = list(c("0 mM", "1 mM")),
#              annotations = "0.0008",
#              y_position = 40)




ggsave(
  filename = "manuscript/pictures/arrest_rebeat_graph.png",
  plot = plot_arrest_rebeat,
  width = 600,
  height = 700,
  units = "px",
  dpi = 300
)


# assemble figure -------------------------------------------------------------

panel_ms <- ggdraw() + draw_image(readPNG("manuscript/pictures/tilt_microscope.png"))

panel_kymograph <- ggdraw() + draw_image(readPNG("manuscript/pictures/balacer_arrest_rebeat_kymograph.png"))

panel_graph <- ggdraw() + draw_image(readPNG("manuscript/pictures/arrest_rebeat_graph.png"))

panel_comparison <- ggdraw() + draw_image(readPNG("manuscript/pictures/comparison.png"))

layout <- "
AAA#B#C
AAA#DDD
AAA#DDD
#######
EEFFFFF
"



Figure4 <- panel_ms + panel_kymograph + 
  panel_graph + panel_comparison +
  plot_layout(design = layout,
              heights = c(1, 1, 1, 0.1, 3),
              widths = c(1, 1, 1, 0.1, 1, 0.1, 1)) + 
  patchwork::plot_annotation(tag_levels = "A") &  
  ggplot2::theme(plot.tag = element_text(size = 12, 
                                         face='plain', color='black'))



ggsave("manuscript/figures/Figure4.png", limitsize = FALSE, 
       units = c("px"), Figure4, width = 2400, height = 1700, bg='white')  


ggsave("manuscript/figures/Figure4.pdf", limitsize = FALSE, 
       units = c("px"), Figure4, width = 2400, height = 1700) 


