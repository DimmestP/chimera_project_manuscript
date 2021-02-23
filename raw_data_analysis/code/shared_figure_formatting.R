library(ggplot2)
library(cowplot)
library(scales)
library(tidyqpcr)

theme_set(
  theme_cowplot(font_size = 12, 
                font_family = "sans",
                rel_small = 11/12,
                rel_tiny = 10/12,
                rel_large = 14/12) %+replace% 
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
          panel.border=element_rect(colour = "grey50",linetype = "solid",size=0.5),
          panel.grid.major = element_line(size = 0.15, linetype = 'solid',colour = "grey"),
          strip.background = element_blank()
    )
)

RPS3_TSA1_colour_scheme <- c( "grey50", "#288a2e","#a84a9a", 
                              "#6f3ba1","#416db0","#CC6666","black")

PIR1_colour_scheme <- c("#971c9c","#1e1c9c","#5d1c9c",
                        "#a84a9a","#9884ab","#416db0",
                        "#CC6666","black")

RNA_relative_abundance_figure_options <- list(
  geom_point(aes(rel_abund_delta_deltacq,UTR3,colour=UTR3), size = 2),
    scale_x_log2nice(omag = seq(-5,5),scilabels=FALSE),
    guides(colour=FALSE),
    theme(axis.text.x=element_text(angle=0,vjust=0.5),
          axis.title.x=element_text(vjust=-2),
          legend.position="bottom",
          legend.box.margin=margin(20,10,10,180)),
  stat_summary(aes(x=rel_abund_delta_deltacq,y=UTR3),
               fun="mean",colour="black",
               geom="crossbar"))

protein_raw_abundance_figure_options <- list(
  geom_point(aes(y=Terminator,x=fluo_per_OD_at_max_gr, colour = Terminator)),
  scale_colour_hue(h = c(0, 360)+20,l=60,c=60),
  geom_vline(xintercept=1, linetype="dotted", color = "grey50"),
  stat_summary(aes(y=Terminator,x=fluo_per_OD_at_max_gr),
               fun="mean",colour="black",
               geom="crossbar",size=0.2, width=0.5),
  theme(legend.position = "none"),
  scale_x_continuous(oob=oob_squish, limits = c(0,NA))
)

protein_relative_abundance_figure_options <- list(
  scale_colour_hue(h = c(0, 360)+20,l=60,c=60),
  geom_vline(xintercept=1, linetype="dotted", color = "grey50"),
  stat_summary(aes(y=Terminator,x=fluo_per_OD_at_max_gr, colour = Terminator),
               fun.data="mean_se",
               geom="errorbarh",
               size=1),
  stat_summary(aes(y=Terminator,x=fluo_per_OD_at_max_gr, colour = Terminator),
               fun="mean",
               geom="point",
               size = 2),
  theme(legend.position = "none"),
  scale_x_continuous(oob=oob_squish, limits = c(0,1.9))
)

update_geom_defaults("point", list(size = 3))
update_geom_defaults("errorbar", list(size = 1))
update_geom_defaults("errorbarh", list(size = 1))
update_geom_defaults("crossbar", list(size=0.5, width=1))


