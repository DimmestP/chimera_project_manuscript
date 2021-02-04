library(ggplot2)
library(cowplot)
library(extrafont)

theme_set(theme_cowplot(font_size=11) %+replace% 
            theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
                  panel.border=element_rect(colour = "grey50",linetype = "solid",size=0.5),
                  panel.grid.major = element_line(size = 0.15, linetype = 'solid',colour = "grey"),
                  plot.title = element_text(family="FreeSans",size=18),
                  strip.background = element_blank(),
                  strip.text.x=element_text(family="FreeSans",size=11),
                  axis.text=element_text(family="FreeSans", size=10),
                  axis.title=element_text(family="FreeSans", size=15),
                  legend.text=element_text(family="FreeSans", size=10),
                  legend.title=element_text(family="FreeSans", size=12)))

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
  stat_summary(aes(rel_abund_delta_deltacq,UTR3),
               fun="mean",colour="black",
               geom="crossbar"))

update_geom_defaults("point", list(size = 3))
update_geom_defaults("errorbar", list(size = 1))
update_geom_defaults("errorbarh", list(size = 1))
update_geom_defaults("crossbar", list(size=0.5, width=1))


