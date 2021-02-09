library(ggplot2)
library(cowplot)
# library(extrafont)

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

update_geom_defaults("point", list(size = 3))
update_geom_defaults("errorbar", list(size = 1))
update_geom_defaults("errorbarh", list(size = 1))

