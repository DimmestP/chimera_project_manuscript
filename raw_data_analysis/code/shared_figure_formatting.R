theme_set(theme_cowplot(font_size=11) %+replace% 
            theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
                  panel.border=element_rect(colour = "grey50",linetype = "solid",size=0.5),
                  panel.grid.major = element_line(size = 0.15, linetype = 'solid',colour = "grey"),
                  plot.title = element_text(family="Tahoma",size=18),
                  strip.background = element_blank(),
                  strip.text.x=element_text(family="Tahoma",size=11),
                  axis.text=element_text(family="Tahoma", size=10),
                  axis.title=element_text(family="Tahoma", size=15),
                  legend.text=element_text(family="Tahoma", size=10),
                  legend.title=element_text(family="Tahoma", size=12)))

update_geom_defaults("point", list(size = 3))
update_geom_defaults("errorbar", list(size = 1))
update_geom_defaults("errorbarh", list(size = 1))

