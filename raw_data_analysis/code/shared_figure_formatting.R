library(ggplot2)
library(cowplot)
library(scales)
library(tidyqpcr)

# set default settings for plots

theme_set(
  theme_cowplot(font_size = 12, 
                font_family = "sans",
                rel_small = 9/12,
                rel_tiny = 8/12,
                rel_large = 14/12) %+replace% 
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
          panel.border=element_rect(colour = "grey50",linetype = "solid",size=0.5),
          panel.grid.major = element_line(size = 0.15, linetype = 'solid',colour = "grey"),
          strip.background = element_blank()
    )
)



update_geom_defaults("point", list(size = 2))
update_geom_defaults("errorbar", list(size = 1))
update_geom_defaults("errorbarh", list(size = 1))

# Set general variables for motif construct names and colour schemes

construct_to_label_dictionary_TSA1_RPS3 <-  
  tibble(construct = c("WT","modC","modE","modD","modA","modB","mod0"), 
         label = c("WT", "mod_NGG", "mod_HTH", "mod_HNH", "mod_NTN", "mod_NAA", "mod_NNN"))

construct_to_label_dictionary_PIR1 <- 
  tibble(construct = c("modG", "modF", "modE", "modD", "modC", "modA", "modB", "WT"), 
         label = c("mod_NTNNN", "mod_ANNNN", "mod_ATNNN", "mod_ATHNH", "mod_ATNHH", "mod_ANHHH", "mod_NTHHH", "WT"))

RPS3_TSA1_colour_scheme <- c( "#6c6c6c", "#7FBF74","#C288BB", 
                              "#71519C","#5978BA","#C05558","black")

PIR1_colour_scheme <- c("#DD76A5","#968BC2","#71519C",
                        "#D97F1D","#4EB0B5","#5978BA",
                        "#C15659","#6c6c6c")

# Create plotting functions for  specific data sets

RNA_relative_abundance_figure_options <- list(
  geom_point(aes(x=rel_abund_delta_deltacq,y=UTR3,colour=UTR3),
             shape = 18, size = 2),
  scale_x_log2nice(omag = seq(-5,5),scilabels=FALSE),
  guides(colour=FALSE),
  geom_vline(xintercept=1, linetype="dotted", color = "grey50"),
  theme(axis.text.x=element_text(angle=0,vjust=0.5),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 10),
        legend.position="bottom",
        legend.box.margin=margin(20,10,10,180)),
  stat_summary(aes(x=rel_abund_delta_deltacq,y=UTR3),
               fun="mean",colour="black",
               geom="crossbar", size=0.3, width=0.9)
  )

protein_raw_abundance_figure_options <- list(
  geom_point(aes(y=Terminator,x=fluo_per_OD_at_max_gr, colour = Terminator),
             shape = 18, size = 2),
  scale_colour_hue(h = c(0, 360)+20,l=60,c=60),
  geom_vline(xintercept=1, linetype="dotted", color = "grey50"),
  stat_summary(aes(y=Terminator,x=fluo_per_OD_at_max_gr),
               fun="mean",colour="black",
               geom="crossbar", size=0.3, width=0.9),
  theme(legend.position = "none",
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 10)),
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
  theme(legend.position = "none",
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 10)),
  scale_x_continuous(oob=oob_squish, limits = c(0,1.9))
)

protein_vs_RNA_figure_options <- list(
    geom_abline(slope = 1, intercept = 0),
    geom_point(aes(y = mean_relative_abundance_protein, x = mean_relative_abundance_mrna, colour = label)),
    geom_errorbar(aes(ymax = mean_relative_abundance_protein + se_relative_abundance_protein,
                      ymin = mean_relative_abundance_protein - se_relative_abundance_protein,
                      x = mean_relative_abundance_mrna, 
                      colour = label)),
    geom_errorbarh(aes(xmax = mean_relative_abundance_mrna + se_relative_abundance_mrna,
                       xmin = mean_relative_abundance_mrna - se_relative_abundance_mrna,
                       y = mean_relative_abundance_protein, 
                       colour = label)),
    facet_wrap(~ promoter, ncol = 1),
    scale_x_log2nice(),
    scale_y_log2nice()
)

