library(magick)
library(patchwork)
library(gridExtra)
library(here)
library(dplyr)
library(readr)
library(stringr)

source(here("raw_data_analysis/code/shared_figure_formatting.R"))

### TSA1 and RPS3 insertion diagram

layout <- "AAA
           BCC"


construct_to_label_dictionary_TSA1_RPS3 = tibble(construct = c("WT","modC","modE","modD","modA","modB","mod0"), 
                                       label = c("WT", "mod_NGG", "mod_HTH", "mod_HNH", "mod_NTN", "mod_NAA", "mod_NNN"),
                                       TSA1_UTR3 = str_c("tTSA1_",construct),
                                       RPS3_UTR3 = str_c("tRPS3_",construct))

TSA1_deltadeltacq <- read_csv(here("raw_data_analysis/data/norm_qpcr/motif_context_dependence/pTSA1_pPGK1_pSRO9_tTSA1_deltadeltacq_platesnorm_summarise.csv")) %>%
  mutate(UTR3 = factor(UTR3,levels = construct_to_label_dictionary_TSA1_RPS3$TSA1_UTR3),
         promoter = factor(promoter, levels = c("pTSA1","pPGK1","pSRO9")))

RPS3_deltadeltacq <- read_csv(here("raw_data_analysis/data/norm_qpcr/motif_context_dependence/pRPS3_pPGK1_pSRO9_tRPS3_deltadeltacq_platesnorm_summarise.csv")) %>%
  mutate(UTR3 = factor(UTR3,levels = construct_to_label_dictionary_TSA1_RPS3$RPS3_UTR3),
         promoter = factor(promoter, levels = c("pRPS3","pPGK1","pSRO9")))

TSA1_deltadeltacq_plot <- ggplot(data = TSA1_deltadeltacq)+
  RNA_relative_abundance_figure_options +
  scale_colour_manual(values=RPS3_TSA1_colour_scheme) +
  theme(axis.text.y=element_text(colour=RPS3_TSA1_colour_scheme)) +
  labs(x ="Fold change in RNA abundance relative to tTSA1 mod_NNN (log2 scale)", y = "") +
  facet_wrap(~promoter,ncol = 3) +
  scale_y_discrete(breaks = construct_to_label_dictionary_TSA1_RPS3$TSA1_UTR3, labels = construct_to_label_dictionary_TSA1_RPS3$label)

RPS3_deltadeltacq_plot <- ggplot(data = RPS3_deltadeltacq) +
  RNA_relative_abundance_figure_options +
  scale_colour_manual(values=RPS3_TSA1_colour_scheme) +
  theme(axis.text.y=element_text(colour=RPS3_TSA1_colour_scheme)) +
  labs(x ="Fold change in RNA abundance relative to tRPS3 mod_NNN (log2 scale)", y = "") +
  facet_wrap(~promoter,ncol = 3) +
  scale_y_discrete(breaks = construct_to_label_dictionary_TSA1_RPS3$RPS3_UTR3, labels = construct_to_label_dictionary_TSA1_RPS3$label)

insertion_layout <- "AC
                     AC
                     BC"

## Can we modularize this plot to produce the tableGrob in one place, then call it from the draw?
insertion_construct_plot <- wrap_elements(full = (wrap_elements(full = (ggdraw() + 
                                                    draw_image(image_read_svg(here("raw_data_analysis/figures/terminator_construct_designs/RPS3_TSA1_motif_mod0_construct_design.svg"))))) +
  wrap_elements(full = tableGrob(tibble("Construct" = c("mod_NNN", "mod_NAA", "mod_NTN",
                                                "mod_HNH", "mod_HTH", "mod_NGG", "WT"), 
                                 "Motif 1" = c("Random", "Random", "Random", 
                                               "HWNCATTWY", "HWNCATTWY", "Random", "NA"), 
                                 "Motif 2" = c("Random", "ATATTC", "TGTAHMNTA", 
                                               "Random", "TGTAHMNTA", "GTATACCTA", "NA"), 
                                 "Motif 3" = c("Random", "ATATTC", "Random", 
                                               "HWNCATTWY", "HWNCATTWY", "GTATACCTA", "NA")),
                          rows = NULL,
                          theme = ttheme_default(base_size=5,
                                                 padding = unit(c(2, 2), "mm"),
                                                 core=list(
                                                   fg_params=list(col = c("black", "#C05558", "#5978BA", "#71519C", "#C288BB", "#7FBF74", "#878787",
                                                                          "black", "black", "black", "#71519C", "#C288BB", "black", "#878787",
                                                                          "black", "#C05558", "#5978BA", "black", "#C288BB", "#7FBF74", "#878787",
                                                                          "black", "#C05558", "black", "#71519C", "#C288BB", "#7FBF74", "#878787")))))) + 
  wrap_elements(full = (ggdraw() + 
                          draw_image(image_read_svg(here("raw_data_analysis/figures/terminator_construct_designs/tRPS3-tTSA1_design.svg"), width = 520, height = 354))))+
    plot_layout(design = insertion_layout))) /
  RPS3_deltadeltacq_plot /
  TSA1_deltadeltacq_plot +
  plot_annotation(tag_levels = "A") +
  plot_layout(heights = c(4,1,1))

ggsave(filename = here("results_chapter/figures/insertion_constructs_design_and_qpcr.png"),
       width = 163,
       height = 200,
       unit = "mm",
       dpi = 300)

### PIR1 deletion diagram

construct_to_label_dictionary_PIR1 = tibble(construct = c("modG", "modF", "modE", "modD", "modC", "modA", "modB", "WT"), 
                                                 label = c("mod_NTNNN", "mod_ANNNN", "mod_ATNNN", "mod_ATHNH", "mod_ATNHH", "mod_ANHHH", "mod_NTHHH", "WT"),
                                                 PIR1_UTR3 = str_c("tPIR1_",construct))

PIR1_deltadeltacq <- read_csv(here("raw_data_analysis/data/norm_qpcr/motif_context_dependence/pPIR1_pPGK1_pSRO9_tPIR1_deltadeltacq_platesnorm_summarise.csv")) %>%
  mutate(UTR3 = factor(UTR3,levels = construct_to_label_dictionary_PIR1$PIR1_UTR3),
         promoter = factor(promoter, levels = c("pPIR1","pPGK1","pSRO9")))

PIR1_deltadeltacq_plot <- ggplot(data = PIR1_deltadeltacq)+
  RNA_relative_abundance_figure_options +
  scale_colour_manual(values=PIR1_colour_scheme) +
  theme(axis.text.y=element_text(colour=PIR1_colour_scheme)) +
  labs(x ="Fold change in RNA abundance relative to tPIR1 mod_WT (log2 scale)", y = "") +
  facet_wrap(~promoter,ncol = 3) +
  scale_y_discrete(breaks = construct_to_label_dictionary_PIR1$PIR1_UTR3, labels = construct_to_label_dictionary_PIR1$label)

deletion_layout <- "AC
                    BC"

## Can we modularize this plot to produce the tableGrob in one place, then call it from the draw?
deletion_construct_plot <- wrap_elements(full = (wrap_elements(full = (ggdraw() + 
                                                                          draw_image(image_read_svg(here("raw_data_analysis/figures/terminator_construct_designs/PIR1_motif_WT_construct_design.svg"))))) +
                                                    wrap_elements(full = tableGrob(tibble("Construct" = rev(c("mod_NTNNN", "mod_ANNNN", "mod_ATNNN", "mod_ATHNH", 
                                                                                                          "mod_ATNHH", "mod_ANHHH", "mod_NTHHH", "WT")), 
                                                                                          "Motif 1" = c("ATATTC", "Random", "ATATTC", 
                                                                                                        "ATATTC", "ATATTC", "ATATTC", "ATATTC", "Random"), 
                                                                                          "Motif 2" = c("TGTAHMNTA", "TGTAHMNTA", "Random", 
                                                                                                        "TGTAHMNTA", "TGTAHMNTA", "TGTAHMNTA", "Random", "TGTAHMNTA"), 
                                                                                          "Motif 3" = c("HWNCATTWY", "HWNCATTWY", "HWNCATTWY", 
                                                                                                          "Random", "HWNCATTWY", "Random", "Random", "Random"),
                                                                                          "Motif 4" = c("HWNCATTWY", "HWNCATTWY", "HWNCATTWY", 
                                                                                                          "HWNCATTWY", "Random", "Random", "Random", "Random"),
                                                                                          "Motif 5" = c("HWNCATTWY", "HWNCATTWY", "HWNCATTWY", 
                                                                                                          "HWNCATTWY", "HWNCATTWY", "Random", "Random", "Random")),
                                                                                   rows = NULL,
                                                                                   theme = ttheme_default(base_size=5,
                                                                                                          padding = unit(c(2, 2), "mm"),
                                                                                                          core=list(
                                                                                                            fg_params=list(col = c("#878787", "#C15659", "#5978BA", "#4EB0B5", "#D97F1D", "#71519C", "#968BC2", "#DD76A5",
                                                                                                                                   "black", "#C15659", "black", "black", "black", "black", "black", "#DD76A5",
                                                                                                                                   "black", "black", "#5978BA", "black", "black", "black", "#968BC2", "black",
                                                                                                                                   "black", "black", "black", "#4EB0B5", "black", "#71519C", "#968BC2", "#DD76A5",
                                                                                                                                   "black", "black", "black", "black", "#D97F1D", "#71519C", "#968BC2", "#DD76A5",
                                                                                                                                   "black", "black", "black", "black", "black", "#71519C", "#968BC2", "#DD76A5")))))) + 
                                                    wrap_elements(full = (ggdraw() + 
                                                                            draw_image(image_read_svg(here("raw_data_analysis/figures/terminator_construct_designs/tPIR1_design.svg"), width = 520, height = 354))))+
                                                    plot_layout(design = deletion_layout))) /
  PIR1_deltadeltacq_plot +
  plot_annotation(tag_levels = "A") +
  plot_layout(heights = c(3,1))

ggsave(here("results_chapter/figures/tPIR1_design_and_qpcr.png"),
       deletion_construct_plot,
       width = 163,
       height = 150,
       units = "mm",
       dpi = 300)
