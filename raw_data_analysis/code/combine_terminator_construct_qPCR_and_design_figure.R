library(magick)
library(patchwork)
library(gridExtra)
library(here)
library(dplyr)
library(readr)

source(here("raw_data_analysis/code/shared_figure_formatting.R"))

### TSA1 and RPS3 insertion diagram

layout <- "AAA
           BCC"


construct_to_label_dictionary_TSA1_RPS3 = tibble(construct = c("WT","modC","modE","modD","modA","modB","mod0"), 
                                       label = c("WT", "mod_NGG", "mod_HTH", "mod_HNH", "mod_NTN", "mod_NAA", "mod_NNN"),
                                       TSA1_UTR3 = str_c("tTSA1_",construct),
                                       RPS3_UTR3 = str_c("tRPS3_",construct))

TSA1_deltadeltacq <- read_csv(here("raw_data_analysis/data/norm_qpcr/motif_context_dependence/pTSA1_pPGK1_pSRO9_tTSA1_deltadeltacq_platesnorm_summarise.csv")) %>%
  mutate(UTR3 = factor(UTR3,levels = construct_to_label_dictionary_TSA1_RPS3$TSA1_UTR3))

RPS3_deltadeltacq <- read_csv(here("raw_data_analysis/data/norm_qpcr/motif_context_dependence/pRPS3_pPGK1_pSRO9_tRPS3_deltadeltacq_platesnorm_summarise.csv")) %>%
  mutate(UTR3 = factor(UTR3,levels = construct_to_label_dictionary_TSA1_RPS3$RPS3_UTR3))

TSA1_deltadeltacq_plot <- ggplot(data = TSA1_deltadeltacq)+
  RNA_relative_abundance_figure_options +
  scale_colour_manual(values=RPS3_TSA1_colour_scheme) +
  theme(axis.text.y=element_text(colour=RPS3_TSA1_colour_scheme)) +
  labs(x ="Fold change in RNA abundance relative to tTSA1-mod0 (log2 scale)", y = "") +
  facet_wrap(~pro_mCh,ncol = 3) +
  scale_y_discrete(breaks = construct_to_label_dictionary_TSA1_RPS3$TSA1_UTR3, labels = construct_to_label_dictionary_TSA1_RPS3$label)

RPS3_deltadeltacq_plot <- ggplot(data = RPS3_deltadeltacq) +
  RNA_relative_abundance_figure_options +
  scale_colour_manual(values=RPS3_TSA1_colour_scheme) +
  theme(axis.text.y=element_text(colour=RPS3_TSA1_colour_scheme)) +
  labs(x ="Fold change in RNA abundance relative to tRPS3-mod0 (log2 scale)", y = "") +
  facet_wrap(~pro_mCh,ncol = 3) +
  scale_y_discrete(breaks = construct_to_label_dictionary_TSA1_RPS3$RPS3_UTR3, labels = construct_to_label_dictionary_TSA1_RPS3$label)

insertion_layout <- "AC
                     AC
                     BC"

insertion_construct_plot <- wrap_elements(full = (wrap_elements(full = (ggdraw() + 
                                                    draw_image(image_read_svg(here("raw_data_analysis/figures/terminator_construct_designs/RPS3_TSA1_motif_mod0_construct_design.svg"))))) +
  wrap_elements(full = tableGrob(tibble("Construct" = c("mod_NNN", "mod_NAA", "mod_NTN",
                                                "mod_HNH", "mod_HTH", "mod_NGG", "WT"), 
                                 "Motif 1" = c("Random", "Random", "Random", 
                                               "HWNCATTWY", "HWNCATTWY", "Random", "NA"), 
                                 "Motif 2" = c("Random", "ATATTC", "TGTAHMNT", 
                                               "Random", "TGTAHMNT", "GTATACCTA", "NA"), 
                                 "Motif 3" = c("Random", "ATATTC", "Random", 
                                               "HWNCATTWY", "HWNCATTWY", "GTATACCTA", "NA")),
                          rows = NULL,
                          theme = ttheme_default(base_size=4,
                                                 padding = unit(c(2, 2), "mm"),
                                                 core=list(
                                                   fg_params=list(col = rev(RPS3_TSA1_colour_scheme)))))) + 
  wrap_elements(full = (ggdraw() + 
                          draw_image(image_read_svg(here("raw_data_analysis/figures/terminator_construct_designs/tRPS3-tTSA1_design.svg"), width = 520, height = 354))))+
    plot_layout(design = insertion_layout))) /
  TSA1_deltadeltacq_plot /
  RPS3_deltadeltacq_plot +
  plot_annotation(tag_levels = "A") +
  plot_layout(heights = c(4,1,1))

ggsave(filename = here("results_chapter/figures/insertion_constructs_design_and_qpcr.png"), width = 163, height = 200, unit = "mm")

### PIR1 deletion diagram

construct_to_label_dictionary_PIR1 = tibble(construct = c("modG", "modF", "modE", "modD", "modC", "modA", "modB", "WT"), 
                                                 label = c("mod_NTNNN", "mod_ANNNN", "mod_ATNNN", "mod_ATHNH", "mod_ATNHH", "mod_ANHHH", "mod_NTHHH", "WT"),
                                                 PIR1_UTR3 = str_c("tPIR1_",construct))

PIR1_deltadeltacq <- read_csv(here("raw_data_analysis/data/norm_qpcr/motif_context_dependence/pPIR1_pPGK1_pSRO9_tPIR1_deltadeltacq_platesnorm_summarise.csv")) %>%
  mutate(UTR3 = factor(UTR3,levels = construct_to_label_dictionary_PIR1$PIR1_UTR3))

PIR1_deltadeltacq_plot <- ggplot(data = PIR1_deltadeltacq)+
  RNA_relative_abundance_figure_options +
  scale_colour_manual(values=PIR1_colour_scheme) +
  theme(axis.text.y=element_text(colour=PIR1_colour_scheme)) +
  labs(x ="Fold change in RNA abundance relative to tPIR1-mod0 (log2 scale)", y = "") +
  facet_wrap(~pro_mCh,ncol = 3) +
  scale_y_discrete(breaks = construct_to_label_dictionary_PIR1$PIR1_UTR3, labels = construct_to_label_dictionary_PIR1$label)

deletion_layout <- "AC
                    BC"

deletion_construct_plot <- wrap_elements(full = (wrap_elements(full = (ggdraw() + 
                                                                          draw_image(image_read_svg(here("raw_data_analysis/figures/terminator_construct_designs/PIR1_motif_WT_construct_design.svg"))))) +
                                                    wrap_elements(full = tableGrob(tibble("Construct" = rev(c("mod_NTNNN", "mod_ANNNN", "mod_ATNNN", "mod_ATHNH", 
                                                                                                          "mod_ATNHH", "mod_ANHHH", "mod_NTHHH", "WT")), 
                                                                                          "ATATTC" = c("", "Random", "", 
                                                                                                        "", "", "", "", "Random"), 
                                                                                          "TGTAHMNTA" = c("", "", "Random", 
                                                                                                        "", "", "", "Random", ""), 
                                                                                          "1st HWNCATTWY" = c("", "", "", 
                                                                                                          "Random", "", "Random", "Random", "Random"),
                                                                                          "2nd HWNCATTWY" = c("", "", "", 
                                                                                                          "", "Random", "Random", "Random", "Random"),
                                                                                          "3rd HWNCATTWY" = c("", "", "", 
                                                                                                          "", "", "Random", "Random", "Random")),
                                                                                   rows = NULL,
                                                                                   theme = ttheme_default(base_size=4,
                                                                                                          padding = unit(c(2, 2), "mm"),
                                                                                                          core=list(
                                                                                                            fg_params=list(col = rev(PIR1_colour_scheme)))))) + 
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
       units = "mm")
