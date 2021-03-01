library(magick)
library(patchwork)
library(gridExtra)

source(here("raw_data_analysis/code/shared_figure_formatting.R"))

layout <- "AAA
           BCC"


construct_to_label_dictionary = tibble(construct = c("WT","modC","modE","modD","modA","modB","mod0"), 
                                       label = c("WT", "mod_NGG", "mod_HTH", "mod_HNH", "mod_NTN", "mod_NAA", "mod_NNN"),
                                       TSA1_UTR3 = str_c("tTSA1_",construct),
                                       RPS3_UTR3 = str_c("tRPS3_",construct))

TSA1_deltadeltacq <- read_csv(here("raw_data_analysis/data/norm_qpcr/motif_context_dependence/pTSA1_pPGK1_pSRO9_tTSA1_deltadeltacq_platesnorm_summarise.csv")) %>%
  mutate(UTR3 = factor(UTR3,levels = construct_to_label_dictionary$TSA1_UTR3))

RPS3_deltadeltacq <- read_csv(here("raw_data_analysis/data/norm_qpcr/motif_context_dependence/pRPS3_pPGK1_pSRO9_tRPS3_deltadeltacq_platesnorm_summarise.csv")) %>%
  mutate(UTR3 = factor(UTR3,levels = construct_to_label_dictionary$RPS3_UTR3))

TSA1_deltadeltacq_plot <- ggplot(data = TSA1_deltadeltacq)+
  RNA_relative_abundance_figure_options +
  scale_colour_manual(values=RPS3_TSA1_colour_scheme) +
  theme(axis.text.y=element_text(colour=RPS3_TSA1_colour_scheme)) +
  labs(x ="Fold change in RNA abundance relative to tTSA1-mod0 (log2 scale)", y = "") +
  facet_wrap(~pro_mCh,ncol = 3) +
  scale_y_discrete(breaks = construct_to_label_dictionary$TSA1_UTR3, labels = construct_to_label_dictionary$label)

RPS3_deltadeltacq_plot <- ggplot(data = RPS3_deltadeltacq) +
  RNA_relative_abundance_figure_options +
  scale_colour_manual(values=RPS3_TSA1_colour_scheme) +
  theme(axis.text.y=element_text(colour=RPS3_TSA1_colour_scheme)) +
  labs(x ="Fold change in RNA abundance relative to tRPS3-mod0 (log2 scale)", y = "") +
  facet_wrap(~pro_mCh,ncol = 3) +
  scale_y_discrete(breaks = construct_to_label_dictionary$RPS3_UTR3, labels = construct_to_label_dictionary$label)

TSA1_simplified_plot <- wrap_plots(ggdraw() + 
                                draw_image(magick::image_read_svg(here("raw_data_analysis/figures/terminator_construct_designs/TSA1_motif_mod0_construct_design.svg"))), 
                              gridExtra::tableGrob(tibble("Constuct" = c("mod_N", "mod_A", "mod_T",
                                                                         "mod_H", "mod_HTH", "mod_G", "WT"), 
                                                          "Motif 1" = c("Random", "Random", "Random", 
                                                                        "HWNCATTWY", "HWNCATTWY", "Random", "NA"), 
                                                          "Motif 2" = c("Random", "ATATTC", "TGTAHMNT", 
                                                                        "Random", "TGTAHMNT", "GTATACCTA", "NA"), 
                                                          "Motif 3" = c("Random", "ATATTC", "Random", 
                                                                        "HWNCATTWY", "HWNCATTWY", "GTATACCTA", "NA")), 
                                                   rows = NULL), 
                              TSA1_deltadeltacq_plot) +
  plot_layout(design = layout)

RPS3_simplified_plot <- wrap_plots(ggdraw() + 
                                draw_image(magick::image_read_svg(here("raw_data_analysis/figures/terminator_construct_designs/RPS3_motif_mod0_construct_design.svg"))), 
                              gridExtra::tableGrob(tibble("Motif 1" = c("Random", "Random", "Random", 
                                                                        "HWNCAUUWY", "HWNCAUUWY", "Random"), 
                                                          "Motif 2" = c("Random", "ATATTC", "UGUAHMNU", 
                                                                        "Random", "UGUAHMNU", "GTATACCTA"), 
                                                          "Motif 3" = c("Random", "ATATTC", "Random", 
                                                                        "HWNCAUUWY", "HWNCAUUWY", "GTATACCTA")), 
                                                   rows = NULL), 
                              RPS3_deltadeltacq_plot) +
  plot_layout(design = layout)

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
                          theme = ttheme_default(base_size=7,
                                                 padding = unit(c(2, 2), "mm"),
                                                 core=list(
                                                   fg_params=list(col = rev(RPS3_TSA1_colour_scheme)))))) + 
  wrap_elements(full = (ggdraw() + 
                          draw_image(here("raw_data_analysis/figures/terminator_construct_designs/design_diagram.png"))))+
    plot_layout(design = insertion_layout))) /
  TSA1_deltadeltacq_plot /
  RPS3_deltadeltacq_plot +
  plot_annotation(tag_levels = "A") +
  plot_layout(heights = c(4,1,1))

ggsave(filename = here("results_chapter/figures/insertion_constructs_design_and_qpcr.png"), width = 183, height = 200, unit = "mm")
