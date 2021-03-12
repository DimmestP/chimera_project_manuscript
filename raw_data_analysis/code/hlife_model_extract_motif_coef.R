library(here)
library(tidyqpcr)
library(latex2exp)
library(dplyr)
library(tibble)
library(stringr)
library(readr)
library(patchwork)
# requires train_hlife_life_linear_model.R to be ran first
#
#
source(here("raw_data_analysis/code/shared_figure_formatting.R"))
load(file = here("raw_data_analysis/data/hlife_model_summary"))

shortlisted_motifs <- chan_motif_coefficients %>%
  filter(term %in% c("ATATTC", "HWNCATTWY", "GTATACCTA", "TGTAHMNTA") | term %in% (sun_motif_coefficients %>% pull(term))) %>%
  pull(term)

# run linear models using only shortlisted motifs (Chan data set)
chan_shortlisted_motif_model <- lm(paste0("log2(hlife)~", str_flatten(codon_no_TTT, "+"), "+ UTR3_length +", str_flatten(shortlisted_motifs, "+")), data = single_count_decay_prediction_dataset_chan)

# run linear models using only shortlisted motifs (Sun data set)
sun_shortlisted_motif_model <- lm(paste0("log2(hlife)~", str_flatten(codon_no_TTT, "+"), "+ UTR3_length +", str_flatten(shortlisted_motifs, "+")), data = single_count_decay_prediction_dataset_sun)

combined_motif_coefficients <- broom::tidy(chan_shortlisted_motif_model) %>% 
  filter(!(term %in% codon_no_TTT), term != "(Intercept)",term != "UTR3_length") %>%
  inner_join(broom::tidy(sun_shortlisted_motif_model) %>% 
              filter(!(term %in% codon_no_TTT), term != "(Intercept)",term != "UTR3_length"),
             by = "term",
             suffix = c("_C", "_S"))

# output motif coefficients graph
model_coefficients <-   ggplot(combined_motif_coefficients, aes(y = estimate_C, x = estimate_S, ymin=estimate_C - std.error_C ,ymax=estimate_C + std.error_C, xmin=estimate_S - std.error_S, xmax=estimate_S + std.error_S, colour=term)) +
  geom_hline(yintercept = 0,size = 0.2) +
  geom_errorbar() +
  geom_errorbarh() +
  theme(axis.text.x=element_text(angle=90,vjust = 0.5),
        panel.grid.minor=element_blank(),
        legend.position = "bottom") +
  labs(x=TeX("$\\Delta \\log_2$ $\\lambda^{1/2}_{Sun}$"),y=TeX("$\\Delta \\log_2$ $\\lambda^{1/2}_{Chan}$"), colour = "Motif") +
  scale_shape_manual(values=1:7) +
  coord_equal() +
  scale_colour_discrete(guide = guide_legend(nrow = 3))

# Create full hlife model summary figure

ggsave2(here("results_chapter/figures/hlife_model_multi_fig.png"), 
        (dataset_comparison + model_coefficients) / 
          (chan_pred_vs_obvs_plot + (sun_pred_vs_obvs_plot + plot_layout(tag_level = "new"))) + 
          plot_annotation(tag_levels = "A"), 
        width = 163, 
        height = 180,
        units = "mm")

# output list of chan motif coefficients with error
write_csv( chan_motif_coefficients %>% select(term, estimate, std.error) %>% rename( "Motif"= "term", "Coefficient"= "estimate"), here("./raw_data_analysis/data/chan_motif_coefficients_with_error.csv"))
