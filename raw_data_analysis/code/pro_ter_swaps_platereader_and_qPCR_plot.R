library(dplyr)
library(tibble)
library(readr)
library(here)
library(latex2exp)

source(here("raw_data_analysis/code/shared_figure_formatting.R"))

# define promoter and terminator levels in desired order
pro6  <- c("pPGK1", "pHSP26", "pRPS3", "pRPS13", "pSRO9", "pCLN2")
pro4  <- c("pPGK1", "pHSP26", "pRPS3", "pRPS13")
ter10 <- rev(c("tPGK1", "tRPS3", "tRPS13",
           "tPAB1", "tHSP26", 
           "tCLN2", "tSRO9", 
           "tTOS6", "tSUN4", "tPMA1"))

# Create unnormalised platereader plot

pro_mTurq_ter_protein <- read_csv(here("raw_data_analysis/data/norm_platereader/promoter_terminator_swaps/mTurq_collection/pro-mTurq-ter_swaps_raw.csv"))%>%
  mutate(Promoter = factor(Promoter, levels = pro6),
         Terminator = factor(Terminator, levels = ter10)) %>% 
  rename(fluo_per_OD_at_max_gr = "mTurq_per_OD_at_max_gr", Protein = "mTurq") %>%
  select(-X1)

pro_mCh_ter_protein <- read_csv(here("raw_data_analysis/data/norm_platereader/promoter_terminator_swaps/mCherry_collection/pro-mCh-ter_swaps_raw.csv"))%>%
  mutate(Promoter = factor(Promoter, levels = pro6),
         Terminator = factor(Terminator, levels = ter10)) %>% 
  rename(fluo_per_OD_at_max_gr = "mCh_per_OD_at_max_gr", Protein = "mCherry") %>%
  select(-X, -sample_id)

pro_mCh_ter_qPCR <- read_csv(here("raw_data_analysis/data/norm_qpcr/promoter_terminator_swaps/pRPS3_pPGK1_pSRO9_tvariable_deltacq_platesnorm_summarise.csv")) %>%
  filter(target_id == "mCh-7") %>%
  group_by(strain) %>%
  summarise(mRNA = mean(deltacq)) %>%
  rename(Strain = "strain")

low_exp_pro_mCh_mTurq <- pro_mCh_ter_protein %>% 
  bind_rows(pro_mTurq_ter_protein) %>%
  filter(Promoter %in% c("pSRO9", "pCLN2"))

high_exp_pro_mCh_mTurq <- pro_mCh_ter_protein %>% 
  bind_rows(pro_mTurq_ter_protein) %>%
  filter(!Promoter %in% c("pSRO9", "pCLN2"))

mCherry_protein_vs_RNA <- pro_mCh_ter_protein %>%
  group_by(Strain) %>%
  summarise(protein = mean(fluo_per_OD_at_max_gr)) %>%
  inner_join(pro_mCh_ter_qPCR, by = "Strain") %>%
  separate(Strain, into = c(NA,"Terminator"), sep = "(?<=y-)", remove = FALSE) %>%
  separate(Strain, into = c("Promoter",NA), sep = "(?=-m)", remove = FALSE) %>%
  mutate(Terminator = factor(Terminator, levels = ter10))

mCherry_protein_vs_RNA_cor <- cor(mCherry_protein_vs_RNA$protein,
                                  mCherry_protein_vs_RNA$mRNA)

mCherry_protein_vs_RNA_cor_text <- paste(TeX(paste0("R = ", signif(abs(mCherry_protein_vs_RNA_cor),digits = 3))))

mCherry_protein_vs_RNA_figure <- ggplot(mCherry_protein_vs_RNA) +
  geom_point(aes(y = protein, x = -mRNA, colour = Terminator, shape = Promoter)) +
  labs(y="mCherry Fluorescence per \n OD at max growth rate", x = TeX("mRNA abundance ($\\Delta$Cq)")) +
  scale_colour_hue(h = c(0, 360)+20,l=60,c=60) +
  guides(colour="none") +
  scale_y_log2nice(omag = seq(-10, 15)) +
  geom_text(aes(label = mCherry_protein_vs_RNA_cor_text),
            y = 1.4,
            x= -0.8, parse = TRUE, size = 4) +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 10)) +
  guides(shape=guide_legend(ncol=1))

high_exp_mCh_mTurq_platereader_raw_figure <- ggplot(high_exp_pro_mCh_mTurq) +
  protein_raw_abundance_figure_options +
  facet_grid(Protein~Promoter) +
  labs(x="Fluorescence per OD at max growth rate") 

low_exp_mCh_mTurq_platereader_raw_figure <- ggplot(low_exp_pro_mCh_mTurq) +
  protein_raw_abundance_figure_options +
  facet_grid(Protein~Promoter) +
  labs(x="Fluorescence per OD \n at max growth rate") +
  theme(legend.position = "bottom") +
  guides(colour=guide_legend(ncol=2))

bottom_row_swap_figure <- 
    plot_grid(low_exp_mCh_mTurq_platereader_raw_figure,
            mCherry_protein_vs_RNA_figure,
            rel_widths = c(1.3,1),
            ncol = 2)

composite_pro_ter_swap_figure <- 
  plot_grid(high_exp_mCh_mTurq_platereader_raw_figure,
            bottom_row_swap_figure,
            ncol = 1,
            rel_heights = c(1,1.3),
            labels = c("A","B","C"))
composite_pro_ter_swap_figure

ggsave(here("results_chapter/figures/pro_ter_swap_protein_and_rna_exp_figure.png"),
       composite_pro_ter_swap_figure,
       width = 165, height = 200, units = "mm", dpi=300)

#####
# Create Normalised Platereader plot and add qPCR of short-vs-long

shortvslong_platesnorm_all_median <- read_csv(here("raw_data_analysis/data/norm_qpcr/short_vs_long/shortvslong_two_exp_rep_deltadeltacq_platesnorm_summarise.csv")) %>%
  separate(UTR3, into = c(NA, "UTR3"), sep = "_") %>%
  mutate(UTR3 = factor(UTR3, levels = rev(c("59bp","86bp","200bp","543bp"))))

normalised_plot_SRO9 <- ggplot(data = shortvslong_platesnorm_all_median %>% filter(Promoter == "pSRO9")) +
  RNA_relative_abundance_figure_options +
  labs(x="Fold change in RNA abundance \n relative to tSRO9_543bp \n (log2 scale)", title = "SRO9", y = NULL) +
  scale_colour_manual(values=c("#CC6666", "black")) +
  theme(axis.text.y=element_text(colour=c("#CC6666", "black")), plot.title = element_text(face = "plain", hjust = 0.5, size = 10)) +
  geom_vline(xintercept=1, linetype="dotted", color = "grey50")

normalised_plot_RPS3 <- ggplot(data = shortvslong_platesnorm_all_median %>% filter(Promoter == "pRPS3")) +
  RNA_relative_abundance_figure_options + 
  labs(x="Fold change in RNA abundance \n relative to tRPS3_200bp \n (log2 scale)", title = "RPS3", y="Terminator \n Length") +
  scale_colour_manual(values=c("#288a2e","#a84a9a","#6f3ba1")) +
  theme(axis.text.y=element_text(colour=c("#288a2e","#a84a9a", "#6f3ba1")), plot.title = element_text(face = "plain", hjust = 0.5, size = 10)) +
  geom_vline(xintercept=1, linetype="dotted", color = "grey50")

norm_pro_mCh_ter_all <- read_csv(here("raw_data_analysis/data/norm_platereader/promoter_terminator_swaps/mCherry_collection/pro-mCh-ter_swaps_summary_PGK1_norm.csv"))

norm_pro_mTurq_ter_all <- read_csv(here("raw_data_analysis/data/norm_platereader/promoter_terminator_swaps/mTurq_collection/pro-mTurq-ter_swaps_summary_PGK1_norm.csv"))

norm_mTurq_and_mCherry <- norm_pro_mTurq_ter_all %>%
  rename(protein = "mTurq", fluo_per_OD_at_max_gr = "norm_fluo_per_OD_at_max_gr") %>%
  select(-mTurq_per_OD_at_max_gr, -remove, -X1) %>%
  bind_rows(
    norm_pro_mCh_ter_all %>%
      rename(protein = "mCherry", fluo_per_OD_at_max_gr = "norm_fluo_per_OD_at_max_gr") %>%
      select(-BioRep, -mCh_per_OD_at_max_gr, -sample_id, -remove, -X)
  ) %>%
  filter(!Promoter %in% c("pCLN2", "pSRO9")) %>%
  mutate(Terminator = factor(Terminator, ter10), Promoter = factor(Promoter, pro4))

norm_mTurq_and_mCherry_figure <- ggplot(data=norm_mTurq_and_mCherry %>% filter(!Promoter %in% c("pSRO9", "pCLN2"))) +
  labs(x="Relative fluorescence per OD relative to tPGK1", y="Terminator") +
  facet_grid(protein~Promoter) +
  protein_relative_abundance_figure_options

UTR_length_panelB <- plot_grid(normalised_plot_RPS3, normalised_plot_SRO9, ncol = 2)

composite_norm_UTRlength_figure <- plot_grid(norm_mTurq_and_mCherry_figure, 
                              UTR_length_panelB, 
                              labels = c("A","B"),
                              rel_heights = c(2,1),
                              ncol = 1)
composite_norm_UTRlength_figure

ggsave(here("./results_chapter/figures/pro_ter_platereader_norm_mTurq_and_mCh.png"), 
       composite_norm_UTRlength_figure,
       width = 165, height = 150, units = "mm", dpi = 300)

