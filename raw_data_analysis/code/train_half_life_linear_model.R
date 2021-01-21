# Takes several minutes to run

library(tidyverse)
library(Biostrings)
library(DECIPHER)
library(latex2exp)
library(here)
library(cowplot)
library(patchwork)

source(here("raw_data_analysis/code/linear_model_functions.R"))
source(here("raw_data_analysis-code/shared_figure_formatting.R"))

# import chan et al half-life data

chan_decay_raw <-read_tsv(here("./raw_data_analysis/data/half_life_data/chan_decay_data.txt"))

# import sun et al 2013 decay rate data
sun_decay_raw <- read_tsv(here("./raw_data_analysis/data/half_life_data/sun_decay_data.txt"), 
                          locale = locale(decimal = ","))
#change weird name of 3rd column and orf to consistent ones geneName
colnames(sun_decay_raw)[3] <- "d_rate"
colnames(sun_decay_raw)[1] <- "transcriptName" 

#convert decayrate in sun et al to hlife and remove other cols
sun_decay_hlife <- sun_decay_raw %>%
  dplyr::select(transcriptName, d_rate) %>%
  mutate(hlife = log(2)/d_rate) %>% 
  dplyr::select(-d_rate)

# rename chan et al columns and filter out missing data
chan_decay_hlife <- chan_decay_raw %>%
  dplyr::rename(transcriptName = gene_id, geneAltId = gene_short_name, hlifeR1 = halflife_160412_r1, hlifeR2 = halflife_160412_r2) %>%
  transmute(transcriptName,hlife = rowMeans(cbind(hlifeR1, hlifeR2), na.rm = TRUE)) %>%
  filter(is.finite(hlife))

# import yeast open reading frame dataset
Scer_ORF <- readDNAStringSet("https://downloads.yeastgenome.org/sequence/S288C_reference/orf_dna/orf_coding_all.fasta.gz")
Scer_ORF_name <- as_tibble(names(Scer_ORF)) %>%
  separate(value,c("transcriptName",NA),extra="drop",sep=" ")

# Gather all yeast 3'UTRs
yeast_3UTRs <- read_csv(here("./raw_data_analysis/data/half_life_data/whole_genome_3UTR.csv"))

# import list of collated 3'UTR motifs
motifs_raw <- scan(here("./raw_data_analysis/data/half_life_data/collated_suspected_decay_motifs.txt"), character())

# import list of codons
codon <- readRDS(here("./raw_data_analysis/data/half_life_data/codons.rds"))

# remove one codon to remove colinearity issue in linear model
codon_no_TTT <- codon[-1]

# Count the number of each codon in each yeast ORF
codon_count <- tibble(geneName = Scer_ORF_name$transcriptName,ORF=as.character(Scer_ORF),length=width(Scer_ORF)) %>%
  filter((length %% 3) == 0) %>%
  group_by(geneName) %>%
  mutate(data = map(ORF,count_codons)) %>%
  unnest(data) %>%
  pivot_wider(names_from = codon,values_from = counts, values_fill=0)

# Convert to relative frequncy of each codon by dividing by gene length 
codon_freq <- codon_count %>%
  pivot_longer(names_to="codon",values_to="number",cols = c(-length,-geneName,-ORF)) %>%
  mutate(number=number/length) %>% 
  spread(key=codon,value=number) %>%
  select(-TTT, -length) %>%
  dplyr::rename(transcriptName = geneName )



# Dictionary for converted IUPAC codes into machine readable regular expressions and converting U -> T
motifs_dictionary <- tibble(motifIUPAC = motifs_raw,motifsRegex = iupac_to_regex(motifs_raw)) %>% 
  group_by(motifIUPAC) %>%
  mutate(motifsStrings = map(motifsRegex,expand_regex_to_redundent_strings)) %>%
  unnest(motifsStrings) 

# check if any smaller motifs appear in other larger motifs
motifs_duplication_match <- motifs_dictionary %>% 
  group_by(motifsStrings) %>%
  mutate(pairedMotifString = list(pairedMotifString=motifs_dictionary$motifsStrings)) %>%
  unnest(pairedMotifString) %>%
  ungroup() %>%
  filter(!(motifsStrings==pairedMotifString))  %>%
  mutate(matchingMotif = str_detect(motifsStrings,pairedMotifString))
  
# remove any versions of larger motifs that contain other motifs
deduplicated_motifs <- motifs_duplication_match %>% 
  group_by(motifsStrings) %>%
  filter(sum(matchingMotif) == 0) %>%
  select(-matchingMotif,-pairedMotifString) %>%
  distinct(motifsStrings,.keep_all = TRUE) %>%
  group_by(motifIUPAC) %>%
  mutate(newMotifIUPAC =as.character(DECIPHER::ConsensusSequence(DNAStringSet(motifsStrings)))) %>%
  ungroup() %>%
  mutate(newMotifsRegex = iupac_to_regex(newMotifIUPAC)) %>%
  select(-motifIUPAC,-motifsRegex)

# Create 'empty' motif count tibble for median length 3'UTRs 
single_count_median_3UTR_motifs_freq <-  yeast_3UTRs %>%
  filter(!is.na(threePrimeUTR))

# Unduplicated motifs
unique_IUPAC <- deduplicated_motifs %>% distinct(newMotifsRegex,.keep_all = TRUE) %>% pull(newMotifIUPAC)

#Search and add frequency of each c(motif) as a column in ref dataset
for (i in 1:length(unique_IUPAC)){
  motif_count <- str_count(single_count_median_3UTR_motifs_freq$threePrimeUTR, str_c(deduplicated_motifs %>% filter(newMotifIUPAC == unique_IUPAC[i]) %>% pull(motifsStrings),collapse = "|"))
  if(sum(motif_count) > 5)  single_count_median_3UTR_motifs_freq <- mutate(single_count_median_3UTR_motifs_freq %>% ungroup(), !!unique_IUPAC[i] := motif_count)
}

# combine motif,codon and chan decay datasets
single_count_decay_prediction_dataset_chan <- single_count_median_3UTR_motifs_freq %>%
  inner_join(chan_decay_hlife) %>%
  inner_join(codon_freq) %>%
  mutate(UTR3_length = str_length(threePrimeUTR))

# combine motif,codon and sun decay datasets
single_count_decay_prediction_dataset_sun <- single_count_median_3UTR_motifs_freq %>%
  inner_join(sun_decay_hlife) %>%
  inner_join(codon_freq) %>%
  mutate(UTR3_length = str_length(threePrimeUTR))

# Find motifs associated with decay for all isoforms without duplicated motif counts for chan et al
single_motif_chan_decay_step_model_chan <- greedy_linear_motif_selection(single_count_decay_prediction_dataset_chan,"hlife")

single_motif_chan_decay_step_model_sun <- greedy_linear_motif_selection(single_count_decay_prediction_dataset_sun,"hlife")


# find motifs that contribute the most the predicting half life in each data set
sun_motif_coefficients <- broom::tidy(single_motif_chan_decay_step_model_sun) %>% 
  filter(!(term %in% codon_no_TTT), term != "(Intercept)",term != "UTR3_length") %>%
  filter(p.value < 0.05) %>% 
  arrange(desc(abs(estimate)))

chan_motif_coefficients <- broom::tidy(single_motif_chan_decay_step_model_chan) %>% 
  filter(!(term %in% codon_no_TTT), term != "(Intercept)",term != "UTR3_length") %>%
  filter(p.value < 0.05) %>% 
  arrange(desc(abs(estimate)))

combined_motif_coefficients <- chan_motif_coefficients %>%
  inner_join(sun_motif_coefficients, by = "term", suffix = c("_S", "_C")) %>%
  bind_rows(chan_motif_coefficients %>%
              left_join(sun_motif_coefficients, by = "term", suffix = c("_C", "_S")) %>%
              filter(is.na(estimate_S), term %in% c("HWNCATTWY", "ATATTC", "TGTAHMNTA", "GTATACCTA")) %>%
              mutate(estimate_S = 0, std.error_S = 0))

# Compare model performanace on both decay datasets
two_data_sets_predictive_power_tibble <- single_count_decay_prediction_dataset_chan %>%
  select(transcriptName,hlife) %>%
  mutate(Measured_Half_Life = hlife, Predicted_Half_Life = 2^predict(single_motif_chan_decay_step_model_chan),Data_Set = "Chan") %>%
  bind_rows(
    single_count_decay_prediction_dataset_sun %>%
      select(transcriptName,hlife) %>% 
      filter(!is.na(hlife)) %>%
      mutate(Measured_Half_Life =hlife, Predicted_Half_Life = 2^predict(single_motif_chan_decay_step_model_sun),Data_Set = "Sun")
  )

chan_step_model_r_squared <- summary(single_motif_chan_decay_step_model_chan)$r.squared

sun_step_model_r_squared <- summary(single_motif_chan_decay_step_model_sun)$r.squared

chan_pred_vs_obvs_plot <- ggplot(two_data_sets_predictive_power_tibble %>% filter(Data_Set == "Chan")) +
  geom_bin2d(aes(x=Predicted_Half_Life,y=Measured_Half_Life), bins = 70) +
  labs(x=TeX("$\\lambda^{1/2}_{pred}$"),y=TeX("$\\lambda^{1/2}_{obs}$"),title = TeX("Chan")) +
  scale_y_log10(limits = c(0.9,200)) +
  scale_x_log10(breaks=c(1,10),minor_breaks = c(3,30),limits = c(1,60)) +
  theme(legend.position = "none",panel.grid.minor.y = element_blank(),plot.title = element_text(hjust=0.5), plot.margin = margin(5,0,5,20)) +
  annotate("text",label = TeX(paste0("$R^2=$",signif(chan_step_model_r_squared,2))),size=7,x=3,y=130)

sun_pred_vs_obvs_plot <- ggplot(two_data_sets_predictive_power_tibble %>% filter(Data_Set == "Sun")) +
  geom_bin2d(aes(x=Predicted_Half_Life,y=Measured_Half_Life), bins = 70) +
  labs(x=TeX("$\\lambda^{1/2}_{pred}$"),y=TeX("$\\lambda^{1/2}_{obvs}$"),title = TeX("Sun")) +
  scale_y_log10(limits = c(0.9,200)) +
  scale_x_log10(breaks=c(1,10),minor_breaks = c(3,30),limits = c(1,60)) +
  theme(legend.position = "none",panel.grid.minor.y = element_blank(),plot.title = element_text(hjust=0.5),axis.title.y=element_blank(), axis.text.y = element_blank(),axis.ticks.y = element_blank(), plot.margin = margin(5,20,5,0)) +
  annotate("text",label = TeX(paste0("$R^2=$",signif(sun_step_model_r_squared,2))),size=7,x=3,y=130)

# output chan vs sun comparison graph
combined_hlife_data_sets <- inner_join(sun_decay_hlife, 
                                       chan_decay_hlife, 
                                       suffix=c("_S","_C"), 
                                       by="transcriptName") %>%
  filter(!is.na(hlife_S), !is.na(hlife_C))

decay_data_set_cor <- cor(combined_hlife_data_sets$hlife_S, combined_hlife_data_sets$hlife_C)

#ggsave2(here("./results_chapter/figures/chan_vs_sun_plot.png"), 
      dataset_comparison <-  ggplot(combined_hlife_data_sets, aes(x=hlife_C, y=hlife_S)) +
           geom_bin2d( bins=70) + 
           theme(panel.grid.minor = element_blank(), legend.position = "none") +
           labs(y = TeX("$\\lambda^^{1/2}_{Sun}$"), 
                x = TeX("$\\lambda^{1/2}_{Chan}$")) +
          annotate("text",label = paste0("R = ",signif(decay_data_set_cor,2)),size=6,x=40,y=2) +
           scale_y_log10(limits = c(0.9,200)) +
           scale_x_log10(breaks=c(1,10,100),limits = c(0.9,200)) +
           coord_fixed()
        #)

# output model predictive power graph

      model_performance_comparison <- cowplot::plot_grid(chan_pred_vs_obvs_plot,sun_pred_vs_obvs_plot, rel_widths = c(0.58,0.42))
      
#ggsave2(here("./results_chapter/figures/model_predictive_power.png"), cowplot::plot_grid(chan_pred_vs_obvs_plot,sun_pred_vs_obvs_plot, rel_widths = c(0.58,0.42)), height = 5, width = 7)


# output motif coefficients graph
#ggsave2(here("./results_chapter/figures/motif_coefficient_comparison.png"),
   model_coefficients <-   ggplot(combined_motif_coefficients, aes(y = estimate_C, x = estimate_S, ymin=estimate_C - std.error_C ,ymax=estimate_C + std.error_C, xmin=estimate_S - std.error_S, xmax=estimate_S + std.error_S, colour=term)) +
     geom_hline(yintercept = 0,size = 0.2) +
        geom_errorbar() +
        geom_errorbarh() +
        theme(axis.text.x=element_text(angle=90,vjust = 0.5),
              panel.grid.minor=element_blank()) + # , plot.margin = margin(5,40,5,80)) +
        labs(x=TeX("$\\Delta \\log_2$ $\\lambda^{1/2}_{Sun}$"),y=TeX("$\\Delta \\log_2$ $\\lambda^{1/2}_{Chan}$"), colour = "Motif") +
        scale_shape_manual(values=1:7) 
      #)

   ((chan_pred_vs_obvs_plot | sun_pred_vs_obvs_plot) / 
       (dataset_comparison | model_coefficients) ) | 
     gridExtra::tableGrob(chan_motif_coefficients %>% select(term, estimate) %>% rename("term" = "Motif", "estimate" = "Coefficient"), rows = NULL)
   
   ggsave2("results_chapter/figures/hlife_model_multi_fig.png", plot_grid(chan_pred_vs_obvs_plot,
             sun_pred_vs_obvs_plot,dataset_comparison,
            model_coefficients, 
             labels = c("A", "", "B", "C"),
             label_size = 20,
             align = "hv",
             axis = "t",
             scale = c(0.9,0.9,1,0.9)), width = 12, height = 9)
   
# output list of chan motif coefficients 
# write_csv( chan_motif_coefficients %>% select(term, estimate) %>% rename("term" = "Motif", "estimate" = "Coefficient"), here("./results_chapter/data/chan_motif_coefficients.csv"))

# output list of chan motif coefficients with error
#write_csv( chan_motif_coefficients %>% select(term, estimate, std.error) %>% rename("term" = "Motif", "estimate" = "Coefficient"), here("./raw_data_analysis/data/chan_motif_coefficients_with_error.csv"))
   
   
# check for co-occurrences of motifs in native 3'UTR
  TGTAHMNTA_co_occurrences <- single_count_median_3UTR_motifs_freq %>% 
  filter(TGTAHMNTA > 0) %>% 
  select(-threePrimeUTR) %>% 
  gather(key = "motif", value = "count", - transcriptName) %>% 
  group_by(motif) %>% 
  filter(count > 0) %>% 
  count(motif)

HWNCATTWY_co_occurrences <- single_count_median_3UTR_motifs_freq %>% 
  filter(HWNCATTWY > 0) %>% 
  select(-threePrimeUTR) %>% 
  gather(key = "motif", value = "count", - transcriptName) %>% 
  group_by(motif) %>% 
  filter(count > 0) %>% 
  count(motif)
