# Takes several minutes to run


library(dplyr)
library(purrr)
library(tibble)
library(stringr)
library(readr)
library(BiocManager)
library(Biostrings)
library(DECIPHER)
library(latex2exp)
library(here)
library(cowplot)
library(patchwork)
library(tidyqpcr)

source(here("raw_data_analysis/code/linear_model_functions.R"))
source(here("raw_data_analysis/code/shared_figure_formatting.R"))

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
print("Counting occurences of each motif in each 3'UTR")
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
print("greedy linear motif selection for Chan et al. This step takes a few minutes to run.")
single_motif_chan_decay_step_model_chan <- greedy_linear_motif_selection(single_count_decay_prediction_dataset_chan,"hlife")

print("greedy linear motif selection for Sun et al. This step takes a few minutes to run.")
single_motif_chan_decay_step_model_sun <- greedy_linear_motif_selection(single_count_decay_prediction_dataset_sun,"hlife")

print("Finished motif selection")

# find motifs that contribute the most the predicting half life in each data set
sun_motif_coefficients <- broom::tidy(single_motif_chan_decay_step_model_sun) %>% 
  filter(!(term %in% codon_no_TTT), term != "(Intercept)",term != "UTR3_length") %>%
  filter(p.value < 0.05) %>% 
  arrange(desc(abs(estimate)))

chan_motif_coefficients <- broom::tidy(single_motif_chan_decay_step_model_chan) %>% 
  filter(!(term %in% codon_no_TTT), term != "(Intercept)",term != "UTR3_length") %>%
  filter(p.value < 0.05) %>% 
  arrange(desc(abs(estimate)))


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

# plot options for pred_vs_obvs_plot
halflife_breaks <- c(1,10,100)
halflife_minor_breaks <- c(2,3,4,5,20,30,40,50)
halflife_limits <- c(0.9,200)

options_half_life_plots <-
  list(
    geom_bin2d(bins = 70),
    scale_y_log10(breaks = halflife_breaks,
                  minor_breaks = halflife_minor_breaks,
                  limits = halflife_limits),
    scale_x_log10(breaks = halflife_breaks,
                  minor_breaks = halflife_minor_breaks,
                  limits = halflife_limits),
    coord_equal(),
    theme(legend.position = "none",
          panel.grid.minor.y = element_blank(),
          plot.title = element_text(hjust=0.5), 
          plot.margin = margin(5,0,5,20))
  )

chan_pred_vs_obvs_plot <- 
  ggplot(two_data_sets_predictive_power_tibble %>% filter(Data_Set == "Chan"),
         aes(x=Predicted_Half_Life,y=Measured_Half_Life)) +
  options_half_life_plots +
  labs(x=TeX("$\\lambda^{1/2}_{pred}$"),y=TeX("$\\lambda^{1/2}_{obs}$"),title = TeX("Chan")) +
  annotate("text",
           label = TeX(paste0("$R^2=$",signif(chan_step_model_r_squared,2))),
           size=6, x=4, y=130)

sun_pred_vs_obvs_plot <- 
  ggplot(two_data_sets_predictive_power_tibble %>% filter(Data_Set == "Sun"),
         aes(x=Predicted_Half_Life,y=Measured_Half_Life)) +
  options_half_life_plots +
  labs(x=TeX("$\\lambda^{1/2}_{pred}$"),y=TeX("$\\lambda^{1/2}_{obs}$"),title = TeX("Sun")) +
  annotate("text",
           label = TeX(paste0("$R^2=$",signif(sun_step_model_r_squared,2))),
           size=6, x=4, y=130)

# output chan vs sun comparison graph
combined_hlife_data_sets <- inner_join(sun_decay_hlife, 
                                       chan_decay_hlife, 
                                       suffix=c("_S","_C"), 
                                       by="transcriptName") %>%
  filter(!is.na(hlife_S), !is.na(hlife_C))

decay_data_set_cor <- cor(combined_hlife_data_sets$hlife_S, combined_hlife_data_sets$hlife_C)


dataset_comparison <-  
  ggplot(combined_hlife_data_sets, aes(x=hlife_C, y=hlife_S)) +
  options_half_life_plots +
  labs(y = TeX("$\\lambda^^{1/2}_{Sun}$"), 
       x = TeX("$\\lambda^{1/2}_{Chan}$")) +
  annotate("text",
           label = paste0("R = ",signif(decay_data_set_cor,2)),
           size = 6, x = 40, y = 2)


# output model predictive power and dataset comparison figures

save(chan_pred_vs_obvs_plot, sun_pred_vs_obvs_plot, codon_no_TTT, dataset_comparison, sun_motif_coefficients, chan_motif_coefficients, single_count_decay_prediction_dataset_chan, single_count_decay_prediction_dataset_sun, file = here("raw_data_analysis/data/hlife_model_summary"))

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
