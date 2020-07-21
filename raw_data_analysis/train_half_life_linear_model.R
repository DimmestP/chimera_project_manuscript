library(tidyverse)
library(Biostrings)
library(DECIPHER)

# import chan et al half-life data

chan_decay_raw <-read_tsv("../methods_chapter/data/new_chan_dr_data.txt")

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
yeast_3UTRs <- read_csv("../methods_chapter/data/whole_genome_3UTR.csv")

# import list of collated 3'UTR motifs
motifs_raw <- scan("../methods_chapter/data/list_motifs.txt", character())

# import list of codons
codon <- readRDS("../methods_chapter/data/codons.rds")

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

# combine motif,codon and decay datasets
single_count_decay_prediction_dataset <- single_count_median_3UTR_motifs_freq %>%
  inner_join(chan_decay_hlife) %>%
  inner_join(codon_freq) %>%
  mutate(UTR3_length = str_length(threePrimeUTR))

# Function to find motifs associated with decay for all isoforms 
greedy_linear_motif_selection <- function(motif_counts_half_life_data,observable){
  motif_counts_half_life_data_glob <<- motif_counts_half_life_data
  
  # predict half life with codon usage only
  base_codon_model <- lm(paste0("log2(",observable,")~",str_flatten(codon_no_TTT,"+"),"+","UTR3_length"),data = motif_counts_half_life_data_glob)
  
  # select motifs to use on top of codons to predrict half life using a greedy AIC method 
  raw_motifs_list <- colnames(single_count_median_3UTR_motifs_freq)[3:length(single_count_median_3UTR_motifs_freq)]
  
  raw_motifs_list <- raw_motifs_list[raw_motifs_list %in% colnames(motif_counts_half_life_data)]
  
  step(base_codon_model,paste0("log2(",observable,")~",str_flatten(codon_no_TTT,"+"),"+",str_flatten(raw_motifs_list,"+"),"+","UTR3_length"),trace=FALSE)
  
}

# Find motifs associated with decay for all isoforms without duplicated motif counts for chan et al
single_motif_chan_decay_step_model <- greedy_linear_motif_selection(single_count_decay_prediction_dataset,"hlife")
