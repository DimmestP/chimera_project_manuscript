# count_codons_all_ORFs.R
# counts the codons on every yeast open reading frame
# warning: this takes a few minutes to run
library(here)
library(tidyr)
library(purrr)
library(Biostrings)
source(here("raw_data_analysis/code/linear_model_functions.R"))

# import yeast open reading frame dataset
Scer_ORF <- readDNAStringSet("https://downloads.yeastgenome.org/sequence/S288C_reference/orf_dna/orf_coding_all.fasta.gz")
Scer_ORF_name <- as_tibble(names(Scer_ORF)) %>%
  separate(value,c("geneName",NA),extra="drop",sep=" ")

# Count the number of each codon in each yeast ORF
# and frequency by dividing by gene length
codon_count_freq_long <- 
    tibble(geneName = Scer_ORF_name$geneName,
           ORF = as.character(Scer_ORF),
           length = width(Scer_ORF)) %>%
  filter((length %% 3) == 0) %>%
  group_by(geneName) %>%
  do(count_codons(.$ORF)) %>%
  mutate(freq = round(counts / sum(counts), digits = 5) )

write_tsv(codon_count_freq_long, 
          here("raw_data_analysis/data/half_life_data/scer_codon_counts_frequency_long.txt") )


# Convert to relative frequncy of each codon by dividing by gene length 
codon_freq_wide <- codon_count_freq_long %>%
  select(-counts) %>%
  pivot_wider(names_from = codon, values_from = freq, values_fill = 0)

write_tsv(codon_freq_wide, 
          here("raw_data_analysis/data/half_life_data/scer_codon_frequency_table.txt") )
