# import and prepare half-life data

chan_decay_raw <-read_tsv("./data/new_chan_dr_data.txt")
chan_decay_hlife <- chan_decay_raw %>%
  dplyr::rename(transcriptName = gene_id, geneAltId = gene_short_name, hlifeR1 = halflife_160412_r1, hlifeR2 = halflife_160412_r2) %>%
  select(-geneAltId) %>%
  transmute(transcriptName,hlife = rowMeans(cbind(hlifeR1, hlifeR2), na.rm = TRUE)) %>%
  filter(is.finite(hlife))

# import list of collated 3'UTR motifs
motifs_raw <- scan("./data/list_motifs.txt", character())

# import yeast open reading frame dataset
Scer_ORF <- readDNAStringSet("https://downloads.yeastgenome.org/sequence/S288C_reference/orf_dna/orf_coding_all.fasta.gz")
Scer_ORF_name <- as_tibble(names(Scer_ORF)) %>%
  separate(value,c("transcriptName",NA),extra="drop",sep=" ")

# Gather all yeast 3'UTRs
yeast_3UTRs <- read_csv("./data/whole_genome_3UTR.csv")

# import list of collated 3'UTR motifs
motifs_raw <- scan("./data/list_motifs.txt", character())

# import list of codons
codon <- readRDS("./data/codons.rds")
codon_no_TTT <- codon[-1]

# Find the frequency of codons for each yeast ORF
codon_freq <- tibble(geneName = Scer_ORF_name$transcriptName,ORF=as.character(Scer_ORF),length=width(Scer_ORF)) %>% 
  mutate(ORF = gsub("([ATCG]{3})([ATCG]{3})",'\\1,\\2,',as.character(ORF))) %>%
  filter((length %% 3) == 0) %>%
  separate_rows(ORF,sep = ",") %>% 
  group_by(geneName,ORF) %>%
  summarise(counts=n()) %>%
  spread(key = ORF,value = counts,fill=0) %>%
  select(-V1) %>%
  inner_join(tibble(geneName = Scer_ORF_name$transcriptName, geneLength = width(Scer_ORF)),by="geneName") %>%
  gather(key=codon,value=number,-geneLength,-geneName) %>%
  mutate(number=number/geneLength) %>% 
  spread(key=codon,value=number) %>%
  select(-TTT, -geneLength) %>%
  dplyr::rename(transcriptName = geneName )

# IUPAC to regex function
iupac_to_regex <- function(iupac_string){
  iupac_string %>% str_replace_all(c("U" = "T", "W" = "(A|T)", "S" = "(C|G)", "M" = "(A|C)", "K" = "(G|T)", "R" = "(A|G)", "Y" = "(C|T)", "B" = "(C|G|T)", "D" = "(A|G|T)", "H" = "(A|C|T)", "V" = "(A|C|G)", "N" = "(A|C|G|T)"))
}

# function to create all alternative version of motifs explicitally as strings
regex_to_string <- function(motif_regex){
  all_regex_as_string <- motif_regex
  while(sum(str_detect(all_regex_as_string,"\\("))>0){
    current_regex_as_string <- vector(mode = "character", length = 0)
    for(i in 1:length(all_regex_as_string)){
      current_regex <- all_regex_as_string[i]
      extracted_alt_nuc <- substring(current_regex,str_locate(current_regex,"\\(")[1,1],str_locate(current_regex,"\\)")[1,1])
      alternative_SNP <- str_remove(str_split(extracted_alt_nuc,"\\|")[[1]],"\\(|\\)")
      extracted_alt_nuc <- str_replace(extracted_alt_nuc,"\\(","\\\\\\(") %>%
        str_replace("\\)","\\\\\\)") %>%
        str_replace_all("\\|","\\\\\\|")
      temp_regex_as_string <- vector(mode="character")
      for(j in 1:length(alternative_SNP)){
        temp_regex_as_string[j] <- str_replace(current_regex,extracted_alt_nuc,alternative_SNP[j])
      }
      current_regex_as_string <- c(current_regex_as_string,temp_regex_as_string)
    }
    all_regex_as_string <- current_regex_as_string
  }
  all_regex_as_string
}

# Dictionary for converted IUPAC codes into machine readable regular expressions and converting U -> T
motifs_dictionary <- iupac_to_regex(motifs_raw) 

# Manipulations to remove double counting of motifs
motifs_all_alternatives <- tibble(motifIUPAC = motifs_raw,motifsRegex = motifs_dictionary) %>% 
  group_by(motifIUPAC) %>%
  mutate(motifsStrings = map(motifsRegex,regex_to_string)) %>%
  unnest(motifsStrings) 

motifs_unique_alternatives <- motifs_all_alternatives %>% 
  group_by(motifsStrings) %>%
  mutate(pairedMotifString = list(pairedMotifString=motifs_all_alternatives$motifsStrings)) %>%
  unnest(pairedMotifString) %>%
  group_by(motifsRegex) %>%
  filter(!(motifsStrings==pairedMotifString)) %>%
  ungroup() %>%
  mutate(matchingMotif = str_detect(motifsStrings,pairedMotifString)) %>%
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
unique_IUPAC <- motifs_unique_alternatives %>% distinct(newMotifsRegex,.keep_all = TRUE) %>% pull(newMotifIUPAC)

#Search and add frequency of each c(motif) as a column in ref dataset
for (i in 1:length(unique_IUPAC)){
  motif_count <- str_count(single_count_median_3UTR_motifs_freq$threePrimeUTR, str_c(motifs_unique_alternatives %>% filter(newMotifIUPAC == unique_IUPAC[i]) %>% pull(motifsStrings),collapse = "|"))
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
