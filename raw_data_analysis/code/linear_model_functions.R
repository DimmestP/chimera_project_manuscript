
library(dplyr)
library(tibble)
library(stringr)
library(readr)
library(testthat)

# list of sense codons
sense_codons <- c("TTT", "TTC", "TTA", "TTG", "TCT", "TCC", "TCA", "TCG", "TAT", "TAC", "TGT", "TGC", "TGG", "CTT", "CTC", "CTA", "CTG", "CCT", "CCC", "CCA", "CCG", "CAT", "CAC", "CAA", "CAG", "CGT", "CGC", "CGA", "CGG", "ATT", "ATC", "ATA", "ATG", "ACT", "ACC", "ACA", "ACG", "AAT", "AAC", "AAA", "AAG", "AGT", "AGC", "AGA", "AGG", "GTT", "GTC", "GTA", "GTG", "GCT", "GCC", "GCA", "GCG", "GAT", "GAC", "GAA", "GAG", "GGT", "GGC", "GGA", "GGG")

# remove one codon to remove colinearity issue in linear model
sense_codons_no_TTT <- sense_codons[-1]

# IUPAC to regex function
iupac_to_regex <- function(iupac_string){
  iupac_string %>% str_replace_all(c("U" = "T", "W" = "(A|T)", "S" = "(C|G)", "M" = "(A|C)", "K" = "(G|T)", "R" = "(A|G)", "Y" = "(C|T)", "B" = "(C|G|T)", "D" = "(A|G|T)", "H" = "(A|C|T)", "V" = "(A|C|G)", "N" = "(A|C|G|T)"))
}
#test_that("IUPAC values are converted to their redundent nucleotides.", {
#  expect_equal(iupac_to_regex("USRN"), "T(C|G)(A|G)(A|C|G|T)")
#})

# function to create all alternative version of motifs explicitly as strings
expand_regex_to_redundent_strings <- function(motif_regex){
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

#test_that("redundent regex expression are expended out into all possible string combinations.", {
#  expect_equal(expand_regex_to_redundent_strings("T(C|G)(A|G)(A|C|G|T)"), c("TCAA","TCAC","TCAG","TCAT","TCGA","TCGC","TCGG","TCGT","TGAA","TGAC","TGAG","TGAT","TGGA","TGGC","TGGG","TGGT"))
#})

# function to count motifs in a vector of sequences
# regex_motifs is a vector of motifs to search for in regex format
# sequence is a vector of DNA or RNA sequences to search through
# min_count_filter is minimal number of times a motif needs to be detected to be included
# alt_motif_name is an optional argument that can be used to change column names associated with the regex motifs
# gene_name is an optional argument that adds the gene name associated with a sequence
# returns a tibble consisting of a column for each given motif containing its counts in the respective sequence
motif_count_function <- function(regex_motifs, sequence, min_count_filter = 5, gene_name = NA){
  motif_freq <- tibble(Sequence = sequence)
  for (i in 1:length(regex_motifs)){
    motif_count <- str_count(sequence, regex_motifs[i])
    if(sum(motif_count) > min_count_filter)  { 
      motif_freq <- mutate(motif_freq, !!regex_motifs[i] := motif_count)
    }
  }
  if(length(gene_name)>1) {
    motif_freq <- motif_freq %>% 
      mutate(geneName = gene_name) %>% 
      relocate(geneName)}
  else{
    if(!is.na(gene_name)) {
      # what does this do?
      motif_freq <- motif_freq %>% 
      mutate(geneName = gene_name) %>% 
      relocate(geneName)
    }
  }
  return(motif_freq)
}
#test_that("motifs are counted correctly and the tibble output is as expected.", {
#  expect_equal(motif_count_function(c("TTT","TTC"),"TTTTTCTTT",0,c("Motif_1","Motif_2"),"My Example"), tibble(geneName = "My Example", transcriptSequence = "TTTTTCTTT", Motif_1 = 2, Motif_2 = 1))
#})

# Function to find motifs associated with decay for all isoforms 
# Warning: hard-coded formulae in this make the functions non-portable
greedy_linear_motif_selection <- function(motif_counts_half_life_data,observable){
  motif_counts_half_life_data_glob <<- motif_counts_half_life_data
  
  # predict half life with codon usage only
  base_codon_model <- lm(paste0("log2(",observable,")~",str_flatten(sense_codons_no_TTT,"+"),"+","UTR3_length"),data = motif_counts_half_life_data_glob)
  
  # select motifs to use on top of codons to predrict half life using a greedy AIC method 
  raw_motifs_list <- colnames(single_count_median_3UTR_motifs_freq)[3:length(single_count_median_3UTR_motifs_freq)]
  
  raw_motifs_list <- raw_motifs_list[raw_motifs_list %in% colnames(motif_counts_half_life_data)]
  
  step(base_codon_model,paste0("log2(",observable,")~",str_flatten(sense_codons_no_TTT,"+"),"+",str_flatten(raw_motifs_list,"+"),"+","UTR3_length"),trace=FALSE)
  
}
