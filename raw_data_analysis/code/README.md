raw_qpcr_analysis folder contains the tidyqpcr analysis files for all of the qpcr data presented here

train_half_life_linear_model.R predicts the half lives of all transcripts in the yeast genome using codon usage and motif presence. It requires the use of functions contained in linear_model_functions.R

qpcr_linear_model.R predicts the contributions of motifs to the abundance of contructs.

protein_vs_RNA_abundance.Rmd produces the manuscript supplimentary figure comparing protein flurescence and mRNA abundance for the motif context constructs.

hlife_motif_pred_vs_qpcr_abund.Rmd produce the figure in the manuscript comparing motif contributions to half life and to mRNA abundence.