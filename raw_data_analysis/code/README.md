raw_qpcr_analysis folder contains the tidyqpcr analysis files for all of the qpcr data presented here. Code for producing Figures 4-6 (The motif qPCR results for the three chosen terminators), and Figure ? (Pro-Term swap qpcr results) can be found here.

train_half_life_linear_model.R predicts the half lives of all transcripts in the yeast genome using codon usage and motif presence. It requires the use of functions contained in linear_model_functions.R. Code making Figure 3 (Results from training the linear model of half life) can be found here.

protein_vs_RNA_abundance.Rmd produces the manuscript supplimentary figure 2 comparing protein fluorescence and mRNA abundance for the motif context constructs.

qpcr_linear_model.R predicts the contributions of motifs to the abundance of contructs.

hlife_motif_pred_vs_qpcr_abund.Rmd produce Figure 7 in the manuscript comparing motif contributions to half life and to mRNA abundence. Requires the files outputted by qpcr_linear_model.R to run.

shared_figure_formatting.R contains several default settings to ensure all figures have similar sizes/fonts/orientations.
