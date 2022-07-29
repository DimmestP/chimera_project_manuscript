## Folder Contents

### raw_qpcr_analysis
Contains the tidyqpcr analysis files for all of the qpcr data presented in raw_data_analysis/data/norm_qpcr.

### raw_platereader_analysis
Contains the omniplate analysis files for all of the platereader data presented in raw_data_analysis/data/platereader_qpcr.

### train_half_life_linear_model.Rmd
Predicts the half lives of all transcripts in the yeast genome using codon usage and motif presence.  Using a greedy algorithm that maximises the AIC to select for the optimum combination of motifs. It requires the use of functions contained in linear_model_functions.R. 

### qpcr_linear_model.R
Predicts the contributions of motifs to the abundance of constructs. Creates the figure in  results_chapter/figures/qPCR_model_coef_and_pred_vs_exp_abund.png

### shared_figure_formatting.R
Contains several default settings to ensure all figures have similar sizes/fonts/orientations.

### linear_model_functions.R
Collection of linear model fitting functions used in train_half_life_linear_model.R

### hlife_model_extract_motif_coef.R
Uses the output of train_half_life_linear_model.R to create the coefficient table in results_chapter/data/chan_motif_coefficients.csv and the summary figure in results_chapter/figures/hlife_model_multi_fig.png

### combine_terminator_construct_qPCR_and_design_figure.R
Uses analysed motif insertion qPCR to make figure results_chapter/figures/insertion_constructs_design_and_qpcr.png and results_chapter/figures/tPIR1_design_and_qpcr.png

Requires the construct design diagram available in;

- raw_data_analysis/figures/terminator_construct_designs/terminator_construct_designs/RPS3_TSA1_motif_mod0_construct_design.svg

- raw_data_analysis/figures/terminator_construct_designs/terminator_construct_designs/PIR1_motif_WT_construct_design.svg

### hlife_motif_pred_vs_qpcr_abund.Rmd
Summarises the results of the qPCR experiments compared to the predicted motif effects from the half life model. Uses results from train_half_life_linear_model.R and qpcr_linear_model.R to make figure results_chapter/figures/hlife_motif_pred_vs_qpcr_abund.Rmd.

### pro_ter_swaps_platereader_and_qPCR_plot.R
Uses results from pro-ter swap qPCR and platereader analysis to create results_chapter/figures/pro_ter_swaps_platereader_and_qPCR_plot.R and results_chapter/figures/pro_ter_platereader_norm_mTurq_and_mCh.png

### count_codons_all_ORFs.R
Counts the codons on every yeast open reading frame and calculates the relative frequency (by dividing by gene length). This is required in the half life model.

### create_motif_position_histograms.R
Produces histograms of the relative position of motifs of interest in their native 3'UTRs. Creates the motif_position_histograms.png figure in supplementary data.

### fps_summary.R
Counts the positions of the 5PSeq reads that are towards the 5'end of a gene (the second read of each pair which is not anchored to polyA). Create the combined_5p_end_read_plot.png present in the supplementary data of the figure.

### polyA_site_usage_plot.Rmd
Create relative polyA site usage plots from the 5PSeq and QuantSeq data sets. Creates the results figure polya_usage_plot_QuantSeq.png and the sup data figures: relative_polya_counts_QuantSeq_and_5PSeq.png, polya_QC_plot_QuantSeq_and_5PSeq.png, plasmid_vs_genome_polyA_combined_plot.png  polya_usage_plot_5PSeq.png

### RNASeq_QC.Rmd
Conducts basic quality control of RNA-seq data. Creates supp data figures sequencing_sample_correlation.png, pRPS3_QuantSeq_5PSeq_comparison.png.



