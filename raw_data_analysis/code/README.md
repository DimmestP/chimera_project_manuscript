## Folder Contents

### raw_qpcr_analysis
Contains the tidyqpcr analysis files for all of the qpcr data presented in raw_data_analysis/data/norm_qpcr.

### raw_platereader_analysis
Contains the omniplate analysis files for all of the platereader data presented in raw_data_analysis/data/platereader_qpcr.

### train_half_life_linear_model.R
Predicts the half lives of all transcripts in the yeast genome using codon usage and motif presence. It requires the use of functions contained in linear_model_functions.R. 

### qpcr_linear_model.R
Predicts the contributions of motifs to the abundance of contructs. Creates the figure in  results_chapter/figures/qPCR_model_coef_and_pred_vs_exp_abund.png

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
Uses results from train_half_life_linear_model.R and qpcr_linear_model.R to make figure results_chapter/figures/hlife_motif_pred_vs_qpcr_abund.Rmd

### pro_ter_swaps_platereader_and_qPCR_plot.R
Uses results from pro-ter swap qPCR and platereader analysis to create results_chapter/figures/pro_ter_swaps_platereader_and_qPCR_plot.R and results_chapter/figures/pro_ter_platereader_norm_mTurq_and_mCh.png


