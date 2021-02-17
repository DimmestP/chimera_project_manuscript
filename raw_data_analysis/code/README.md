## Folder Contents

### raw_qpcr_analysis
Contains the tidyqpcr analysis files for all of the qpcr data presented in raw_data_analysis/data/norm_qpcr.

### raw_platereader_analysis
Contains the omniplate analysis files for all of the platereader data presented in raw_data_analysis/data/platereader_qpcr.

### train_half_life_linear_model.R
Predicts the half lives of all transcripts in the yeast genome using codon usage and motif presence. It requires the use of functions contained in linear_model_functions.R. 

### qpcr_linear_model.R
Predicts the contributions of motifs to the abundance of contructs. Creates the figure in  results_chapter/figures/qPCR_model_coef_and_pred_vs_exp_abund.

### shared_figure_formatting.R
Contains several default settings to ensure all figures have similar sizes/fonts/orientations.

### linear_model_functions.R
Collection of linear model fitting functions used in train_half_life_linear_model.R

### hlife_model_extract_motif_coef.R
Uses the output of train_half_life_linear_model.R to create the coefficient table in results_chapter/data/chan_motif_coefficients.csv and the summary figure in results_chapter/figures/hlife_model_multi_fig.png
