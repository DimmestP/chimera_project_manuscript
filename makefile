#### Define variables ####
MANUSCRIPT = _book/chimera_project_manuscript

# In this make file .Rout and .md files are made as temporary files to confirm .R and .Rmd have ran successfully

# Only chapters dependent on R analysis code are listed
CHAPTERS_RMD = results_chapter/results.Rmd supplementary_data_chapter/supplementary_data.Rmd
CHAPTERS_MD = $(CHAPTERS_RMD:.Rmd=.md)

ANALYSIS_FOLDER = raw_data_analysis/code

MOTIF_QPCR_ANALYSIS_FILES_RMD := $(wildcard $(ANALYSIS_FOLDER)/raw_qpcr_analysis/motif_context_dependence/*.Rmd)
MOTIF_QPCR_ANALYSIS_FILES_MD := $(MOTIF_QPCR_ANALYSIS_FILES_RMD:.Rmd=.md)

CHIMERA_PLATEREADER_ANALYSIS_FILES_RMD := $(wildcard $(ANALYSIS_FOLDER)/platereader/promoter_terminator_swaps/*.Rmd)
CHIMERA_PLATEREADER_ANALYSIS_FILES_MD := $(CHIMERA_PLATEREADER_ANALYSIS_FILES_RMD:.Rmd=.md)

LINEAR_MODELLING_RMD_FILES := $(wildcard $(ANALYSIS_FOLDER)/*.Rmd)
LINEAR_MODELLING_R_FILES := $(wildcard $(ANALYSIS_FOLDER)/*.R)

#### Define function to run Rmd files ####
define KNIT_RMD =
	Rscript -e "knitr::knit(\"$<\",output = \"$@\")"
endef

#### Make rules ####

# Check if the manuscript PDF is older than chapter files, if so remake
$(MANUSCRIPT).pdf : $(CHAPTERS_MD)
	Rscript -e "bookdown::render_book(\"abstract.Rmd\")"

# Check if results chapter is older than relevent analysis files, if so remake
results_chapter/results.md : $(MOTIF_QPCR_ANALYSIS_FILES_MD) $(CHIMERA_PLATEREADER_ANALYSIS_FILES_MD)\
$(ANALYSIS_FOLDER)/raw_qpcr_analysis/shortvslong_two_exp_rep.md\
hlife_model_extract_motif_coef.Rout $(ANALYSIS_FOLDER)/hlife_motif_pred_vs_qpcr_abund.md\
$(ANALYSIS_FOLDER)/qpcr_linear_model.md results_chapter/results.Rmd
	$(KNIT_RMD)
	
# Check if supplementary data chapter is older than relevent analysis files, if so remake
supplementary_data_chapter/supplementary_data.md : supplementary_data_chapter/supplementary_data.Rmd $(ANALYSIS_FOLDER)/platereader/motif_context_dependence/RPS3_PIR1_TSA1_protein_vs_RNA_abund.md $(ANALYSIS_FOLDER)/raw_qpcr_analysis/promoter_terminator_swaps/pRPS3_pPGK1_pSRO9_tvariable_three_bio_rep.md
	$(KNIT_RMD)

# Check if analysis files need to be updated
$(ANALYSIS_FOLDER)/raw_qpcr_analysis/promoter_terminator_swaps/pRPS3_pPGK1_pSRO9_tvariable_three_bio_rep.md\
$(ANALYSIS_FOLDER)/raw_qpcr_analysis/shortvslong_two_exp_rep.md\
$(CHIMERA_PLATEREADER_ANALYSIS_FILES_MD)\
$(ANALYSIS_FOLDER)/platereader/motif_context_dependence/RPS3_PIR1_TSA1_protein_vs_RNA_abund.md\
$(MOTIF_QPCR_ANALYSIS_FILES_MD)\
$(ANALYSIS_FOLDER)/qpcr_linear_model.md : %.md : %.Rmd
	$(KNIT_RMD)

$(ANALYSIS_FOLDER)/hlife_motif_pred_vs_qpcr_abund.md : $(ANALYSIS_FOLDER)/hlife_motif_pred_vs_qpcr_abund.Rmd $(ANALYSIS_FOLDER)/train_half_life_linear_model.Rout
	Rscript -e "knitr::knit(\"$(ANALYSIS_FOLDER)/hlife_motif_pred_vs_qpcr_abund.Rmd\",output = \"$@\")"
	
hlife_model_extract_motif_coef.Rout : $(ANALYSIS_FOLDER)/train_half_life_linear_model.Rout $(ANALYSIS_FOLDER)/hlife_model_extract_motif_coef.R
	R CMD BATCH $(ANALYSIS_FOLDER)/hlife_model_extract_motif_coef.R hlife_model_extract_motif_coef.Rout

$(ANALYSIS_FOLDER)/train_half_life_linear_model.Rout :  %.Rout : %.R
	R CMD BATCH $< $@

# Analysis/chapter files should already be there, nothing needs to be done
$(MOTIF_QPCR_ANALYSIS_FILES) :
	
$(LINEAR_MODELLING_RMD_FILES) :

$(LINEAR_MODELLING_R_FILES) :

$CHIMERA_PLATEREADER_ANALYSIS_FILES_RMD:

$(CHAPTERS_RMD) :

$(ANALYSIS_FOLDER)/raw_qpcr_analysis/shortvslong_two_exp_rep.Rmd :

$(ANALYSIS_FOLDER)/platereader/motif_context_dependence/RPS3_PIR1_TSA1_protein_vs_RNA_abund.Rmd :

$(CHAPTERS_RMD) :

(ANALYSIS_FOLDER)/raw_qpcr_analysis/promoter_terminator_swaps/pRPS3_pPGK1_pSRO9_tvariable_three_bio_rep.Rmd:

clean:
	rm -fv $(LINEAR_MODELLING_R_FILES:.R=.Rout)
	rm -fv $(LINEAR_MODELLING_RMD_FILES:.Rmd=.md)
	rm -fv $(MOTIF_QPCR_ANALYSIS_FILES_MD)
	rm -fv $(CHAPTERS_MD)
	rm -fv $(CHIMERA_PLATEREADER_ANALYSIS_FILES_MD)
	rm -fv $(ANALYSIS_FOLDER)/platereader/motif_context_dependence/RPS3_PIR1_TSA1_protein_vs_RNA_abund.md
	rm -fv $(ANALYSIS_FOLDER)/raw_qpcr_analysis/promoter_terminator_swaps/pRPS3_pPGK1_pSRO9_tvariable_three_bio_rep.md
	rm -fv $(ANALYSIS_FOLDER)/raw_qpcr_analysis/shortvslong_two_exp_rep.md
	
.PHONY : clean
