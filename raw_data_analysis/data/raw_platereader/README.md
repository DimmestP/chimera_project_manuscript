## Platereader file name nomenclature

All the raw platereader files held in this folder have the same nomenclature, for example;

- 20200929-pro-mCherry-ter_swaps_n1

- 20201202-pPIR1-tPIR1mod_n1

First is the date the platereader was ran, then the experimental name the results refer to and finally which experimental replicate the results are part of.

In addition, there are three files for each experiment, namely; 

- .xlsx
- contents.xls
- ODcorrection_Glucose_Haploid.txt

.xlsx holds the fluorescence and OD values for all the wells in the platereader across the entire timecourse and is created by the platereader machine

contents.xls holds the plateplan for the experiment with identifying descriptions for the contents of each well. Multiple contents files may be provided for different experimental replicates.

ODcorrection_Glucose_Haploid.txt is a file required by the Omniplate normalisation function to account for different changes in OD levels due to chromosome copy numbers.

## Folder Contents

### promoter_terminator_swaps
Contains two folders, one for mCh and one for mTurq assays, each containing 6 experimental replicates of the 6 promoter and 10 terminators swap assay.

### motif_context_dependence
Contains two folders, one for each motif host promoter-terminator pairing, each containing 2 experimental replicates of the motif context dependence constructs.
