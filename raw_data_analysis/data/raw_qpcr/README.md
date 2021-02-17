## RT-qPCR file name nomenclature

All the raw qPCR files held in this folder have the same nomenclature, for example;

- JA_20200830-shortvslong_n1

- JA_20201212-pPGK1-tPIR1mod-n1

First two letters state the initials of the person who ran the experiment, then the date the qPCR was ran, the experimental name the results refer to and finally which experimental replicate the results are part of.

In addition, there are four files for each experiment, namely; 

- .txt
- .ixo
- melt.txt
- ct.txt

.ixo is a propritary file type created by LightCycler qPCR machines and contains all the entire script and results ran by the machine.

ct.txt holds the final No of cycles each sample needed to reach the threshold value.

.txt holds the full time course values for plotting amplification curves.

melt.txt  holds the full time course values for plotting melt curves.

## Folder Contents

### motif_context_dependence
Contains one folder for each of the 9 promoter-terminator pairings used to test motif context. Within each subfolder there are 8 raw qPCR files (4 qPCR files for 2 experimental replicates)

### promoter_terminator_swaps
Contains qPCR reults of the promoter-terminator swaps. There are 5 lots of the set of 4 qPCR files; 3 file sets hold separate experimental replicates for a plate combining pRPS3 and pPGK1, the other 2 file sets are a sole pSRO9 run and a double pSRO9 run (so that 3 bioreps in total are available). 

### RPS3_SRO9_short_vs_long
Contains two experimental replicates for the short vs long terminator length assay. 
