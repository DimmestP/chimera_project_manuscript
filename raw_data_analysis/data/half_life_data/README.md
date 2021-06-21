## Folder Contents

### chan_decay_data.txt
Contains the transcript wide decay rate dataset as measured by Chan et al 2018.

### sun_decay_data.txt
Contains the transcript wide decay rate dataset as measured by Sun et al 2013.

### collated_suspected_decay_motifs.txt
Contains all suspected RNA binding protein motifs as presented in Hogan et al 2008, Cheng et al 2017 and Shalgi et al 2005.

### scer_codon_counts_frequency_long.txt
A long-format tidy table of the count and proportion of each sense codon within each (non-frameshifted) S. cerevisiae ORF.
It has missing entires, not zero entires, if a codon is not present.
This is calculated by `count_codons_all_ORFs.R`.

### scer_codon_frequency_table.txt
A wide-format frequency table of the proportion of each sense codon within each (non-frameshifted) S. cerevisiae ORF.
This is calculated by `count_codons_all_ORFs.R`.
It contains the same frequencies as `scer_codon_counts_frequency_long.txt`, only reshaped into a wide table.

### whole_genome_3UTR.csv
Contains the median length 3'UTR sequence for all genes in the yeast genome according to Pelechano et al 2013
