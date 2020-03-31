# Modular effects of yeast 3'UTRs and cis-regulartory elements inside yeast 3'UTRs on mRNA abundance
The development of a manuscript investigating the relationship between mRNA 3'UTR motifs effecting half-lives and their modular effects on transcript abundance.

The data and code (including ALL analysis referred to in the paper) required to create the manuscript is included in this GitHub page. The R package Bookdown, created by Yihui Xie, can be used to recreate the manuscript and re-calculate all of the analysis. 

To recreate the manuscript:

``` 
    git clone <this repo>
    open -a RStudio

    # In console at the bottom of RStudio
    setwd( <local repo address> )
    source("./library_setup.R") # check all required libraries are installed
    bookdown::render_book("abstract.Rmd") # this will take a little while
```

The pdf version of the manuscript will be created in a _book folder.

