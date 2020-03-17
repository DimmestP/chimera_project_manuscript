# Prediction of cell-to-cell variability with stability motifs
The development of a manuscript investigating the relationship between mRNA transcript half-lives and cell-to-cell variability in transcript numbers.

The data and code (including ALL analysis referred to in the paper) required to create the manuscript is included in this GitHub page. The R package Bookdown, created by Yihui Xie, can be used to recreate the manuscript and re-calculate all of the analysis. 

To recreate the manuscript:

``` 
    git clone <this repo>
    open -a RStudio

    # In R console
    setwd( " <local repo address> " )
    install.packages("bookdown") # if necessary
    bookdown::render_book("abstract.Rmd")
```

The pdf version of the manuscript will be created in a _book folder.

