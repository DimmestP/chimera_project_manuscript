# Modular effects of yeast 3'UTRs and cis-regulartory elements inside yeast 3'UTRs on mRNA abundance
The development of a manuscript investigating the relationship between mRNA 3'UTR motifs effecting half-lives and their modular effects on transcript abundance.

The data and code (including ALL analysis referred to in the paper) required to create the manuscript is included in this GitHub page. The R package Bookdown, created by Yihui Xie, can be used to recreate the manuscript and re-calculate all of the analysis. 
## Getting Started


### Making the Manuscript

``` 
    git clone <this repo>
    open -a RStudio

    # In console at the bottom of RStudio
    setwd( <local repo address> )
    bookdown::render_book("abstract.Rmd") # this will take a little while
    
    # If you get an error similar to 
    Error: Input files not all in same directory, please supply explicit wd
    # then run this code in your console and and try to render book again
    options(bookdown.render.file_scope = FALSE)

    # The pdf version of the manuscript will be created in a _book folder.
```
### Using Latex

The way bookdown works is that it runs all the code in all files, catches any code/graphs outputed then inserts it all into one massive latex document. Therefore, all of the text between code chunks is parsed through a latex compiler, so all latex functionality is availiable. There are many good latex tutorials out there but I recommend (this)[https://www.youtube.com/playlist?list=PLnC5h3PY-znyDQKn3knfXfekZLgWyL7QW] YouTube series for beginners, (this)[https://www.overleaf.com/learn/latex/Main_Page] wiki for more details and finally (this)[https://wch.github.io/latexsheet/] cheatsheet to refresh your memory on syntax.

## Structure
There are five main .Rmd files that hold the chapter text and analysis code;


## Tips 

When adding new lines of code, especially for importing data, remember it will be ran from the home directory (where ever abstract.Rmd) is not where your current *.Rmd is.

Manage citations using (Bibtex)[http://www.bibtex.org/] and the chimeraProject.bib file. Insert citations using [@<citeCode>] where needed.



