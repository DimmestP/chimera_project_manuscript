# Probing 3'UTRs as Modular Regulators of Gene Expression

Modular effects of yeast 3'UTRs and cis-regulatory elements on mRNA abundance.
This repository contains the data and analysis files for the manuscript:

> Probing 3'UTRs as Modular Regulators of Gene Expression.
> Jamie Auxillos, Samuel Haynes, Abhishek Jain, Clemence Alibert, Weronika Danecka, Edward W.J. Wallace

The data and code (including ALL analysis referred to in the paper) required to create the manuscript is included in this GitHub page. The R package [Bookdown](https://bookdown.org/) can be used to recreate the manuscript and a makefile is provide to re-calculate most of the analysis. 

## Overview 

There are three main stages to the analysis code presented in this repo, there are;

- Manuscript text and formatting

- Raw data analysis and normalisation

- Linear model fitting and results plotting

### Manuscript text and formatting

All of the text, formatting code, processed data and figures relating to the manuscript are presented in the home folder of the repository.
It is entirely self-contained from the analysis files producing the processed data/figures which are found in the raw_data_analysis folder (explored below).
The main text of the manuscript is held in five folders, denoted by chapter title; 

- intro_chapter

- results_chapter

- methods_chapter

- discussion_chapter 

- supplimentary_data_chapter

Each of these folders is self-contained; each chapter can be rendered separately and does not call on any file outside that folder. 
They hold an .Rmd file containing the text of that chapter and folders of all of the figures/data presented in that chapter.
Each data/figure folder also contains a README file that explains which analysis file in the raw_data_analysis folder recreates that figure/data.
Other than the text held in the chapter folders, the final sections of the manuscript found in the home repository are; 

- abstract.Rmd

- author.tex 

- references/chimeraProject.bib

Finally, For rendering the entire manuscript the _bookdown.yaml file tells the bookdown::render_book function where each chapter file is and the formatting folder holds the latex files hold the code to render the manuscript in the biorxiv format.

### Raw data analysis and normalisation
The raw_data_analysis folder contains all the raw data, analysis code and figure making code used in the manuscript.
It is split into code and data folders.
The data folder contains all of the raw data, normalised data and intermediate results produced/used by the analysis code.
The code folder contains all the R, python Jupyter notebook and Rmd files required to provide the results in the manuscript.
Both of these folders are further split by whether they contain/used platereader or qpcr experimental results (See respective folder README files for more details).
Apart from the Omniplate python code provided, all other libraries need to be downloaded by the user.
The sessionInfo.txt file in the repo's home folder contains all the libraries (and their versions) used to run the analysis code to help the user check they have the required R libraries.

### Linear model fitting and results plotting
Apart from figures directly plotting the normalised qpcr/platereader results that are contained in the respective raw analysis code file, the majority of the figure plotting and detailed data analysis files can be found in the raw_data_analysis/code folder. 
See the folder README file for full descriptions of what each file analysis and which figures it plots.
If you want to rerun the analysis files yourself they need to be ran in a certain order to provide dependences, please see the makefile or README file for more details.

## Getting Started


### Making the Manuscript

- Download file locally

``` 
    git clone https://github.com/DimmestP/chimera_project_manuscript.git
```
    
- Open RStudio

- Click File -> Open Project and find the 'chimera_project_manuscript.Rproj' file

- Make manuscript
```
    # In console at the bottom of RStudio
    install.packages(bookdown)
    
    bookdown::render_book("abstract.Rmd")

    # The pdf version of the manuscript will be created in a _book folder.
```
All of the figures required to make the manuscript from chapter Rmd file are provided in the repo. 
However, if you want to rerun the analysis from raw data you can use the makefile to automatically run the analysis in the right order.

- Check you have all the required R libraries using the renv package. 
```
# In RStudio Command line

install.packages("renv")

renv::restore()
```

- Use GNU make to run R scripts (GNU make is available by default in Mac/Linux terminals but you'll need to download a GNU terminal for Windows)

- Open terminal

```
	cd <chimera_project_manuscript folder>
	make
```

This will take a few minutes. All analysis files except the jupyter notebook files used to analyse the platereader data are ran. (Platereader analysis has been left to run manually because they can take hours).

### Using Latex

The way bookdown works is that it runs all the code in all files, catches any code/graphs outputed then inserts it all into one massive latex document. Therefore, all of the text between code chunks is parsed through a latex compiler, so all latex functionality is availiable. There are many good latex tutorials out there but I recommend (this)[https://www.youtube.com/playlist?list=PLnC5h3PY-znyDQKn3knfXfekZLgWyL7QW] YouTube series for beginners, (this)[https://www.overleaf.com/learn/latex/Main_Page] wiki for more details and finally (this)[https://wch.github.io/latexsheet/] cheatsheet to refresh your memory on syntax.

## Folder Contents

### intro_chapter, results_chapter, materials_and_methods_chapter, discussion_chapter, supplementary_data_chapter
Folders containing the text and data/figures required to create that chapter. Please read [above](#Manuscript text and formatting).

### formatting
Folder containing the TeX files required for LaTeX to render an article using the biorxiv template.

### references
Folder containing a bib file that holds the citing information for all articles reference in the manuscript.

### raw_data_analysis
Folder containing all the data and analysis code described in the manuscript. Please read [above](#Raw data analysis and normalisation).

### .gitignore
File used by git to ensure any temporary files created locally are not detected by the version control software (avoids errors in uploading wrong files to github)

### _bookdown.yml
File used by the bookdown R library to completely render the manuscript with all chapters figures references.

### abstract.Rmd 
File containing the text present in the abstract of the manuscript. It is also the anchor file required by bookdown to detect where all the chapter files are located according to the .yml file.

### author.tex
Contact and institution information for all the authors of the paper.

### chimera_project_manuscript.Rproj
Key file used by RStudio to automatically set up the working directory for R console (and enable all relative file paths to work)

### sessionInfo.txt
Useful R library information to aid with running analysis scripts on other computers. It is a list of all active libraries and their versions used the last time the script ran successfully.

## Tips 

When adding new lines of code, especially for importing data, remember it will be ran from the home directory (where ever abstract.Rmd) is not where your current *.Rmd is.

Manage citations using (Bibtex)[http://www.bibtex.org/] and the chimeraProject.bib file. Insert citations using [@<citeCode>] where needed.



