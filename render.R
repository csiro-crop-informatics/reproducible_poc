#!/usr/bin/env Rscript

#required modules if executed on the cluster: R/3.4.4 pandoc/1.12.3 

if(!require(rmarkdown)){
    install.packages("rmarkdown")
    library(rmarkdown)
}

if(!require(revealjs)){
    location <- "~/local/R_libs/"
    dir.create(location, recursive = TRUE) 
    install.packages("revealjs", lib=location, repos='https://cran.csiro.au')
    library(revealjs, lib.loc=location)
}

render("doc.Rmd")
render("show.Rmd")

