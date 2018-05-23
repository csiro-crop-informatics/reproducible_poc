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

#rmarkdown::render("show.Rmd", output_format = "all")
rmarkdown::render("show.Rmd", output_format = "revealjs::revealjs_presentation", output_file = "show.html" )
rmarkdown::render("show.Rmd", output_format = "html_document", output_file = "docu.html" )
#rmarkdown::render("show.Rmd", output_format = "slidy_presentation", output_file = "slidy.html" )
#rmarkdown::render("show.Rmd", output_format = "ioslides_presentation", output_file = "iso.html" )
#rmarkdown::render("show.Rmd", output_format = "beamer_presentation", output_file = "show.pdf")
