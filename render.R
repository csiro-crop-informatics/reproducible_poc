#!/usr/bin/env Rscript

library(rmarkdown)
rmarkdown::render("doc.Rmd", output_format = "all")
rmarkdown::render("show.Rmd", output_format = "all")

