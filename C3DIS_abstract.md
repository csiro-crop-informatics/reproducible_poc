Reproducible publication for data science
-----------------------------------------

## Introduction

Reproducibility of experiments is an essential part of the scientific method. Many scientific studies are difficult or impossible to replicate, so the correctness of their results cannot always be verified. In data science, the widespread use of often proprietary, graphical interfaces is generally detrimental to reproducibility of data processing and analyses. Replication attempts are hampered by the usually limited accuracy and completeness of description of steps undertaken. In contrast, the rise of free and open source software coupled with advances in computational process management, version control, containerisation and automation, brings reproducibility within reach.
  
## Methods 

BioKanga is a suite of bioinformatics tools developed at CSIRO. To evaluate BioKanga's sequence alignment module against other state-of-the-art tools, we developed workflows using two popular automation frameworks, Snakemake and Nextflow. Individual tasks are executed in Cloud or HPC environment using dedicated, lightweight (Docker/Singularity) containers for each tool being evaluated. Containers are also used for other tasks including extensive quality-control of the input as well as the interim data and final results. Use of containers provides a reproducible software environment which contributes to the replicability of the results. Given a pre-defined experimental set-up, raw data is acquired from a public source, processed and analysed. Final document including dynamically generated figures and tables is compiled from R Markdown or LaTeX/knitr. 

## Conclusion

This work is a proof of concept and provides an open source template for automated generation of a fully reproducible data science publication. 




