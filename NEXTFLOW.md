Job execution using [![Nextflow](https://www.nextflow.io/img/nextflow2014_no-bg.png)](https://www.nextflow.io/)
--------------------------------------------------------

We use [nextflow](https://www.nextflow.io/) to handle compute. One way to make nextflow available on your system: 

`curl -s https://get.nextflow.io | bash && mkdir -p ~/bin && mv nextflow ~/bin && PATH+=":~/bin"`

Current set-up executes a simple test pipeline composed of the following processes

* `curl` 
  * assembly download
* [biokanga](https://github.com/csiro-crop-informatics/biokanga) 
  * reference indexing 
  * short read simulation
  * short read alignment 
* [samtools](http://www.htslib.org/) 
  * (re-)indexing of the output BAM file 


The pipeline logic is largely separated from the execution environment configuration. 
Currentl to keep the set-up modular, we opt for container per tool, but if required for efficiency these could potentially be packaged in a single container.
Below, we list several alternative ways of executing the same pipeline. 

### Alternative ways of running the pipeline

Note, for any of the below using `-resume` prevents processes from being re-run if nighter input no scripts have changed.

Local/interactive with required software assumed to be available:

``` nextflow pipeline.nf```

Local/interactive with required software loaded via `module` directives:

```nextflow pipeline.nf -profile modules```

In container(s) using docker:

```nextflow pipeline.nf -profile docker```


In container(s) using singularity:

```nextflow pipeline.nf -profile singularity```

On a SLURM cluster with modules:

```nextflow pipeline.nf -profile slurm,modules```

* `slurm` profile sets some SLRUM defaults and ensures processes are submitted using `sbatch`
* `modules` profile facilitates loading of required software modules

In container(s) using singularity on a SLURM cluster, first ensuring singularity is avaiable on head/login node


```
module load singularity
nextflow pipeline.nf -profile slurm,singularity,singularitymodule
```

* `slurm` profile sets some SLRUM defaults and ensures processes are submitted using `sbatch`
* `singularity` profile facilitates pulling and converting appropriate docker images 
* `singularitymodule` profile ensures singularity module is loaded on each execution node

### Execution summary report

After execution of the pipeline a summary [report](report.html) is crated in the working directory.


