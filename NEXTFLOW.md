Job execusion using [![Nextflow](https://www.nextflow.io/img/nextflow2014_no-bg.png)](https://www.nextflow.io/)
--------------------------------------------------------

We use [nextflow](https://www.nextflow.io/) to handle compute. One way to make nextflow available on your system: 

`curl -s https://get.nextflow.io | bash && mkdir -p ~/bin && mv nextflow ~/bin && PATH+=":~/bin"`

Current set-up executes a test pipeline composed of three processes. Alternative ways of running the pipeline include:

* starndard `nextflow kangalign.nf` - the required software, i.e. [biokanga](https://github.com/csiro-crop-informatics/biokanga) is assumed to be available 
* in container(s) using docker `nextflow kangalign.nf -profile docker`
* in container(s) using singularity `nextflow kangalign.nf -profile singularity`
* on a SLURM cluster `nextflow kangalign.nf -profile slurm`
* in container(s) using singularity on a SLURM cluster `module load singularity && kangalign.nf -profile slurm,singularity`


