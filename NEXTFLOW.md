Job execusion using [![Nextflow](https://www.nextflow.io/img/nextflow2014_no-bg.png)](https://www.nextflow.io/)
--------------------------------------------------------

We use [nextflow](https://www.nextflow.io/) to handle compute. One way to make nextflow available on your system: 

`curl -s https://get.nextflow.io | bash && mkdir -p ~/bin && mv nextflow ~/bin && PATH+=":~/bin"`

Current set-up allows for a test script to be executed 

* locally `nextflow kangalign.nf`
* on a SLURM cluster `nextflow kangalign.nf -profile cluster`
* in container(s) (using singularity) `nextflow kangalign.nf -profile containers`
* using singularity on a SLURM cluster `module load singularity && kangalign.nf -profile cluster,containers`


