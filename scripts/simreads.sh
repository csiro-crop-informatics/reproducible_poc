#!/usr/bin/env bash

IN='data/reference/Oryza_sativa.IRGSP-1.0.dna.fa'
SEQERRS=1.5
OUTDIR='data/simreads'

mkdir -p ${OUTDIR}
OUT=${OUTDIR}/${IN##*/}
OUT=${OUT%.fa}_ERR_${SEQERRS}

module load biokanga/4.3.9

biokanga simreads --pegen --seqerrs ${SEQERRS} --threads ${SLURM_CPUS_PER_TASK} \
--in ${IN} \
--out ${OUT}.R1.fa \
--outpe ${OUT}.R2.fa
