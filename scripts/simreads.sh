#!/usr/bin/env bash

#IN='data/reference/Oryza_sativa.IRGSP-1.0.dna.fa'
IN='data/reference/Arabidopsis_thaliana.TAIR10.dna.fa'
SEQERRS=1.5
OUTDIR='data/simreads'
NREADS=100000


mkdir -p ${OUTDIR}
OUT=${OUTDIR}/${IN##*/}
OUT=${OUT%.fa}_ERR_${SEQERRS}_NREADS_${NREADS}

module load biokanga/4.3.9

biokanga simreads \
--pegen \
--seqerrs ${SEQERRS} \
--threads ${SLURM_CPUS_PER_TASK:-1} \
--in ${IN} \
--nreads ${NREADS} \
--out ${OUT}.R1.fa \
--outpe ${OUT}.R2.fa 
