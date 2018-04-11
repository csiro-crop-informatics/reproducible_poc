#!/usr/bin/env bash
#SBATCH --time=2:00:00
#SBATCH --job-name=kangalign
#SBATCH --nodes=1
#SBATCH --cpus-per-task 4
#SBATCH --mem=8G 
#SBATCH --output jobs/%j_kangalign.out

module load biokanga/4.3.9

INDEX=${3%.fa}.sfx
if [[ ! -f ${INDEX} ]]; then
    biokanga index -i ${3} -o ${INDEX} -r blah -t ${SLURM_CPUS_PER_TASK:-1}
fi 

OUTDIR='partial/aligned'
mkdir -p ${OUTDIR} 
OUT=${OUTDIR}/${1##*/}_vs_${INDEX##*/}

echo "Kanga aligning ${1} ${2} vs ${3}"

time biokanga align \
-i ${1} \
-u ${2} \
--sfx ${INDEX} \
-o ${OUT}.bam \
--log ${OUT}.log \
--pemode 2 \
--threads ${SLURM_CPUS_PER_TASK:-1} \
--substitutions 3
