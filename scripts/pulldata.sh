#!/usr/bin/env bash

#REFSOURCE='ftp://ftp.ensemblgenomes.org/pub/plants/release-38/fasta/oryza_sativa/dna/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa.gz'
REFSOURCE='ftp://ftp.ensemblgenomes.org/pub/plants/release-38/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz'
REFBASE=${REFSOURCE##*/}
REFBASE=${REFBASE%.toplevel.fa.gz}.fa
REFDIR='data/reference'

mkdir -p ${REFDIR}
curl --progress-bar ${REFSOURCE} | gunzip --stdout > ${REFDIR}/${REFBASE}

