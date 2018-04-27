params.url = "ftp://ftp.ensemblgenomes.org/pub/plants/release-38/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"
url = params.url

//READ SIMULATION PARAMS
params.seqerrs = 1.5
params.nreads = 100000
seqerrs = params.seqerrs
nreads = params.nreads


process fetchRef {
  input:
    val url
    
  output:
    file ref into refs

    """
    curl ${url} | gunzip --stdout > ref
    """
}

process kangaSimReads {
  input:
    file ref from refs
  
  output:
    file r1 into simreads1
    file r2 into simreads2
  
  """
  biokanga simreads \
  --pegen \
  --seqerrs ${seqerrs} \
  --in ${ref} \
  --nreads ${nreads} \
  --out r1 \
  --outpe r2
  """
}

process kangaIndex {
  input:
    file ref from refs

  output:
    file db

    """
    biokanga index \
    -i ${ref} \
    -o db \
    --ref ${ref}
    """
}

process kangaAlign {
  input:
    file r1 from simreads1
    file r2 from simreads2
    file db

  output:
    file 'out.bam' into bams

    """
    biokanga align \
    -i ${r1} \
    -u ${r2} \
    --sfx ${db} \
    -o out.bam \
    --pemode 2 \
    --substitutions 3 
    """
}


process bamReIndex {
  input:
    file bam from bams
    
    """
    samtools index $bam
    """
}


