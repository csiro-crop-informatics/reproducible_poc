//READ SIMULATION PARAMS
seqerrs = params.seqerrs
nreads = params.nreads
url = params.url

def helpMessage() {
    log.info"""
    ===========================================================
    csiro-crop-informatics/reproducible_poc  ~  version ${params.version}
    ===========================================================
    Usage:
    
    nextflow run csiro-crop-informatics/reproducible_poc 
    
    """.stripIndent()
}

// Show help emssage
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

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
    file r1 into R1
    file r2 into R2
  
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
    file kangadb

    """
    biokanga index \
    -i ${ref} \
    -o kangadb \
    --ref ${ref}
    """
}

process hisat2Index {
  input:
    file ref from refs
  
  output:
    file hisat2db 
    
    """
    hisat2 -h > hisat2db
    """
}

process hisat2Align {
  input:
    file r1 from R1
    file r2 from R2
    file hisat2db


    """
    hisat2 -h 
    """
    
//    hisat2 -x $index_base \\
//      -1 ${reads[0]} \\
//      -2 ${reads[1]} \\
//      $rnastrandness \\
//      --known-splicesite-infile $alignment_splicesites \\
//      --no-mixed \\
//      --no-discordant \\
//      -p ${task.cpus} \\
//      --met-stderr \\
//      --new-summary \\
//      --summary-file ${prefix}.hisat2_summary.txt \\
//      | samtools view -bS -F 4 -F 8 -F 256 - > ${prefix}.bam
//    """
}

process kangaAlign {
  input:
    file r1 from R1
    file r2 from R2
    file kangadb

  output:
    file 'out.bam' into bams

    """
    biokanga align \
    -i ${r1} \
    -u ${r2} \
    --sfx ${kangadb} \
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


