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

//process fasta2mockFASTQ {
//  input:
//    file r1 from R1
//    file r2 from R2
//    
//  output:
//    file q1 into FASTQ1
//    file q2 into FASTQ2
//    
//    """
//    sed 's/^>/@/' \$r1 | awk -vOFS="\n" '{print $1,$2,"+";gsub(".","a",$2);print $2}' > q1
//    sed 's/^>/@/' \$r2 | awk -vOFS="\n" '{print $1,$2,"+";gsub(".","a",$2);print $2}' > q2
//    """

//}

//process fastQC {

////  tag "$name"
////  publishDir "${params.outdir}/fastqc", mode: 'link',
////        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

//  input:
//    file r1 from FASTQ1 //R1
//    file r2 from FASTQ2 //R2


//  output:
//    file "*_fastqc.{zip,html}" into fastqc_results

//    """
//    fastqc -q ${r1} ${r2}
//    """


//}

//process multiQC {
//  input: 
//  
//  output:
//    file 
//}

process hisat2Index {
  input:
    file ref from refs
  
  output:
    file 'hisat2db*' into hisat2dbs
    
    """
    hisat2-build ${ref} hisat2db -p 8
    """
}

process hisat2Align {
  input:
    file r1 from R1
    file r2 from R2
    file hisat2db from hisat2dbs
  output:
    stdout sam 

    """
    hisat2 -x hisat2db -f -1 ${r1} -2 ${r2} 
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

process sam2bam {
  input: 
    stdin sam

  output:
    file bam
    
    """
    samtools view -bS -F 4 -F 8 -F 256 > bam
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


