//READ SIMULATION PARAMS
seqerrs = params.seqerrs
//nreads = params.nreads
nreadsarr = [1000,2000,4000,8000];
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
  tag {nreads}
  input:
    file ref from refs
    each nreads from nreadsarr
    
  output:
    set val(nreads), file ("*r1") into kangaR1, hisat2R1, fa2fqR1
    set val(nreads), file ("*r2") into kangaR2, hisat2R2, fa2fqR2
//    file "*r1" into kangaR1, hisat2R1, fa2fqR1
//    file "*r2" into kangaR2, hisat2R2, fa2fqR2
//    val nreads into fnames


  """
  biokanga simreads \
  --pegen \
  --seqerrs ${seqerrs} \
  --in ${ref} \
  --nreads ${nreads} \
  --out "${nreads}_r1" \
  --outpe "${nreads}_r2"
  """
}

process fasta2mockFASTQ {
  tag {nreads}
  input:
    set val(nreads), file("*r1*") from fa2fqR1
    set val(nreads), file("*r2*") from fa2fqR2
//    val nreads from fnames
    
  output:
    set val(nreads), file ("*q1") into FASTQ1
    set val(nreads), file ("*q2") into FASTQ2
//    val nreads into fqnames
    
    """
    fasta2fastqDummy.sh *r1 > ${nreads}_r1_q1
    fasta2fastqDummy.sh *r2 > ${nreads}_r2_q2
    """

}

process fastQC {
  tag {nreads}
////  tag "$name"
////  publishDir "${params.outdir}/fastqc", mode: 'link',
////        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

  input:
    set val(nreads), file(reads1) from FASTQ1 //R1
    set val(nreads), file(reads2) from FASTQ2 //R1
    //file "*r2" from FASTQ2 //R2
//    val nreads from fqnames

  output:
    file "*_fastqc.{zip,html}" into fastqc_results

    """
    fastqc -q $reads1 $reads2 //*r1 *r2
    """


}

process multiQC {    
  publishDir "${params.outdir}/MultiQC", mode: 'link'

  input: 
    file f from fastqc_results.collect()
  
  output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"
    
    """
    pwd
    multiqc . -f
    """
}

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
  tag {nreads}
  input:
    set val(nreads), file(r1) from hisat2R1
    set val(nreads), file(r2) from hisat2R2
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
  publishDir "${params.outdir}/hisat2BAMs", mode: 'link'
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
  tag {nreads}
  publishDir "${params.outdir}/kangaBAMs", mode: 'link'
  tag {nreads}
  input:
    set val(nreads), file(r1) from kangaR1
    set val(nreads), file(r2) from kangaR2
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
  publishDir "${params.outdir}/kangaBAMs", mode: 'link'
  input:
    file bam from bams
    
    """
    samtools index $bam
    """
}


