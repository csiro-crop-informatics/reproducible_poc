//READ SIMULATION PARAMS
seqerrs = params.seqerrs
//nreads = params.nreads
nreadsarr = [1000,2000,4000];
urls = params.urls

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
  tag {name}
  input:
    set val(name), val(url) from urls
    
  output:
    set val(name), file(ref) into kangaRefs, hisat2Refs, simReadsRefs

    """
    curl ${url} | gunzip --stdout > ref
    """
}

process kangaSimReads {


  tag {name+"_"+nreads}
  input:
    set val(name), file(ref) from simReadsRefs
    each nreads from nreadsarr
    
  output:
    set val(nametag), file ("*r1") into kangaR1, hisat2R1, fa2fqR1
    set val(nametag), file ("*r2") into kangaR2, hisat2R2, fa2fqR2

  script:
  nametag = name+"_"+nreads
  """
  biokanga simreads \
  --pegen \
  --seqerrs ${seqerrs} \
  --in ${ref} \
  --nreads ${nreads} \
  --out "${nametag}_r1" \
  --outpe "${nametag}_r2"
  """
}

process fasta2mockFASTQ {
  tag {nametag}
  input:
    set val(nametag), file("*r1*") from fa2fqR1
    set val(nametag), file("*r2*") from fa2fqR2
    
  output:
    set val(nametag), file ("*q1") into FASTQ1
    set val(nametag), file ("*q2") into FASTQ2
    
    """
    fasta2fastqDummy.sh *r1 > ${nametag}_r1_q1
    fasta2fastqDummy.sh *r2 > ${nametag}_r2_q2
    """

}

process fastQC {
  tag {nametag}

  input:
    set val(nametag), file(reads1) from FASTQ1 //R1
    set val(nametag), file(reads2) from FASTQ2 //R1

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
  tag{name}
  input:
    set val(name), file(ref) from hisat2Refs
  
  output:
    file 'hisat2db*' into hisat2dbs
    
    """
    hisat2-build ${ref} hisat2db -p 8
    """
}

process hisat2Align {
  tag {nametag}
  input:
    set val(nametag), file(r1) from hisat2R1
    set val(nametag), file(r2) from hisat2R2
    file hisat2db from hisat2dbs
  output:
    val nametag into samname
    stdout sam

    """
    hisat2 -x hisat2db -f -1 ${r1} -2 ${r2} 
    """
}

process sam2bam {
  publishDir "${params.outdir}/BAMs/hisat2", mode: 'link'  
  tag {nametag}
  input: 
    stdin sam
    val nametag from samname

  output:
    set val(nametag), file("*bam") into BAMs
    
    """
    samtools view -bS -F 4 -F 8 -F 256 > ${nametag}.bam
    """
}

process kangaIndex {
  tag{name}
  input:
    set val(name), file(ref) from kangaRefs

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
  publishDir "${params.outdir}/BAMs/biokanga", mode: 'link'
  tag {nametag}
  input:
    set val(nametag), file(r1) from kangaR1
    set val(nametag), file(r2) from kangaR2
    file kangadb

  output:
    set val(nametag), file("*bam") into kangaBAMs

    """
    biokanga align \
    -i ${r1} \
    -u ${r2} \
    --sfx ${kangadb} \
    -o ${nametag}_kanga.bam \
    --pemode 2 \
    --substitutions 3 
    """
}


//process bamReIndex {
//  tag {nametag}
//  publishDir "${params.outdir}/BAMs/biokanga", mode: 'link'
//  input:
//    set val(nametag),file("*bam") from bams
////    each tool from 
//    
//  output:
//    file "*.bai"
//    
//    """
//    samtools index *bam
//    """
//}


