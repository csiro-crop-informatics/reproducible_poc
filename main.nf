//READ SIMULATION PARAMS
seqerrs = params.seqerrs
nreadsarr = params.nreadsarr
url = params.url
name = params.name

def helpMessage() {
    log.info"""
    ===========================================================
    csiro-crop-informatics/reproducible_poc  ~  version ${params.version}
    ===========================================================
    Usage:
    
    nextflow run csiro-crop-informatics/reproducible_poc -r develop
    
    Default params:
    seqerrs    : ${params.seqerrs}
    nreadsarr  : ${params.nreadsarr}
    url        : ${params.url}
    name       : ${params.name}
    outdir     : ${params.outdir}
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
    val url
    val name
    //set val(name), val(url) from namedurl
    
  output:
    set val(name), file(ref) into kangaRefs, hisat2Refs, simReadsRefs

    """
    #head -1000 ${baseDir}/data/ref1k > ref
    curl ${url} | gunzip --stdout > ref
    """
}

process kangaSimReads {
  tag {nametag}
  input:
    set val(name), file(ref) from simReadsRefs
    each nreads from nreadsarr
    
  output:
    set val(nametag),file(r1),file(r2) into kangaReads, hisat2reads, fa2fqreads

  script:
  nametag = name+"_"+nreads
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

process fasta2mockFASTQ {
  tag {nametag}
  input:
    set val(nametag),file(r1),file(r2) from fa2fqreads
    
  output:
    set val(nametag), file ("*.q1"), file("*.q2") into FASTQ
    
    """
    fasta2fastqDummy.sh r1 > "${nametag}_r1.q1"
    fasta2fastqDummy.sh r2 > "${nametag}_r2.q2"
    """
}

process fastQC {
  tag {nametag}

  input:
    set val(nametag), file(reads1), file(reads2) from FASTQ 

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
    set val(name), file("hisat2db.*.ht2") into hisat2dbs

    
    """
    hisat2-build ${ref} hisat2db -p 8
    """
}

process hisat2Align {
  tag {name+" vs "+dbname}
  input:
    set val(name), file(r1),file(r2) from hisat2reads
    set val(dbname), file("hisat2db.*.ht2") from hisat2dbs

  output:
    val tag into samname
    stdout sam

  script:
    tag = name+"_vs_"+dbname
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
    set val(name), file(kangadb) into kangadbs
    
    """
    biokanga index \
    -i ${ref} \
    -o kangadb \
    --ref ${ref}
    """
}

process kangaAlign {
  publishDir "${params.outdir}/BAMs/biokanga", mode: 'link'
  tag {nametag+" vs "+dbname}
  input:
    set val(nametag),file(r1),file(r2) from kangaReads
    set val(dbname),file(kangadb) from kangadbs

  output:
    set val(nametag), file("${nametag}_vs_${dbname}.bam") into kangaBAMs

    """
    biokanga align \
    -i ${r1} \
    -u ${r2} \
    --sfx ${kangadb} \
    -o "${nametag}_vs_${dbname}.bam" \
    --pemode 2 \
    --substitutions 3 
    """
}





