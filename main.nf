//READ SIMULATION PARAMS
seqerrs = params.seqerrs.toString().tokenize(",")
nreadsarr = params.nreads.toString().tokenize(",")
nrepeat = params.nrepeat
//INPUT GENOME PARAMS
url = params.url
name = params.name


//TODO add (conditional?) simulation of transcript reads using a feature file biokanga simreads --featfile 

def helpMessage() {
    log.info"""
    ===========================================================
    csiro-crop-informatics/reproducible_poc  ~  version ${params.version}
    ===========================================================
    Usage:
    
    nextflow run csiro-crop-informatics/reproducible_poc -r develop
    
    Default params:
    seqerrs    : ${params.seqerrs}
    nreads     : ${params.nreads} - this can be a comma-delimited set e.g. 100,20000,400
    nrepeat   : ${params.nrepeat} 
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
    each seqerr from seqerrs
    each rep from 1..nrepeat

  output:
    set val(nametag),file("r1.gz"),file("r2.gz") into kangaReads, hisat2reads, fa2fqreads //simReads 

  script:
    println(rep)
    println(nrepeat)
    nametag = nrepeat == 1 ? name+"_"+nreads : name+"_"+nreads+"_"+rep
    """
    biokanga simreads \
    --pegen \
    --seqerrs ${seqerr} \
    --in ${ref} \
    --nreads ${nreads} \
    --out r1 \
    --outpe r2 \
    && pigz --fast r1 r2
    """
}

process fasta2mockFASTQ {
  tag {nametag}
  input:
    set val(nametag),file(r1),file(r2) from fa2fqreads
    
  output:
    set val(nametag), file ("*.q1.gz"), file("*.q2.gz") into FASTQ
    
    """
    zcat ${r1} | fasta2fastqDummy.sh | pigz --fast --stdout > "${nametag}.q1.gz"
    zcat ${r2} | fasta2fastqDummy.sh | pigz --fast --stdout > "${nametag}.q2.gz"
    """
}

process fastQC {
  tag {nametag}

  input:
    set val(nametag), file("${nametag}.q1.gz"), file("${nametag}.q2.gz") from FASTQ 

  output:
    file "*_fastqc.{zip,html}" into fastqc_results

    """
    fastqc -q "${nametag}.q1.gz" "${nametag}.q2.gz"
    """


}

process multiQC {    
  publishDir "${params.outdir}/MultiQC", mode: 'copy'

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
  publishDir "${params.outdir}/BAMs/hisat2", mode: 'copy'
  tag {name+" vs "+dbname}
  input:
    set val(name), file(r1),file(r2) from hisat2reads
    set val(dbname), file("hisat2db.*.ht2") from hisat2dbs

  output:    
    set val(tag), file("${tag}.bam") into hisat2BAMs
//    val tag into samname
//    stdout sam into ssam

  script:
    tag = name+"_vs_"+dbname+".hisat2"
    """
    hisat2 -x hisat2db -f -1 ${r1} -2 ${r2} \
    | samtools view -bS -F 4 -F 8 -F 256 - > ${tag}.bam
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
  publishDir "${params.outdir}/BAMs/biokanga", mode: 'copy'
  tag {name+" vs "+dbname}
  input:
    set val(name),file(r1),file(r2) from kangaReads
    set val(dbname),file(kangadb) from kangadbs

  output:
    set val(tag), file("${tag}.bam") into kangaBAMs

  script:
    tag = name+"_vs_"+dbname+".biokanga"
    """
    biokanga align \
    -i ${r1} \
    -u ${r2} \
    --sfx ${kangadb} \
    --threads ${task.cpus} \
    -o "${tag}.bam" \
    --pemode 2 \
    --substitutions 3 
    """
}


process extractStatsFromBAMs {
  tag {nametag}
  input: 
    set val(nametag), file("${nametag}*.bam") from kangaBAMs.mix(hisat2BAMs)
  
  output:
    set val(nametag), file(statsFile) into statsFiles
  
//  exec:
//    println "Placeholder for extracting stats from ${nametag}" 
   
  script:
    """
    echo -n "${nametag}\t" > statsFile
    samtools view ${nametag}.bam | extractStatsFromBAM.sh >> statsFile
    """
}


//| awk -vOFS="\t" '{split($1,sim,"|");if(sim[4]==$3 && sim[5]==$4-1){count++}};END{print count,NR,count/NR}' >> statsFile

process MOCK_generateFigures {
  tag {nametag}
  label "MOCK_PROCESS"
  input: 
    set val(nametag), file(statsFile) from statsFiles
  
  output:
    file figure into figures
    
  script:
    """
    echo "${nametag}" > figure
    """
}


