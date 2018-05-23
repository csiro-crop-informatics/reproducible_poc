//READ SIMULATION PARAMS
seqerrs = params.seqerrs
nreadsarr = params.nreads.toString().tokenize(",")
//INPUT GENOME PARAMS
url = params.url
name = params.name
rmarkdownfile = file(params.rmarkdown)

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
    url        : ${params.url}
    name       : ${params.name}
    outdir     : ${params.outdir}
    rmarkdown  : ${params.rmarkdown}
    """.stripIndent()
}

// Show help message
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

  output:
    set val(name), file(ref) into kangaRefs, hisat2Refs, simReadsRefs

    """
    curl ${url} | gunzip --stdout > ref
    """
}

process kangaSimReads {
  tag {nametag}
  input:
    set val(name), file(ref) from simReadsRefs
    each nreads from nreadsarr

  output:
    set val(nametag),file("r1.gz"),file("r2.gz") into kangaReads, hisat2reads, fa2fqreads //simReads

  script:
    nametag = name+"_"+nreads
    """
    biokanga simreads \
    --pegen \
    --seqerrs ${seqerrs} \
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
    file "*_data" into multiqc_data

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

process MOCK_extractStatsFromBAMs {
  tag {nametag}
  label "MOCK_PROCESS"
  input:
    set val(nametag), file("${nametag}*.bam") from kangaBAMs.mix(hisat2BAMs)

  output:
    set val(nametag), file(statsFile) into statsFiles, statsFilesForFigures

//  exec:
//    println "Placeholder for extracting stats from ${nametag}"

  script:
    """
    echo "${nametag}" > statsFile
    """

}

process MOCK_generateFigures {
  tag {nametag}
  label "MOCK_PROCESS"
  input:
    set val(nametag), file(statsFile) from statsFilesForFigures

  output:
    set file(metadata), file(figure) into figures

  script:
    """
    echo "${nametag}" > figure
    echo "${nametag}" > metadata
    """
}

process generateReport {
  publishDir "${params.outdir}", mode: 'copy'
  
  input:
    set val(nametag), file(statsFile) from statsFiles
    set val(metadata), file(figure) from figures
    file "*multiqc_report.html" from multiqc_report
    file "*_data" from multiqc_data

  output:
    file "show.html"
    file "docu.html"

  script:
    """
    #!/usr/bin/env Rscript

    #required modules if executed on the cluster: R/3.4.4 pandoc/1.12.3

    #Install packages if absent
    if(!require(rmarkdown)){
        install.packages("rmarkdown")
        library(rmarkdown)
    }
    if(!require(revealjs)){
        location <- "~/local/R_libs/"
        dir.create(location, recursive = TRUE)
        install.packages("revealjs", lib=location, repos='https://cran.csiro.au')
        library(revealjs, lib.loc=location)
    }
    rmarkdown::render("${rmarkdownfile}", output_format = "revealjs::revealjs_presentation", output_file = "show.html" )
    """
}
