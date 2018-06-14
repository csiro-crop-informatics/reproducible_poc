//READ SIMULATION PARAMS
seqerrs = params.seqerrs.toString().tokenize(",")
nsimreadsarr = params.nsimreads.toString().tokenize(",")*.toInteger()
nrepeat = params.nrepeat
//INPUT GENOME PARAMS
url = params.url
name = params.name
//INPUT READS PARAMS
reads1url = params.realreads1
reads2url = params.realreads2

docheader = file(params.docheader)

def helpMessage() {
    log.info"""
    ===========================================================
    csiro-crop-informatics/reproducible_poc  ~  version ${params.version}
    ===========================================================
    Usage:

    nextflow run csiro-crop-informatics/reproducible_poc -r develop

    Default params:
    seqerrs     : ${params.seqerrs}
    nsimreads   : ${params.nsimreads} - this can be a comma-delimited list e.g. 100,20000,400
    nrepeat     : ${params.nrepeat} 
    url         : ${params.url}
    name        : ${params.name}
    outdir      : ${params.outdir}
    publishmode : ${params.publishmode} use copy or move if working across filesystems
    """.stripIndent()
}

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

/*
 * Create a channel for (local) input read files
 */
//Channel
//    .fromFilePairs( params.reads, size: 2 )
//    .ifEmpty { exit 1, "Cannot find reads matching: ${params.reads}\nNB: Path must contain at least one * wildcard and be enclosed in quotes." }
//    .set { local_read_files }

process fetchRef {
  tag {name}
  input:
    val url
    val name

  output:
    set val(name), file(ref) into kangaRefs, hisat2Refs, simReadsRefs

  script:
    """
    curl ${url} | gunzip --stdout | head -100000 > ref
    """
}


process fetchReads {

  input: 
    val reads1url
    val reads2url

  output:
    set val(longtag), val(nametag),file("r1.gz"), file("r2.gz") into FASTQ, hisat2FASTQ, kangaFASTQ

  script:
    nametag = "tmpTAG"
    longtag = ["name":"real", "nreads":"10k", "seqerr":"unk", "rep":"na", "format":"fq"]
    """
    curl ${reads1url} | gunzip --stdout | head -n 40000 | pigz --fast > r1.gz
    curl ${reads2url} | gunzip --stdout | head -n 40000 | pigz --fast > r2.gz
    """

}

process kangaSimReads {
  label 'biokanga'
  tag {longtag}
  input:
    set val(name), file(ref) from simReadsRefs
    each nsimreads from nsimreadsarr
    each seqerr from seqerrs
    each rep from 1..nrepeat

  output:
    set val(longtag), val(nametag),file("r1.gz"),file("r2.gz") into kangaReads, hisat2reads, fa2fqreads //simReads

  when:
    nsimreads > 0 
  
  script:
    nametag = name+"_"+nsimreads+"_"+seqerr+"_"+rep
    longtag = ["name":name, "nreads":nsimreads, "seqerr":seqerr, "rep":rep, "format":"fa"]
    """
        biokanga simreads \
        --pegen \
        --seqerrs ${seqerr} \
        --in ${ref} \
        --nreads ${nsimreads} \
        --out r1 \
        --outpe r2 \
        && pigz --fast r1 r2
    """
}

process fasta2mockFASTQ {
  tag {longtag}
  input:
    set val(longtag),val(nametag),file(r1),file(r2) from fa2fqreads

  output:
    set val(longtag), val(nametag), file ("*.q1.gz"), file("*.q2.gz") into MockFASTQ

    """
    zcat ${r1} | fasta2fastqDummy.sh | pigz --fast --stdout > "${nametag}.q1.gz"
    zcat ${r2} | fasta2fastqDummy.sh | pigz --fast --stdout > "${nametag}.q2.gz"
    """
}

process fastQC {
  tag {longtag}
  input:
    set val(longtag), val(nametag), file("${nametag}.q1.gz"), file("${nametag}.q2.gz") from MockFASTQ.mix(FASTQ)

  output:
    file "*_fastqc.{zip,html}" into fastqc_results

    """
    fastqc -q "${nametag}.q1.gz" "${nametag}.q2.gz"
    """
}

process multiQC {
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
  label 'hisat2'
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
  label 'hisat2'
  label 'samtools'
  tag {longtag}
  input:
    set val(longtag0), val(name), file(r1),file(r2) from hisat2reads.mix(hisat2FASTQ)
    set val(dbname), file("hisat2db.*.ht2") from hisat2dbs

  output:
    set val(longtag), val(tag), file("${tag}.bam") into hisat2BAMs

  script:
    tag = name+"_vs_"+dbname+".hisat2"
    longtag = longtag0.clone() //deepCopy(longtag0)
    longtag.ref = dbname
    longtag.aligner = "HISAT2"
    format = longtag["format"]=="fq"?"-q":"-f"
    """
    hisat2 -x hisat2db ${format} -1 ${r1} -2 ${r2} \
    | samtools view -bS -F 4 -F 8 -F 256 - > ${tag}.bam
    """
}

process kangaIndex {
  label 'biokanga'
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
  label 'biokanga'
  tag {longtag}
  input:
    set val(longtag0), val(name),file(r1),file(r2) from kangaReads.mix(kangaFASTQ)
    set val(dbname),file(kangadb) from kangadbs

  output:
    set val(longtag), val(tag), file("${tag}.bam") into kangaBAMs

  script:
    tag = name+"_vs_"+dbname+".biokanga"
    longtag = longtag0.clone() //otherwise modifying orginal map, triggering re-runs with -resume
    longtag.ref = dbname
    longtag.aligner = "BioKanga"
    """
    biokanga align \
    -i ${r1} \
    -u ${r2} \
    --sfx ${kangadb} \
    --threads ${task.cpus} \
    -o "${tag}.bam" \
    --pemode 2 \
    """
}

process extractStatsFromBAMs {
  label 'samtools'
  tag {longtag}
  input: 
    set val(longtag), val(nametag), file("${nametag}*.bam") from kangaBAMs.mix(hisat2BAMs)

  output:
    file statsFile into statsFiles
    val longtag into longtags

  script:
    statsPrefix = longtag.values().join("\t")+"\t"
    """
    echo -ne "${statsPrefix}" > statsFile
    samtools view ${nametag}.bam | extractStatsFromBAM.sh >> statsFile
    """
}

process combineStats {
  input:
    file("statsFile*") from statsFiles.collect()
    val longtag from longtags.first()
    
  output:
    file allStats into allStatsForFigs, allStatsForDoc
    
  script:
  statsHeader = longtag.keySet().join("\t")+"\t"+"Matches\tAlignments\tMatchRate"
    """
      cat <(echo -e "${statsHeader}") statsFile* >> allStats
    """
}

process MOCK_generateFigures {
  label "MOCK_PROCESS"
  input:
    file allStats from allStatsForFigs

  output:
    file("*.figure") into figures
    
  script: 
  """
    cat allStats > one.figure
    cat allStats > another.figure
  """
//    set file("*${nametag}.metadata"), file("*${nametag}.figure") into figures

//  script:
//    """
//    echo "${nametag}" > "${nametag}.metadata"
//    echo "${nametag}" > "${nametag}.figure"
//    """
}

process MOCK_generateReportMatter {
  label "MOCK_PROCESS"
  input:
    file "*figure" from figures.collect()
    file allStats from allStatsForDoc
    //set file(metadata), file(figure) from figures.collate(2)
    //set val(nametag), file(statsFile) from statsFiles.collate(2)
    file "*multiqc_report.html" from multiqc_report
    file "*_data" from multiqc_data
    file docheader

  script:
    """
    
    """
//    echo "---" > "${writeup}"
//    cat "${docheader}" >> "${writeup}"
//    echo -e "---\n" >> "${writeup}"
//    echo "# Stats\n" >> "${writeup}"
//    """
}


//process tagLocalReads {
//  input: 
//    set val(name),file(reads) from local_read_files
//  
//  output:
//    set val(longtag), val(nametag),file(r1), file(r2) into FASTQlocal //, hisat2FASTQlocal, kangaFASTQlocal

//  exec:
//    nametag = name
//    longtag = ["name":"local", "nreads":"unk", "seqerr":"unk", "rep":"na", "format":"fq"]
//}
