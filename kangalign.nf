params.r1 = "$PWD/data/simreads/Oryza_sativa.IRGSP-1.0.dna_ERR_1.5_NREADS_1000000.R1.fa"
params.r2 = "$PWD/data/simreads/Oryza_sativa.IRGSP-1.0.dna_ERR_1.5_NREADS_1000000.R2.fa"
params.ref = "$PWD/data/reference/Oryza_sativa.IRGSP-1.0.dna.fa"

r1 = file(params.r1)
r2 = file(params.r2)
ref = file(params.ref)

process kangaIndex {
  input:
  file ref

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
    file r1
    file r2
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

//process samToBam {
//    input: 
//    file sam
//    
//    output:
//    file outbam 
//    
//    """
//    samtools view -b sam > outbam
//    """
//}


//process parseLog {
//    input:
//    file kangalog
//    
//    output: 
//    stdout outchannel2
//    
//    """
//    head kangalog
//    """
//}

//outchannel1.subscribe { print "Indexing..  $it" }
//outchannel2.subscribe { print "Parselog..  $it" }
