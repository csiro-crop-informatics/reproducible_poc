params.r1 = "$PWD/data/simreads/Oryza_sativa.IRGSP-1.0.dna_ERR_1.5_NREADS_1000000.R1.fa"
params.r2 = "$PWD/data/simreads/Oryza_sativa.IRGSP-1.0.dna_ERR_1.5_NREADS_1000000.R2.fa"
params.db = "$PWD/data/reference/Oryza_sativa.IRGSP-1.0.dna.sfx"
params.out = "$PWD/partial/aligned/Oryza_sativa.IRGSP-1.0.dna_ERR_1.5_NREADS_1000000.bam"

db = file(params.db)
r1 = file(params.r1)
r2 = file(params.r2)
out = file(params.out)

//biokanga requires .bam extension otherwise produces plain text SAM
 

process kangaAlign {
    input:
    file r1
    file r2

    output:
    file kangalog 

    kanga = file("kanga.bam")
    
    """
    time biokanga align \
    -i ${r1} \
    -u ${r2} \
    --sfx ${db} \
    -o ${out} \
    --log kangalog \
    --pemode 2 \
    --substitutions 3 
    

    """
}

process parseLog {
    input:
    file kangalog
    
    output: 
    stdout outchannel
    
    """
    head kangalog
    """
}

outchannel.subscribe { print "I say..  $it" }
