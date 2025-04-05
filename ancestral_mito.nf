#! /usr/bin/env nextflow
nextflow.enable.dsl=2

if (!params.out){
    error("Output directory is required.")
}
if (!params.mito1 || !params.mito2 || !params.mito_outgroup){
    error("mito1, mito2, and mito_outgroup are required.")
}

process pad_seqs{
    time '1h'
    
    input:
    tuple path(mito1), path(mito2), path(mito3)
    
    output:
    tuple path("mito1_mod.fa"), path("mito2_mod.fa"), path("mito3_mod.fa")
    
    script:
    """
    samtools faidx ${mito1}
    samtools faidx ${mito2}
    samtools faidx ${mito3}
    
    if [ \$( cat ${mito1}.fai | wc -l ) -gt 1 ]; then
        >&2 echo "ERROR: ${mito1} contains multiple seqs"
        exit 1
    elif [ \$( cat ${mito2}.fai | wc -l ) -gt 1 ]; then
        >&2 echo "ERROR: ${mito2} contains multiple seqs"
        exit 1
    elif [ \$( cat ${mito3}.fai | wc -l ) -gt 1 ]; then
        >&2 echo "ERROR: ${mito3} contains multiple seqs"
        exit 1
    fi
    
    m1len=\$( cat ${mito1}.fai | head -1 | cut -f2 )
    m2len=\$( cat ${mito2}.fai | head -1 | cut -f2 )
    m3len=\$( cat ${mito3}.fai | head -1 | cut -f2 )
    
    m1estart=\$(( \$m1len - ${params.padding} ))
    m2estart=\$(( \$m2len - ${params.padding} ))
    m3estart=\$(( \$m3len - ${params.padding} ))

    m1name=\$( cat ${mito1}.fai | head -1 | cut -f1 )
    m2name=\$( cat ${mito2}.fai | head -1 | cut -f1 )
    m3name=\$( cat ${mito3}.fai | head -1 | cut -f1 )
    
    echo -e "\${m1name}\t\${m1estart}\t\${m1len}" | bedtools getfasta -bed stdin \
        -fi ${mito1} -fo m11.fa
    echo -e "\${m2name}\t\${m2estart}\t\${m2len}" | bedtools getfasta -bed stdin \
        -fi ${mito2} -fo m21.fa
    echo -e "\${m3name}\t\${m3estart}\t\${m3len}" | bedtools getfasta -bed stdin \
        -fi ${mito3} -fo m31.fa
    
    echo -e "\${m1name}\t0\t1000" | bedtools getfasta -bed stdin -fi ${mito1} \
        -fo m12.fa
    echo -e "\${m2name}\t0\t1000" | bedtools getfasta -bed stdin -fi ${mito2} \
        -fo m22.fa
    echo -e "\${m3name}\t0\t1000" | bedtools getfasta -bed stdin -fi ${mito3} \
        -fo m32.fa
    
    echo ">seq1" > mito1_mod.fa
    echo ">seq2" > mito2_mod.fa
    echo ">seq3" > mito3_mod.fa
    
    cat m11.fa | grep -v ">" | tr -d '\n' >> mito1_mod.fa
    cat m21.fa | grep -v ">" | tr -d '\n' >> mito2_mod.fa
    cat m31.fa | grep -v ">" | tr -d '\n'  >> mito3_mod.fa
    
    cat ${mito1} | grep -v ">" | tr -d '\n' >> mito1_mod.fa
    cat ${mito2} | grep -v ">" | tr -d '\n' >> mito2_mod.fa
    cat ${mito3} | grep -v ">" | tr -d '\n' >> mito3_mod.fa
    
    cat m12.fa | grep -v ">" | tr -d '\n' >> mito1_mod.fa
    cat m22.fa | grep -v ">" | tr -d '\n' >> mito2_mod.fa
    cat m32.fa | grep -v ">" | tr -d '\n' >> mito3_mod.fa
    """
}

process aln_padded{
    time '6h'
    
    input:
    tuple path(mito1), path(mito2), path(mito3)
    
    output:
    path("mito_aln.fa")
    
    script:
    """
    cat ${mito1} ${mito2} ${mito3} > mito_all.fa
    muscle -align mito_all.fa -output mito_aln.fa   
    """
}

process make_consensus{
    time '10m'
    
    input:
    path(aln)
    
    publishDir "${params.out}", mode: 'copy'
    
    output:
    path("chrM_consensus.fa")
    
    script:
    """
    ${baseDir}/scripts/make_consensus.py -a ${aln} -p ${params.padding} > chrM_consensus.fa
    """
}

workflow{
    
    seqs_padded = Channel.fromPath(params.mito1).combine(Channel.fromPath(params.mito2)).combine(Channel.fromPath(params.mito_outgroup)) | pad_seqs
    
    aln = aln_padded(seqs_padded)    
    make_consensus(aln)
}

