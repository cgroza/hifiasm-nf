params.design = 'design.csv'
params.ref = 'genome.fa.gz'


design = Channel.fromPath(params.design).splitCsv(header : true).multiMap{
    row ->
    mat: [row.sample, "maternal", file(row.maternal_1), file(row.maternal_2)]
    pat: [row.sample, "paternal", file(row.paternal_1), file(row.paternal_2)]
    hifi: [row.sample, file(row.hifi)]
}

ref_ch = Channel.fromPath(params.ref)

yac_fastq_ch = design.mat.concat(design.pat)

process yak_kmers {
    cpus 40
    memory '170GB'
    time '12h'

    input:
    set val(sample), val(parent), file(fastq_1), file(fastq_2) from yac_fastq_ch

    output:
    set val(sample), file("${parent}.yak") into yac_kmers_ch

    script:
    """
    yak count -b37 -t40 -o ${parent}.yak <(cat ${fastq_1} ${fastq_2}) <(cat ${fastq_1} ${fastq_2})
    """
}

hifiasm_ch = design.hifi.join(yac_kmers_ch.groupTuple(by: 0).map{
    it -> [it[0], it[1][0], it[1][1]]}, by: 0).view()

process hifiasm_denovo {
    cpus 40
    memory '170GB'
    time '24h'

    publishDir 'assemblies'

    input:
    set val(sample), file(hifi_fasta), file(yak1), file(yak2) from hifiasm_ch

    output:
    set val(sample), file("${sample}_diploid") into hifiasm_denovo_ch

    script:
    """
    module load samtools
    module load bcftools
    mkdir ${sample}_asm
    hifiasm -o ${sample}_asm/${sample}.asm -t40 -1 ${yak1} -2 ${yak2} ${hifi_fasta}

    mkdir ${sample}_diploid
    awk '/^S/{print ">"\$2;print \$3}' ${sample}_asm/${sample}.asm.dip.hap1.p_ctg.gfa | bgzip > ${sample}_diploid/${sample}_hap1.fa.gz
    awk '/^S/{print ">"\$2;print \$3}' ${sample}_asm/${sample}.asm.dip.hap2.p_ctg.gfa | bgzip > ${sample}_diploid/${sample}_hap2.fa.gz
    samtools faidx ${sample}_diploid/${sample}_hap1.fa.gz
    samtools faidx ${sample}_diploid/${sample}_hap2.fa.gz

    rm -r ${sample}_asm
    """
}

process dipcall_variants {
    cpus 2
    memory '20GB'
    time '6h'
    publishDir 'variants'

    input:
    set val(sample), file(asm), file(ref) from hifiasm_denovo_ch.combine(ref_ch)

    output:
    set val(sample), file("${sample}.dip.vcf.gz")

    script:
    """
    module load samtools

    samtools faidx ${ref}

    ~/dipcall.kit/run-dipcall ${sample} ${ref} ${asm}/${sample}_hap1.fa.gz ${asm}/${sample}_hap2.fa.gz > ${sample}.mak
    make -j2 ${sample}.mak
    """
}
