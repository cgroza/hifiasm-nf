params.design = 'design.csv'


design = Channel.fromPath(params.design).splitCsv(header : true).multiMap{
    row ->
    mat: [row.sample, "maternal", file(row.maternal_1), file(row.maternal_2)]
    pat: [row.sample, "paternal", file(row.paternal_1), file(row.paternal_2)]
    hifi: [row.sample, file(row.hifi)]
}

yac_fastq_ch = design.mat.concat(design.pat)

process yak_kmers {
    cpus 16
    memory '170GB'
    time '12h'

    input:
    set val(sample), val(parent), file(fastq_1), file(fastq_2) from yac_fastq_ch

    output:
    set val(sample), file("${parent}.yak") into yac_kmers_ch

    script:
    """
    yak count -b37 -t16 -o ${parent}.yak <(cat ${fastq_1} ${fastq_2}) <(cat ${fastq_1} ${fastq_2})
    """
}

hifiasm_ch = design.hifi.join(yac_kmers_ch.groupTuple(by: 0).map{
    it -> [it[0], it[1][0], it[1][1]]}, by: 0).view()

process hifiasm_denovo {
    cpus 40
    memory '170GB'
    time '48h'

    publishDir 'assemblies'

    input:
    set val(sample), file(hifi_fasta), file(yak1), file(yak2) from hifiasm_ch

    output:
    file("${sample}_asm") into hifiasm_denovo_ch

    script:
    """
    mkdir ${sample}_asm
    hifiasm -o ${sample}_asm/${sample}.asm -t40 -1 ${yak1} -2 {yak2} ${hifi_fasta}
    """
}
