params.design = 'design.csv'
params.ref = 'genome.fa.gz'
params.hifi_parents = true
params.dipcall = true

ref_ch = Channel.fromPath(params.ref)

if(params.hifi_parents) {
    design = Channel.fromPath(params.design).splitCsv(header : true).multiMap {
        row ->
        mat: [row.sample, "maternal", file(row.maternal)]
        pat: [row.sample, "paternal", file(row.paternal)]
        hifi: [row.sample, file(row.hifi)]
    }
} else {
    design = Channel.fromPath(params.design).splitCsv(header : true).multiMap {
        row ->
        mat: [row.sample, "maternal", file(row.maternal_1), file(row.maternal_2)]
        pat: [row.sample, "paternal", file(row.paternal_1), file(row.paternal_2)]
        hifi: [row.sample, file(row.hifi)]
    }
}

yac_fastq_ch = design.mat.concat(design.pat)

if(params.hifi_parents) {
    process yak_kmers_single {
        cpus 40
        memory '170GB'
        time '12h'

        input:
        set val(sample), val(parent), file(fastq) from yac_fastq_ch

        output:
        set val(sample), file("${parent}.yak") into yac_kmers_ch

        script:
        """
        yak count -b37 -t40 -o ${parent}.yak ${fastq}
        """
    }
} else {
    process yak_kmers_paired {
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

if(params.dipcall) {
    process dipcall_variants {
        cpus 40
        memory '160GB'
        time '24h'
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
        make -j2 -f ${sample}.mak
        """
    }
} else{
    process svim_asm_variants {
        cpus 40
        memory '160GB'
        time '24h'
        publishDir 'variants'

        input:
        set val(sample), file(asm), file(ref) from hifiasm_denovo_ch.combine(ref_ch)

        output:
        set val(sample), file("${sample}.dip.vcf.gz"), file("${sample}.snvs.dip.vcf.gz")

        script:
        """
        module load samtools
        module load bcftools
        module load python/3.9.6

        samtools faidx ${ref}

        mkdir ${sample}
        minimap2 -a -x asm5 --cs -r2k -t 40 ${ref} ${asm}/${sample}_hap1.fa.gz | samtools sort -m4G -@4 -o hap1.sorted.bam -
        minimap2 -a -x asm5 --cs -r2k -t 40 ${ref} ${asm}/${sample}_hap2.fa.gz | samtools sort -m4G -@4 -o hap2.sorted.bam -

        samtools index hap1.sorted.bam
        samtools index hap2.sorted.bam

        mkdir snvs
        (seq 1 22; echo X) | parallel -j24 'python3 ${projectDir}/parse_snvs.py --min_quality 40 --reference ${ref} --hap1 hap1.sorted.bam --hap2 hap2.sorted.bam --region chr{} --vcf_out snvs/chr{}.vcf.gz --vcf_template ${projectDir}/header.vcf.gz --sample ${sample}'
        bcftools concat snvs/*.vcf.gz -Oz -o ${sample}.snvs.dip.vcf.gz

        svim-asm diploid --min_sv_size 1 --types DEL,INS --sample ${sample} ${sample}/ hap1.sorted.bam hap2.sorted.bam ${ref}
        mv ${sample}/variants.vcf ${sample}.dip.vcf
        bgzip ${sample}.dip.vcf
        """
    }

}
