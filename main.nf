params.design = null
params.ref = null
params.hifi_parents = false
params.trio = false
params.assemblies = false
params.dipcall = false
params.svim_asm = false
params.cpus = 40
params.memory = '160G'


process yak_kmers_single {
  cpus params.cpus
  memory params.memory
  time '12h'

  input:
  tuple val(sample), val(parent), path(fastq)

  output:
  tuple val(sample), path("${parent}.yak")

  script:
  """
  yak count -b37 -t${params.cpus} -o ${parent}.yak ${fastq}
  """
}

process yak_kmers_paired {
  cpus params.cpus
  memory params.memory
  time '12h'

  input:
  tuple val(sample), val(parent), path(fastq_1), path(fastq_2)

  output:
  tuple val(sample), path("${parent}.yak")

  script:
  """
  yak count -b37 -t${params.cpus} -o ${parent}.yak <(cat ${fastq_1} ${fastq_2}) <(cat ${fastq_1} ${fastq_2})
  """
}

process hifiasm_trio_denovo {
  cpus params.cpus
  memory params.memory
  time '24h'

  publishDir 'assemblies'

  input:
  tuple val(sample), path(hifi_fasta), path(yak1), path(yak2)

  output:
  tuple val(sample), path("${sample}_diploid/${sample}_hap1.fa.gz"), path("${sample}_diploid/${sample}_hap2.fa.gz")

  script:
  """
  mkdir ${sample}_asm
  hifiasm -o ${sample}_asm/${sample}.asm -t${params.cpus} -1 ${yak1} -2 ${yak2} ${hifi_fasta}

  mkdir ${sample}_diploid
  awk '/^S/{print ">"\$2;print \$3}' ${sample}_asm/${sample}.asm.dip.hap1.p_ctg.gfa | bgzip > ${sample}_diploid/${sample}_hap1.fa.gz
  awk '/^S/{print ">"\$2;print \$3}' ${sample}_asm/${sample}.asm.dip.hap2.p_ctg.gfa | bgzip > ${sample}_diploid/${sample}_hap2.fa.gz

  rm -r ${sample}_asm
  """
}

process hifiasm_hifi_denovo {
  cpus params.cpus
  memory params.memory
  time '24h'

  publishDir 'assemblies'

  input:
  tuple val(sample), path(hifi_fasta)

  output:
  tuple val(sample), path("${sample}_diploid/${sample}_hap1.fa.gz"), path("${sample}_diploid/${sample}_hap2.fa.gz")

  script:
  """
  mkdir ${sample}_asm
  hifiasm -o ${sample}_asm/${sample}.asm -t${params.cpus} ${hifi_fasta}

  mkdir ${sample}_diploid
  awk '/^S/{print ">"\$2;print \$3}' ${sample}_asm/${sample}.asm.bp.hap1.p_ctg.gfa | bgzip > ${sample}_diploid/${sample}_hap1.fa.gz
  awk '/^S/{print ">"\$2;print \$3}' ${sample}_asm/${sample}.asm.bp.hap2.p_ctg.gfa | bgzip > ${sample}_diploid/${sample}_hap2.fa.gz

  rm -r ${sample}_asm
  """
}

process dipcall_variants {
  cpus params.cpus
  memory params.memory
  time '24h'
  publishDir 'variants'

  input:
  tuple val(sample), path(hap1), path(hap2), path(ref)

  output:
  tuple val(sample), path("${sample}.dip.vcf.gz")

  script:
  """
  samtools faidx ${ref}

  ~/dipcall.kit/run-dipcall ${sample} ${ref} ${asm}/${sample}_hap1.fa.gz ${asm}/${sample}_hap2.fa.gz > ${sample}.mak
  make -j2 -f ${sample}.mak
  """
}

process svim_asm_variants {
  cpus params.cpus
  memory params.memory
  time '24h'
  publishDir 'variants'

  input:
  tuple val(sample), path(hap1), path(hap2), path(ref) from hifiasm_denovo_ch.combine(ref_ch)

  output:
  tuple val(sample), path("${sample}.dip.vcf.gz"), path("${sample}.snvs.dip.vcf.gz")

  script:
  """
  parallel -j3 'samtools faidx {}' ::: ${ref} ${hap1} ${hap2}

  mkdir ${sample}
  minimap2 -a -x asm5 --cs -r2k -t ${params.cpus} ${ref} ${hap1} | samtools sort -m4G -@4 -o hap1.sorted.bam -
  minimap2 -a -x asm5 --cs -r2k -t ${params.cpus} ${ref} ${hap2} | samtools sort -m4G -@4 -o hap2.sorted.bam -

  parallel -j2 'samtools index {}' ::: hap1.sorted.bam hap2.sorted.bam

  mkdir snvs
  (seq 1 22; echo X) | parallel -j24 'python3 ${projectDir}/parse_snvs.py --min_quality ${params.cpus} --reference ${ref} --hap1 hap1.sorted.bam --hap2 hap2.sorted.bam --region chr{} --vcf_out snvs/chr{}.vcf.gz --vcf_template ${projectDir}/header.vcf.gz --sample ${sample}'
  bcftools concat snvs/*.vcf.gz -Oz -o ${sample}.snvs.dip.vcf.gz

  svim-asm diploid --min_sv_size 1 --types DEL,INS --sample ${sample} ${sample}/ hap1.sorted.bam hap2.sorted.bam ${ref}
  mv ${sample}/variants.vcf ${sample}.dip.vcf
  bgzip ${sample}.dip.vcf
  """
}

workflow {
  if(!params.assemblies) {
    if(params.trio && params.hifi_parents) {
      design = Channel.fromPath(params.design).splitCsv(header : true).multiMap {
        row ->
        mat: [row.sample, "maternal", file(row.maternal)]
        pat: [row.sample, "paternal", file(row.paternal)]
        hifi: [row.sample, file(row.hifi)]
      }
    } else if(params.trio && !params.hifi_parents) {
      design = Channel.fromPath(params.design).splitCsv(header : true).multiMap {
        row ->
        mat: [row.sample, "maternal", file(row.maternal_1), file(row.maternal_2)]
        pat: [row.sample, "paternal", file(row.paternal_1), file(row.paternal_2)]
        hifi: [row.sample, file(row.hifi)]
      }
    } else if(!params.trio) {
        design = Channel.fromPath(params.design).splitCsv(header : true).map {
          row -> [row.sample, file(row.hifi)]
        }
      }

      if(params.trio) {
        yac_fastq_ch = design.mat.concat(design.pat)

        yac_kmers_ch = channel.empty()
        if(params.hifi_parents) {
          yac_kmers_ch = yak_kmers_single(yac_fastq_ch)
        } else {
          yac_kmers_ch = yak_kmers_paired(yac_fastq_ch)
        }

        hifiasm_ch = design.hifi.join(yac_kmers_ch.groupTuple(by: 0).map{
          it -> [it[0], it[1][0], it[1][1]]}, by: 0).view()

        hifiasm_denovo_ch = hifiasm_trio_denovo(hifiasm_ch)
      } else if(!params.trio) {
        hifiasm_denovo_ch = hifiasm_hifi_denovo(design)
      }
    } else if(params.assemblies){
      hifiasm_denovo_ch = Channel.fromPath(params.assemblies).splitCsv(header : true).map {
        row -> [row.sample, file(row.hap1), file(row.hap2)]
      }
    }

    if(params.dipcall) {
      ref_ch = Channel.fromPath(params.ref)
      dipcall_variants(hifiasm_denovo_ch.combine(ref_ch))
    } else if(params.svim_asm) {
      ref_ch = Channel.fromPath(params.ref)
      svim_asm_variants(hifiasm_denovo_ch.combine(ref_ch))
    }
}
