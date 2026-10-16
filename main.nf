params.design = null
params.ref = null
params.assemblies = false
params.assemblies_outdir = 'assemblies'
params.dipcall = false
params.svim_asm = false
params.cpus = 40
params.memory = '160G'


process bam_to_fastq {
  cpus params.cpus
  memory params.memory
  time '12h'

  input:
  tuple val(sample), val(family), val(role), path(bam)

  output:
  tuple val(sample), val(family), val(role), path("${sample}.fastq.gz")

  script:
  """
  samtools fastq -@ ${task.cpus} -n ${bam} | gzip > ${sample}.fastq.gz
  """

  stub:
  """
  touch ${sample}.fastq.gz
  """
}

process yak_kmers_single {
  cpus params.cpus
  memory params.memory
  time '12h'

  input:
  tuple val(sample), val(family), val(role), path(fastq)

  output:
  tuple val(family), val(role), path("${sample}.yak")

  script:
  """
  yak count -b37 -t${task.cpus} -o ${sample}.yak ${fastq}
  """

  stub:
  """
  touch ${sample}.yak
  """
}

process hifiasm_trio_denovo {
  cpus params.cpus
  memory params.memory
  time '24h'

  input:
  tuple val(sample), path(hifi_fasta), path(pat_yak), path(mat_yak)

  output:
  tuple val(sample), path("${sample}_diploid/${sample}_hap1.fa.gz"), path("${sample}_diploid/${sample}_hap2.fa.gz")

  script:
  """
  mkdir ${sample}_asm
  hifiasm -o ${sample}_asm/${sample}.asm -t${task.cpus} -1 ${pat_yak} -2 ${mat_yak} --dual-scaf ${hifi_fasta}

  mkdir ${sample}_diploid
  awk '/^S/{print ">"\$2;print \$3}' ${sample}_asm/${sample}.asm.dip.hap1.p_ctg.gfa | bgzip > ${sample}_diploid/${sample}_hap1.fa.gz
  awk '/^S/{print ">"\$2;print \$3}' ${sample}_asm/${sample}.asm.dip.hap2.p_ctg.gfa | bgzip > ${sample}_diploid/${sample}_hap2.fa.gz

  rm -r ${sample}_asm
  """

  stub:
  """
  mkdir ${sample}_diploid
  touch ${sample}_diploid/${sample}_hap1.fa.gz ${sample}_diploid/${sample}_hap2.fa.gz
  """
}

process hifiasm_hifi_denovo {
  cpus params.cpus
  memory params.memory
  time '24h'

  input:
  tuple val(sample), path(hifi_fasta)

  output:
  tuple val(sample), path("${sample}_diploid/${sample}_hap1.fa.gz"), path("${sample}_diploid/${sample}_hap2.fa.gz")

  script:
  """
  mkdir ${sample}_asm
  hifiasm -o ${sample}_asm/${sample}.asm -t${task.cpus} --dual-scaf ${hifi_fasta}

  mkdir ${sample}_diploid
  awk '/^S/{print ">"\$2;print \$3}' ${sample}_asm/${sample}.asm.bp.hap1.p_ctg.gfa | bgzip > ${sample}_diploid/${sample}_hap1.fa.gz
  awk '/^S/{print ">"\$2;print \$3}' ${sample}_asm/${sample}.asm.bp.hap2.p_ctg.gfa | bgzip > ${sample}_diploid/${sample}_hap2.fa.gz

  rm -r ${sample}_asm
  """

  stub:
  """
  mkdir ${sample}_diploid
  touch ${sample}_diploid/${sample}_hap1.fa.gz ${sample}_diploid/${sample}_hap2.fa.gz
  """
}

process rename_contigs_pansn {
  cpus 1
  memory '1G'
  time '2h'

  publishDir params.assemblies_outdir, mode: 'copy'

  input:
  tuple val(sample), path(hap1), path(hap2)

  output:
  tuple val(sample), path("${sample}_diploid/${sample}_hap1.fa.gz"), path("${sample}_diploid/${sample}_hap2.fa.gz")

  script:
  """
  mkdir ${sample}_diploid
  gzip -cdf ${hap1} | awk -v prefix='${sample}#1#' '/^>/ { if (index(\$1, "#")) { print "Contig identifiers must not contain # before PanSN renaming" > "/dev/stderr"; exit 1 } sub(/^>/, ">" prefix) } { print }' | bgzip > ${sample}_diploid/${sample}_hap1.fa.gz
  gzip -cdf ${hap2} | awk -v prefix='${sample}#2#' '/^>/ { if (index(\$1, "#")) { print "Contig identifiers must not contain # before PanSN renaming" > "/dev/stderr"; exit 1 } sub(/^>/, ">" prefix) } { print }' | bgzip > ${sample}_diploid/${sample}_hap2.fa.gz
  """

  stub:
  """
  mkdir ${sample}_diploid
  touch ${sample}_diploid/${sample}_hap1.fa.gz ${sample}_diploid/${sample}_hap2.fa.gz
  """
}

process dipcall_variants {
  cpus params.cpus
  memory params.memory
  time '24h'
  publishDir 'variants', mode: 'copy'

  input:
  tuple val(sample), path(hap1), path(hap2), path(ref)

  output:
  tuple val(sample), path("${sample}.dip.vcf.gz")

  script:
  """
  samtools faidx ${ref}

  /root/dipcall/run-dipcall ${sample} ${ref} ${hap1} ${hap2} > ${sample}.mak
  make -j2 -f ${sample}.mak
  """

  stub:
  """
  touch ${sample}.dip.vcf.gz
  """
}

process svim_asm_variants {
  cpus params.cpus
  memory params.memory
  time '24h'
  publishDir 'variants', mode: 'copy'

  input:
  tuple val(sample), path(hap1), path(hap2), path(ref)

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
  (seq 1 22; echo X) | parallel -j24 'parse_snvs.py --min_quality ${params.cpus} --reference ${ref} --hap1 hap1.sorted.bam --hap2 hap2.sorted.bam --region chr{} --vcf_out snvs/chr{}.vcf.gz --vcf_template ${projectDir}/header.vcf.gz --sample ${sample}'
  bcftools concat snvs/*.vcf.gz -Oz -o ${sample}.snvs.dip.vcf.gz

  svim-asm diploid --min_sv_size 1 --types DEL,INS --sample ${sample} ${sample}/ hap1.sorted.bam hap2.sorted.bam ${ref}
  mv ${sample}/variants.vcf ${sample}.dip.vcf
  bgzip ${sample}.dip.vcf
  """

  stub:
  """
  touch ${sample}.dip.vcf.gz ${sample}.snvs.dip.vcf.gz
  """
}

workflow {
  if (!params.assemblies && !params.design) {
    error 'Provide --design with a BAM,sample_id,sample_family,role CSV, or --assemblies with existing haplotypes.'
  }
  if ((params.dipcall || params.svim_asm) && !params.ref) {
    error 'Variant calling requires --ref.'
  }
  if(!params.assemblies) {
    design = channel.fromPath(params.design, checkIfExists: true)
      .splitCsv(header: true)
      .toList()
      .flatMap { rows ->
        if (!rows) {
          error 'The design CSV must contain at least one sample.'
        }
        def required = ['BAM', 'sample_id', 'sample_family', 'role']
        def samples = rows.collect { row ->
          required.each { column ->
            if (!row[column]?.trim()) {
              error "Each design row requires a non-empty '${column}' column."
            }
          }
          def sample = row.sample_id.trim()
          def family = row.sample_family.trim()
          def role = row.role.trim()
          [sample, family].each { name ->
            if (!(name ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) {
              error "Invalid sample or family identifier '${name}': use letters, numbers, periods, underscores, and hyphens, starting with a letter or number."
            }
          }
          if (!(role in ['mat', 'pat', 'proband'])) {
            error "Invalid role '${role}' for '${sample}': expected mat, pat, or proband."
          }
          [sample, family, role, file(row.BAM.trim(), checkIfExists: true)]
        }
        samples.groupBy { it[0] }.each { sample, members ->
          if (members.size() != 1) {
            error "Duplicate sample_id '${sample}': provide one BAM per sample."
          }
        }
        samples.groupBy { it[1] }.each { family, members ->
          ['mat', 'pat'].each { role ->
            def count = members.count { it[2] == role }
            if (count > 1 || (members.any { it[2] == 'proband' } && count != 1)) {
              error "Family '${family}' requires exactly one '${role}' parent when a proband is present, and cannot contain duplicate parents."
            }
          }
        }
        return samples
      }

    reads = bam_to_fastq(design).branch {
      parents: it[2] in ['mat', 'pat']
      probands: it[2] == 'proband'
    }

    parent_assemblies = hifiasm_hifi_denovo(reads.parents.map {
      sample, family, role, fastq -> [sample, fastq]
    })
    parent_kmers = yak_kmers_single(reads.parents).branch {
      maternal: it[1] == 'mat'
      paternal: it[1] == 'pat'
    }
    parental_pairs = parent_kmers.paternal.map { family, role, yak -> [family, yak] }
      .join(parent_kmers.maternal.map { family, role, yak -> [family, yak] })

    trio_inputs = reads.probands.map {
      sample, family, role, fastq -> [family, sample, fastq]
    }.combine(parental_pairs, by: 0).map {
      family, sample, fastq, pat_yak, mat_yak -> [sample, fastq, pat_yak, mat_yak]
    }
    proband_assemblies = hifiasm_trio_denovo(trio_inputs)
    hifiasm_denovo_ch = parent_assemblies.mix(proband_assemblies)
  } else {
    hifiasm_denovo_ch = channel.fromPath(params.assemblies, checkIfExists: true).splitCsv(header: true).map {
      row ->
      if (!(row.sample ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) {
        error "Invalid assembly sample identifier '${row.sample}': use letters, numbers, periods, underscores, and hyphens, starting with a letter or number."
      }
      [row.sample, file(row.hap1, checkIfExists: true), file(row.hap2, checkIfExists: true)]
    }
  }

  pansn_assemblies = rename_contigs_pansn(hifiasm_denovo_ch)

  if(params.dipcall) {
    ref_ch = channel.fromPath(params.ref, checkIfExists: true)
    dipcall_variants(pansn_assemblies.combine(ref_ch))
  } else if(params.svim_asm) {
    ref_ch = channel.fromPath(params.ref, checkIfExists: true)
    svim_asm_variants(pansn_assemblies.combine(ref_ch))
  }
}
