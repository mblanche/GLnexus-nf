params.input = null
params.outdir = null
params.study_name = 'study'

workflow  {
  Channel.fromPath(params.input)
    .splitCsv(header:true)
    .map { row-> 
            tuple(row.sample, file(row.gVCF))
        }
    .set { gVCFs }
    // .map { _sample, gvcf -> gvcf}
    // .collect()
    // .set{ filePaths }

  getChromsomeInfo(gVCFs.take(1))
    .splitCsv()
    .set { chrSize }

  indexVCF = indexVCF(gVCFs)

  splitVCF = splitVCF(chrSize
    .map{ chr, _size -> chr}
    .combine(indexVCF.map{_id, vcf, csi -> [vcf,csi]}))
    .groupTuple()

  glnexus = glnexus(chrSize.join(splitVCF))

  bcf2vcf(glnexus.bcf)

  /**/
}

process  indexVCF {
  cpus 2
  memory '4 GB'
  container 'community.wave.seqera.io/library/bcftools:1.21--374767bf77752fc2'

  input:
  tuple val(id), path(vcf)
  
  output:
  tuple val(id), path(vcf), path('*.csi'), emit: vcf

  script:
  """
  bcftools index --threads ${task.cpus} ${vcf}
  """
}

process splitVCF {
  cpus 1
  memory '3 GB'
  container 'community.wave.seqera.io/library/bcftools:1.21--374767bf77752fc2'
  
  input:
  tuple val(chr), path(vcf), path(csi)

  output:
  tuple val(chr), path('*.vcf.gz')

  script:
  """
  bcftools view ${vcf} ${chr} \\
  | bgzip -@ ${task.cpus} -c > ${vcf.baseName}-${chr}.vcf.gz
  """
}

process getChromsomeInfo {
  cpus 1
  memory '2 GB'
  container 'community.wave.seqera.io/library/bcftools:1.21--374767bf77752fc2'

  input: 
  tuple val(id), path(vcf)

  output:
  path('chrInfo.txt'), emit: chrInfo

  script:
  """
  bcftools view -h ${vcf} \
  | awk -v OFS=',' '/chr([1-9]|1[0-9]|2[0-2]|[XY]),/ {
  gsub(/##contig=<|>/,"")
  split(\$0,a,",")
  split(a[1], chr, "=")
  split(a[2], len, "=")
  print chr[2], len[2]
  }' \
  | tee chrInfo.txt
  """
}

//TODO Need to code the memory size
process glnexus {
  cpus 16
  memory "72 GB"
  container 'community.wave.seqera.io/library/glnexus_jemalloc:f7d379f09441d9b8'

  input:
  tuple val(chr), val(size), path(vcfs)

  output:
  val(chr), emit: chr
  path('*.bcf'), emit: bcf

  script:
  """
  LD_PRELOAD=\$MAMBA_ROOT_PREFIX/lib/x86_64-linux-gnu/libjemalloc.so \\
  glnexus_cli -t ${task.cpus} -m 64 --config DeepVariant \\
  --bed <(echo -e '${chr}\t1\t${size}')  \\
  ${vcfs} > ${params.study_name}-${chr}.bcf
  """
}

process bcf2vcf {
  cpus 1
  memory "2 GB"
  container 'community.wave.seqera.io/library/bcftools:1.21--374767bf77752fc2'

  publishDir "${params.outdir}/joint-genotyping"

  input:
  path(bcf)

  output:
  path('*.vcf.gz')

  script:
  """
  bcftools view --threads ${task.cpus} ${bcf} \\
  | bgzip -@ ${task.cpus} -c > ${bcf.baseName}.vcf.gz
  """

}

