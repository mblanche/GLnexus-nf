params.input = null
params.outDir = null
params.study_name = 'study'

workflow  {
  Channel.fromPath(params.input)
    .splitCsv(header:true)
    .map { row-> 
            tuple(row.sample, file(row.gVCF))
        }
    .tap { gVCFs }
    .map { sample, gvcf -> gvcf}
    .collect()
    .set{ filePaths }

  getChromsomeInfo(gVCFs.take(1))
    .splitCsv()
    .map{ chr, size ->   //TODO remove after testing
      if (chr == 'chr22') {
        tuple(chr,size)
        }
    }
    .set { chrSize }
    
  glnexus(chrSize, filePaths)
  bcf2vcf(glnexus.out)
  
}

process getChromsomeInfo {
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

process glnexus {
  cpus 16
  memory "252 GB" //fdasfas
  container 'community.wave.seqera.io/library/glnexus_jemalloc:f7d379f09441d9b8'

  input:
  tuple val(chr), val(size)
  path(vcfs)

  output:
  path('*.bcf'), emit: bcf

  script:
  """
  LD_PRELOAD=\$MAMBA_ROOT_PREFIX/lib/x86_64-linux-gnu/libjemalloc.so \\
  glnexus_cli -t ${task.cpus} --config DeepVariant \\
  --bed <(echo -e '${chr}\t1\t${size}')  \\
  ${vcfs} > ${params.study_name}-${chr}.bcf
  """

  stub:
  """
  touch ${params.study_name}-${chr}.bcf
  """
}

process bcf2vcf {
  cpus 4
  memory "16 GB"
  container 'community.wave.seqera.io/library/bcftools:1.21--374767bf77752fc2'

  publishDir "${params.outDir}/joint-genotyping"

  input:
  path(bcf)

  output:
  path('*.vcf.gz')

  script:
  """
  bcftools view --threads ${task.cpus} ${bcf} \\
  | bgzip -@ ${task.cpus} -c > ${bcf.baseName}.vcf.gz
  """

  stub:
  """
  touch ${params.study_name}-${chr}.vcf
  """
}

