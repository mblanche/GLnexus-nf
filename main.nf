params.input = null
params.study_name = 'study'
// params.chrBed
// params.glnexusConfig


workflow  {
  Channel.fromPath(params.input)
    .splitCsv(header:true)
    .map { row-> 
            tuple(row.sample, file(row.gVCF))
        }
    .tap { gVCFs }
    .map { sample, gvcf -> gvcf}
    .take(2)
    .collect()
    .set{ filePaths }

  getChromsomeInfo(gVCFs.take(1))
    .splitCsv()
    .map{ chr, size ->
    if (chr == 'chr22') {
      tuple(chr,size)
      }
    }
    .set { chrSize }
    
  glnexus(chrSize, filePaths)
}


process getChromsomeInfo {
  container 'community.wave.seqera.io/library/bcftools:1.21--374767bf77752fc2'
  //container '096672585100.dkr.ecr.us-west-1.amazonaws.com/bcftools:latest'

  input: 
  tuple val(id), path(vcf)

  output:
  path('chrInfo.txt')

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
  memory "16G"
  //container 'ghcr.io/dnanexus-rnd/glnexus:v1.4.1'
  container 'community.wave.seqera.io/library/glnexus_jemalloc:f7d379f09441d9b8'

  input:
  tuple val(chr), val(size)
  path(vcfs)

  output:

  script:
  """
  glnexus_cli -t ${task.cpus} --config DeepVariant \\
  --bed <(echo -e '${chr}\t1\t${size}')  \\
  ${vcfs} > ${params.study_name}-${chr}.bcf
  """
}

