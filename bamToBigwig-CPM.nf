
/*

Nextflow pipeline to create normalized bigwig files for RNA-seq

Note: publishDir results path is relative to the directory of this script

*/

nextflow.enable.dsl=2

params.bamFiles = "/fs03/ha66/graham/pipeline-mapped/Pool-seq/2023-03-01-G1E/aligned/mm9/*.bam"
params.name = "merged"




/*

CPU and memory management. This allows you to request the amount of CPU from the command line as a percentage. E.g. to use half the available CPUs
nextflow run main.nf --requestedCPU 50

You can do the same for memory. CURRENTLY DISABLED

*/

params.requestedCPU = 80
params.requestedMem = 50

maxcpus = Runtime.runtime.availableProcessors()
requestedCpus = Math.round((params.requestedCPU * maxcpus) / 100)

// maxmem = Runtime.runtime.totalMemory() / 10241024
// requestedMem = (Math.round((params.requestedMem * maxmem) / 100))

process mergeBAMfiles {

  label 'samtools'

  module 'samtools'

  cpus = requestedCpus

  publishDir "${baseDir}/results", mode: 'copy'

  input:
  tuple val(name), path(bamFiles)

  output:
  tuple val(name), path("${name}_merged.bam")

  script:
  """
  samtools merge ${name}_merged.bam \
  ${bamFiles}
  """
}


process bamCoverage {

  label 'deeptools'

  module 'deeptools/3.1.3'
  module 'samtools'

  cpus = requestedCpus

  publishDir "${baseDir}/results", mode: 'copy'

  input:
  tuple val(name), path(merged)

  output:
  path("${name}.bw")

  script:
  """
  samtools index ${merged}

  bamCoverage \
  -b ${merged} \
  -o ${name}.bw \
  --binSize 10 \
  -p ${task.cpus}  \
  --normalizeUsing CPM
  """
}

workflow {
  //Channel.fromPath(params.bamFiles) \
  Channel.fromFilePairs('/fs03/ha66/graham/pipeline-mapped/Pool-seq/2023-03-01-G1E/aligned/mm9/*{1,2,3,4,5,6}_Aligned.sortedByCoord.out.bam', size:6)
  | mergeBAMfiles
  | bamCoverage
}



