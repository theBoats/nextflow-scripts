

/*

Deeptools plotting pipeline using DSL2 syntax

*/

nextflow.enable.dsl=2


params.bamFiles = "/home/gmaj0001/gm/pipeline-mapped/ChIP/K1/E2f4/mm9/samfilter/*.bam"
//params.bamFiles = (["/fs03/qk33/Neville/G-string_Chips/K119ub_DMS/mm10/samfilter/ChIP_K119UB__F3_DMSO.filt.bam", "/fs03/qk33/Neville/G-string_Chips/K119ub_SGC/mm10/samfilter/ChIP_K119UB__F3_SGC.filt.bam", "/fs03/qk33/Neville/G-string_Chips/K119ub_VTP/mm10/samfilter/ChIP_K119UB__F3_VTP.filt.bam"])
// params.outdir = "results"


/*

CPU and memory management. This allows you to request the amount of CPU from the command line as a percentage. E.g. to use half the available CPUs
nextflow run main.nf --requestedCPU 50

*/

params.requestedCPU = 80
params.requestedMem = 50

maxcpus = Runtime.runtime.availableProcessors()
requestedCpus = Math.round((params.requestedCPU * maxcpus) / 100)

// maxmem = Runtime.runtime.totalMemory() / 10241024
// requestedMem = (Math.round((params.requestedMem * maxmem) / 100))


process bamCoverage_without_normalization {

  label 'deeptools'

  module 'deeptools/3.1.3'
  module 'samtools'

  cpus = requestedCpus

  publishDir "${baseDir}/results", mode: 'copy'

  input:
  path bamFile 

  output:
  tuple val("no_norm"), path("${bamFile.baseName}.no.normalization.bw"), emit: norm_ch

  script:
  """
  samtools index ${bamFile}

  bamCoverage \
  -b $bamFile \
  -o ${bamFile.baseName}.no.normalization.bw \
  --binSize 10 \
  -p ${task.cpus}  \
  --normalizeUsing None
  """
}

process bamCoverage_with_CPM_normalization {

  label 'deeptools'

  module 'deeptools/3.1.3'
  module 'samtools'

  cpus = requestedCpus

  publishDir "${baseDir}/results", mode: 'copy'

  input:
  path bamFile 

  output:
  tuple val('CPM'), path("${bamFile.baseName}.CPM.normalization.bw"), emit: cpm_ch

  script:
  """
  samtools index ${bamFile}

  bamCoverage \
  -b $bamFile \
  -o ${bamFile.baseName}.CPM.normalization.bw \
  --binSize 10 \
  -p ${task.cpus}  \
  --normalizeUsing CPM
  """
}

process bamCoverage_with_RPGC_normalization {

  label 'deeptools'

  module 'deeptools/3.1.3'
  module 'samtools'

  cpus = requestedCpus

  publishDir "${baseDir}/results", mode: 'copy'

  input:
  path bamFile 

  output:
  tuple val('RPGC'), path("${bamFile.baseName}.RPGC.normalization.bw"), emit: rpgc_ch

  script:
  """
  samtools index ${bamFile}

  bamCoverage \
  -b $bamFile \
  -o ${bamFile.baseName}.RPGC.normalization.bw \
  --binSize 10 \
  -p ${task.cpus} \
  --normalizeUsing RPGC \
  --effectiveGenomeSize 2150570000
  """
}


process computeMatrix {

  label 'computeMatrix'

  module 'deeptools/3.1.3'

  cpus = requestedCpus

  publishDir "${baseDir}/results", mode: 'copy'
  
  input:
  tuple val(name), path(bigwigFiles)

  output:
  path "${name}.matrix.mat.gz"

  script:
  """
  computeMatrix scale-regions \
  -S ${bigwigFiles} \
  -R /home/gmaj0001/gm/Analysis/Deeptools/bed/mm10_UCSC_knownGene.bed \
  --beforeRegionStartLength 3000 \
  --regionBodyLength 5000 \
  --afterRegionStartLength 3000 \
  --skipZeros \
  -p max/2 \
  -o "${name}.matrix.mat.gz"
  """
}

process plotProfile {

    label 'deeptools'

    module 'deeptools/3.1.3'

    publishDir "${baseDir}/results", mode: 'copy'

    input:
    path matrixFile

    output:
    path "plotProfile-K119Ub-${matrixFile.baseName}-normalized.png"

    script:
    """
    plotProfile \
    -m ${matrixFile} \
    --perGroup \
    -out plotProfile-K119Ub-${matrixFile.baseName}-normalized.png \
    --plotTitle "plot-${matrixFile.baseName}"
    """
}


workflow {

  bam_ch = Channel.fromPath(params.bamFiles)

  (bamCoverage_without_normalization(bam_ch) | groupTuple).mix(bamCoverage_with_CPM_normalization(bam_ch) | groupTuple, bamCoverage_with_RPGC_normalization(bam_ch) | groupTuple) \
  | computeMatrix \
  | plotProfile
  
  

 