#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// multilane files
// Lanes is supplied as parameter; default 4, should catch errors where less than four lanes exist
// Only match Read 1
// Fastq files are supplied at folder level as a parameter

// params.lanes = 4


// params.fastqFolder = "/mnt/lustre/working/lab_kateg/grahamM/temp/nf-practice/dummy-data-multi-lane" // empty files
// params.fastqFolder = "/mnt/lustre/working/lab_kateg/grahamM/temp/nf-practice/data-multi-lane"


/*

FEATURES TO ADD:
Multi genome support
Consider parameterizing STAR threads

*/

// m3 test data
params.lanes = 4
params.fastqFolder = "/fs03/ha66/graham/temp/temp_fastq"
params.genome = "/fs03/ha66/references/STAR-index2/iGenomes/Mus_musculus/UCSC/mm9"

fastqFiles = params.fastqFolder + "/*L00${(1..params.lanes)}*_R1*.fastq.gz"

process COPY_RAW_DATA {

  // DOCSTRING: copies raw fastqfiles to output directory

  publishDir 'sequence/fastq_files', mode: 'copy'

  input:
  path(file)

  output:
  path(file)

  script:
  """
  cp ${file} "${file}_raw"
  mv "${file}_raw" ${file}
  """
}


process CONCAT_FASTQ {
  // concatenates multiple fastq files from separate lanes

  // files are published relative to script directory
  publishDir 'sequence/concatenated_fastq', mode: 'copy'

  input:
  tuple val(sampleID), path(files)

  output:
  tuple val(sampleID), path("${sampleID}.concatenated.fastq.gz")

  script:
  """
  cat ${files} > ${sampleID}.concatenated.fastq.gz
  """
}

process FASTQC {

  module 'fastqc/0.12.1'

  publishDir 'sequence/fastqc', mode: 'copy'
  
  input:
  tuple val(sampleID), path(concatenated_fastq)

  output:
  path("${sampleID}.concatenated_fastqc.html")
  path("${sampleID}.concatenated_fastqc.zip")

  script:
  """
  fastqc ${concatenated_fastq}
  """
}


// Variables for TRIM
polyA="/usr/local/bbmap/38.81/resources/polyA.fa.gz"
polyG="/fs03/ha66/graham/common_scripts/bbmap/polyG.fa"
truseqpolyA="/usr/local/bbmap/38.81/resources/truseq.fa.gz"
adapters="/usr/local/bbmap/38.81/resources/adapters.fa"
poolseq_oligos="/fs03/ha66/graham/common_scripts/bbmap/poolseq-oligos.fa"


process TRIM {

  // on m3 bbmap loaded as module. At QIMR local install.
  module 'bbmap/38.81'

  // files are published relative to script directory
  publishDir 'sequence/trimmed_fastq', mode: 'copy'

  input:
  tuple val(sampleID), path(fastqfile)

  output:
  tuple val(sampleID), path("${sampleID}.trimmed.fastq.gz")

  script:
  """
  bbduk.sh in=${fastqfile} \
  out=${sampleID}.trimmed.fastq \
  ref=${polyA},${adapters},${truseqpolyA},${polyG},${poolseq_oligos}

  gzip ${sampleID}.trimmed.fastq
  """
}


// STAR variables

process ALIGN {

  module 'star/2.5.2b'

  // !! fix this to genome variable in future
  publishDir 'aligned/mm9', mode: 'copy'

  input:
  tuple val(sampleID), path(fastqfile)

  output:
  path("${sampleID}_Aligned.sortedByCoord.out.bam")
  path("${sampleID}_Log.final.out")
  path("${sampleID}_Log.out")
  path("${sampleID}_Log.progress.out")
  path("${sampleID}_SJ.out.tab")
  path("${sampleID}_ReadsPerGene.out.tab"), emit: gene_counts 

  script:
  """
  STAR \
  --runThreadN 6 \
  --genomeDir ${params.genome} \
  --readFilesIn ${fastqfile} \
  --alignSJDBoverhangMin 1 \
  --outFilterMismatchNoverLmax 0.1 \
  --alignIntronMax 1000 \
  --outFileNamePrefix "${sampleID}_" \
  --outSAMtype BAM SortedByCoordinate \
  --readFilesCommand zcat \
  --quantMode GeneCounts
  """

}


process WRANGLE_COUNTS_PYTHON {

  // module 'python/3.8.5'
  // Source env outside of running nextflow source ~/gm/python_envs/py387_stat_env/bin/activate

  publishDir 'summary_counts', mode: 'copy'

  input:
  path(count_files)

  output:
  path("summary-counts-pool.csv")

  script:
  """
  #!/usr/bin/env python

  import pandas as pd
  import os

  files = "${count_files}".split(" ")

  # pull gene ids from first file
  gene_ids = []

  with open(files[0]) as f:
    lines = f.readlines()[4:]
    for line in lines:
      gene_ids.append(line.split('\t')[0])

  summary_dict = {'gene_id' : gene_ids}


  # put counts for each sample into dictionary
  for sample in files:
    
    sample_name = os.path.basename(sample).split('_')[0]
    
    counts = []
    with open(sample) as f:
      lines = f.readlines()[4:]
      for line in lines:
        counts.append(line.split('\t')[1])

    summary_dict[sample_name] = counts

  # create data frame from summary dictionary
  summary_df = pd.DataFrame(summary_dict)

  summary_df.to_csv("./summary-counts-pool.csv")
  """
}



workflow {
  
  // create copy of raw data in output directory
  raw_data_ch = Channel.fromPath(params.fastqFolder + "/*.fastq.gz")
  COPY_RAW_DATA(raw_data_ch)

  
  // MAIN
  fastq_ch = Channel.fromFilePairs(fastqFiles, size: params.lanes) 
  concat_fastq_ch = CONCAT_FASTQ(fastq_ch)

  // pass concatenated fastq into two streams
  FASTQC(concat_fastq_ch)
  TRIM(concat_fastq_ch) | ALIGN
  
  // Pass only gene counts from ALIGN
  ALIGN.out.gene_counts | collect | WRANGLE_COUNTS_PYTHON
}


workflow.onComplete {
  println "\nPipeline completed at: $workflow.complete"
  println "Execution status: ${ workflow.success ? 'OK' : 'failed' } \n"

}

workflow.onError {
  println "Pipeline error'd. Message: ${workflow.errorMessage}"
}




