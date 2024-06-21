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
Fix genome publishDir relative to a param
Consider parameterizing STAR threads

*/

// m3 test data
params.lanes = 4
params.fastqFolder = "/fs03/ha66/graham/temp/temp_fastq"
params.genome = "/fs03/ha66/references/STAR-index2/iGenomes/Mus_musculus/UCSC/mm9"

// Need to explicitly code working and publish directories on PBS !! NOT SURE WHY
// NOT FIXED. STILL WRITING TO HOME FOLDER

// global params
fastqFiles = params.fastqFolder + "/**/*L00${(1..params.lanes)}*_R1*.fastq.gz"
params.help = false

// limit input files for testing
params.dev = false
params.number_of_inputs = 2



if (params.help) {
  log.info """

  Nextflow requires java version >= 11 and <= 18. Check available modules: `module avail java`
  This script requires a version of python with pandas installed.
  Output is written relative to directory script is executed in.

  On QIMR HPC:
  Load interactive session.

  Load the following modules:
  `module load java/17.0.4.1`
  `module load python/3.9.13`

  Make sure nextflow is in your PATH or execute from:
  `/mnt/lustre/working/lab_kateg/grahamM/software/nextflow`


  Usage:
    nextflow run pool-seq-nextflow.nf --lanes 4 --genome 'STARindex/../mm9' --fastqFolder "/fs03/ha66/graham/temp/temp_fastq"

  Input:
    * --lanes: Number of lanes run on Illumina sequencer, between 1 and 4
    * --genome: Path to STAR genome index.
    * --fastqFolder: Path to a directory containing all of the fastqFiles
    * --help: Show this message


  """
  exit 0
}



// PROCESSES


process COPY_RAW_DATA {

  // DOCSTRING: copies raw fastqfiles to output directory

  label 'small'

  publishDir "sequence/fastq_files", mode: "link", overwrite: true

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

  label 'small'

  // files are published relative to script directory
  publishDir "sequence/concatenated_fastq", mode: "link", overwrite: true

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

  label 'small'

  module 'fastqc/0.12.1'

  publishDir "sequence/fastqc", mode: "link", overwrite: true
  
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


process FASTP {

  label 'small'

  module 'fastp/0.23.2'

  publishDir "sequence/fastp", mode: "link", overwrite: true

  input:
  tuple val(sampleID), path(fastqfile)

  output:
  tuple val(sampleID), path("${sampleID}.fastp.fastq.gz"), emit: fastp_files
  path("${sampleID}.fastp.html")
  path("${sampleID}.fastp.json")

  script:
  """
  fastp --json=${sampleID}.fastp.json --html=${sampleID}.fastp.html --report_title=${sampleID} -i ${fastqfile} -o ${sampleID}.fastp.fastq.gz
  """


}


// STAR variables

process ALIGN {

  label 'star_align'


  // try with default module
  module 'STAR'

  // !! fix this to genome variable in future
  publishDir "aligned/mm10", mode: "link", overwrite: true

  input:
  tuple val(sampleID), path(fastqfile)

  output:
  tuple val(sampleID), path("${sampleID}_Aligned.sortedByCoord.out.bam"), emit: bam_files
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

  label 'small'

  module 'python'

  publishDir "summary_counts", mode: "link", overwrite: true

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

process SORT_FOR_PRESEQ {

  label 'small'

  module 'bedtools'

  publishDir "bed_files", mode: 'link'

  input:
  tuple val(sampleID), path(bamfile)

  output:
  tuple val(sampleID), path("${sampleID}.sorted.bed"), emit: sorted_bed

  script:
  """
  bedtools bamtobed -i ${bamfile} > ${sampleID}.bed
  sort -k 1,1 -k 2,2n -k 3,3n -k 6,6 ${sampleID}.bed > ${sampleID}.sorted.bed
  """
}

process PRESEQ_C_CURVE {

  /*
   *  Generate the complexity curve
   */

  label 'small'

  publishDir 'preseq', mode: 'link'

  input:
  tuple val(sampleID), path(bedFile)

  output:
  tuple val(sampleID), path("${sampleID}_output.txt"), emit: c_curve_output

  script:
  """
  preseq c_curve -o ${sampleID}_output.txt ${bedFile}
  """
}

process PRESEQ_LC_EXTRAP {

  label 'small'

  publishDir 'preseq', mode: 'link'

  input:
  tuple val(sampleID), path(bedFile)

  output:
  path("*")

  script:
  """
  preseq lc_extrap -o ${sampleID}_future_yield.txt ${bedFile}
  """
}

process PRESEQ_BOUND_POP {

  label 'small'

  publishDir 'preseq', mode: 'link'

  input:
  tuple val(sampleID), path(bedFile)

  output:
  path("*")

  script:
  """
  preseq bound_pop -o ${sampleID}_species_richness.txt ${bedFile}
  """
}

process PLOT_C_CURVE {

  label 'small'

  publishDir "preseq/plots", mode: 'link'

  module 'python/3.9.13'

  input:
  tuple val(sampleID), path(c_curve_output)

  output:
  path("*")

  script:
  """
  #!/usr/bin/env python

  import matplotlib.pyplot as plt 
  import pandas as pd

  df = pd.read_csv("${c_curve_output}", header = 'infer', delim_whitespace = True)

  plt.plot(df["total_reads"], df["total_reads"], label="hypothetical")
  plt.plot(df["total_reads"], df["distinct_reads"], label="${sampleID}")
  plt.xlabel("Total reads (M)")
  plt.ylabel("Distinct reads (M)")
  plt.legend(loc="upper left")
  plt.title("Complexity analysis")
  plt.savefig("${sampleID}_c_curve.png")
  """
}

// process MULTIQC {

//   debug true

//   label 'small'

//   publishDir 'sequence/multiqc'

//   module 'singularity/3.7.1'

//   input:
//   path(output_directory)

//   output:
//   path('*')

//   script:
//   """
//   singularity run /mnt/lustre/working/lab_kateg/grahamM/software/singularity/multiqc-1.20.sif ${output_directory}
//   """

// }



workflow {
  
  // create copy of raw data in output directory
  raw_data_ch = Channel.fromPath(params.fastqFolder + "/**/*.fastq.gz")
  .take( params.dev ? params.number_of_inputs : -1 )
  
  COPY_RAW_DATA(raw_data_ch)

  
  // MAIN
  fastq_ch = Channel.fromFilePairs(fastqFiles, size: params.lanes) 
   .take( params.dev ? params.number_of_inputs : -1 )


  concat_fastq_ch = CONCAT_FASTQ(fastq_ch)

  // pass concatenated fastq into two streams
  FASTQC(concat_fastq_ch)
  // TRIM(concat_fastq_ch) | ALIGN
  
  FASTP(concat_fastq_ch) 
  FASTP.out.fastp_files | ALIGN
  
  // Pass only gene counts from ALIGN
  ALIGN.out.gene_counts | collect | WRANGLE_COUNTS_PYTHON

  ALIGN.out.bam_files | SORT_FOR_PRESEQ | (PRESEQ_C_CURVE & PRESEQ_LC_EXTRAP & PRESEQ_BOUND_POP)

  PRESEQ_C_CURVE.out.c_curve_output | PLOT_C_CURVE

  // // MultiQC
  // Channel.fromPath(params.outDir) { it.getParent() } | MULTIQC
}


workflow.onComplete {
  println "\nPipeline completed at: $workflow.complete"
  println "Execution status: ${ workflow.success ? 'OK' : 'failed' } \n"

}

workflow.onError {
  println "Pipeline error'd. Message: ${workflow.errorMessage}"
}





