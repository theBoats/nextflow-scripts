#!/usr/bin/nextflow

params.samplesheet = "/fs03/ha66/graham/temp/truly_temp/sample-sheet.csv"


process PARSE {
  debug true
  
  input:
  tuple val(sampleID), val(read1), val(read2)

  output:
  tuple val(sampleID), val(read1), val(read2)

  script:
  """
  echo "+Processing " ${sampleID} 
  echo "+Read1: " ${read1} 
  echo "+Read2: " ${read2}
  """
}

process CONCATENATE_LANES {

  debug true

  publishDir 'sequence/concatenated_fastq', mode: 'link'


  input:
  tuple val(sampleID), val(read1), val(read2)

  output:
  tuple val(sampleID), path("${sampleID}_R1.fastq.gz"), path("${sampleID}_R2.fastq.gz")


  script:
  """
  echo Concatenating lanes on ${sampleID}
  cat ${read1} > ${sampleID}_R1.fastq.gz
  cat ${read2} > ${sampleID}_R2.fastq.gz
  """
}


// process COPY_RAW_DATA {
//   debug true

//   publishDir 'test_output', mode: 'copy'
  
//   input:
//   tuple val(sampleID), path(read1), path(read2)

//   output:
//   tuple val(sampleID), path(read1), path(read2)

//   script:
//   """
//   cp ${read1} "${read1}_raw"
//   mv "${read1}_raw" ${read1}

//   cp ${read2} "${read2}_raw"
//   mv "${read2}_raw" ${read2}
//   """
// }


process FASTQC {

  module 'fastqc/0.12.1'

  publishDir 'sequence/fastqc', mode: 'link'
  
  input:
  tuple val(sampleID), path(read1), path(read2)

  output:
  path("*")

  script:
  """
  echo Running fastqc on ${sampleID}
  fastqc ${read1}
  fastqc ${read2}
  """
}

process FASTQ_SCREEN {
  debug true

  publishDir 'sequence/fastq_screen', mode: 'link'

  module 'bowtie2/2.3.5'

  input:
  tuple val(sampleID), path(read1), path(read2)

  output:
  path("*")

  script:
  """
  echo +Running fastq screen on ${sampleID} 

  fastq_screen --conf /scratch/ha66/local/bin/fastq_screen_v0.12.0/fastq_screen.conf \
  --aligner bowtie2 \
  --threads 6 \
  ${read1} \
  ${read2}
   """
}


process FASTP {

  module 'fastp'

  publishDir 'sequence/fastp', mode: 'link'

  input:
  tuple val(sampleID), path(read1), path(read2)

  output:
  path("*.fastq.gz") , emit: fastp_files
  path("*.html")
  path("*.json")

  // tuple val(sampleID), path("${sampleID}.fastp.fastq.gz"), emit: fastp_files
  // path("${sampleID}.fastp.html")
  // path("${sampleID}.fastp.json")

  script:
  """
  echo +Running fastp on ${sampleID}

  fastp \
  -l 35 \
  -W 3 \
  -5 \
  -3 \
  --thread 6 \
  --json=${sampleID}.fastp.json \
  --html=${sampleID}.fastp.html \
  --report_title=${sampleID} \
  -i ${read1} \
  -o ${read1}_R1.fastp.fastq.gz \
  -I ${read2} \
  -O ${read2}_R2.fastp.fastq.gz
  """
}




workflow {
  Channel.fromPath(params.samplesheet) \
    | splitCsv(header:true) \
    | map { row -> tuple(row.sampleID, file(row.read1), file(row.read2)) } \
    | PARSE \
    | CONCATENATE_LANES \
    | (FASTQC & FASTQ_SCREEN & FASTP)
}