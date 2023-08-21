/*
 * Trim paired-end reads using default options in trim-galore (CUTADAPT)
 * Scope: GLOBAL
 * Input: PE Fastq files with grouping ID
 * Emits: trimmed reads
 * Emits: cutadapt trimming report
 * Feeds: fastqc
 * Feeds: Mapping step (Kall, Bt2)
 */


process TRIM {
  tag "trimming $id"
  label 'med_mem'

  input:
  tuple val(id), path(reads)

  output:
  path("${id}*.txt"), emit: trim_report
  tuple val(id), path("${id}*.fq.gz"), emit: trimmed_reads

  script:
  """
  trim_galore --paired --basename ${id} -j $task.cpus ${reads[0]} ${reads[1]}
  """
}
