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
  tuple val(id), path("${id}*.txt"), emit: tupled_report

  script:
  def isPaired = params.SE ? '' : '--paired'
  """
  trim_galore $isPaired --basename ${id} -j $task.cpus ${reads.join(' ')}
  """
}
