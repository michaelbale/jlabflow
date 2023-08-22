/*
 * Get basic sequence QC Stats (fastqc)
 * Scope: GLOBAL
 * Input: trimmed reads with grouping ID
 * Emits: Log folder containing Fastqc output
 * Feeds: MultiQC
 */
 
process FASTQC {
  tag "FASTQC on ${id}"
  label 'med_mem'
  
  
  input:
  tuple val(id), path(reads)
  
  output:
  path("fastqc_${id}_logs"), emit: fqc_logs
  
  script:
  """
  mkdir fastqc_${id}_logs
  fastqc -o fastqc_${id}_logs -f fastq -t $task.cpus -q ${reads}
  """
}
