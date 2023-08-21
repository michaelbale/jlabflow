/*
 * Parse bam file for map rates along chromosomes prior to chrM removal
 * Scope: global
 * Input: Deduped bam file with group ID
 * Emits: idxstats log file
 * Feeds: multiqc
 */

process IDXSTATS {
  tag "Getting mapping stats for ${id}"
  label 'small_mem'
  
  input:
  tuple val(id), path(bam)
  
  output:
  path("${id}_idxStats.log"), emit: idx_log
  
  script:
  """
  samtools index ${bam}
  samtools idxstats ${bam} > ${id}_idxStats.log
  """
}