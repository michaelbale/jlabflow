/* Align qPRO data to index
 * Scope: qPRO
 * Input: bwa-mem2 index
 * Input: Trimmed fastq
 * Emits: Sorted Bam
 * Feeds: UMI dedup
 * Feeds: remove rDNA
 */
 
process BWAMEM2 {
  tag "mapping to ${idx} for ${params.genome}"
  label 'med_mem'
  
  input:
  tuple(val), path(reads)
  val(idx)
  
  output:
  path(