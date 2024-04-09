/*
 * Deduplicate reads from sequences including UMI dogtags
 * Scope: QPRO
 * Input: Indexed initial bam file with grouping ID
 * Emits: deduplicated bam file
 * Emits: Dedup log file
 * Feeds: Final Filter step
 */
 
 
process UMIDUPES {
  tag "Performing UMI-aware deduplication on ${id}"
  label 'med_mem'
  
  input:
  tuple val(id), path(bam), path(index)
  
  output:
  tuple val(id), path("${id}_rmDup.bam"), emit: dedup_bam
  path("${id}_dedup.log"), emit: dedup_log
  
  script:
  """
  umi_tools dedup \
    -I $bam \
	--umi-separator=":" \
	-L ${id}_dedup.log \
	-S ${id}_rmDup.bam
	
  """
}
