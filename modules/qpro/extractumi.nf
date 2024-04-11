/*
 * Process UMI-containing reads prior to mapping
 * Scope: qPRO
 * Input: trimmed reads with grouping ID
 * Emits: Processed reads
 * Feeds: Bowtie2 mapping
 */

process EXTRACTUMI {
  tag "Aligning $id to ${params.genome}"
  label 'med_mem'
  
  input:
  tuple val(id), path(reads)
  
  output:
  tuple val(id), path("${id}_umiTools-processed.fq.gz"), emit: umi_extract_reads
  
  script:
  """
  umi_tools \
    extract \
	--umi-separator=':' \
	-p NNNNNN \
	-I ${reads} \
	-S ${id}_umiTools-processed.fq.gz
  """
}
