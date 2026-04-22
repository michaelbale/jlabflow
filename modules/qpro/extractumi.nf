/*
 * Process UMI-containing reads prior to mapping
 * Scope: qPRO
 * Input: trimmed reads with grouping ID
 * Emits: Processed reads
 * Feeds: Bowtie2 mapping
 */

process EXTRACTUMI {
  tag "Extracting UMI from $id"
  label 'med_mem'
  
  input:
  tuple val(id), path(reads)
  
  output:
  tuple val(id), path("${id}*.fq.gz"), emit: umi_extract_reads
  
  script:
  def inputArgs = params.SE 
    ? "-I ${reads} \
       -S ${id}_umiTools-processed.fq.gz" 
    : "-I ${reads[0]} \
       -S ${id}_R1_umiTools-processed.fq.gz \
       --read2-in ${reads[1]} \
       --read2-out ${id}_R2_umiTools-processed.fq.gz"
  def bcPattern = params.SE ? '' : "--bc-pattern2 NNNNNNX"
  """
  umi_tools extract \
    --umi-separator=':' \
    -p "NNNNNNX" \
    $bcPattern \
    $inputArgs
  """
}
