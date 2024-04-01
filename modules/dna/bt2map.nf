/*
 * Map trimmed fastqc to genome (bowtie2)
 * Scope: DNA
 * Input: trimmed reads with grouping ID
 * Input: Pre-indexed Genome from command line
 * Emits: Mapping statistics log
 * Emits: Sorted bam file
 * Emits: Index file for bam
 * Feeds: Duplication Removal
 */

process BOWTIE2MAP {
  tag "Aligning $id to ${params.genome}"
  label 'big_mem'
  
  input:
  val(idx)
  tuple val(id), path(reads)
  
  output:
  path("${id}_bt2.log"), emit: bt2_logs
  tuple val(id), file("${id}_iSort.bam"), file("${id}_iSort.bam.bai"), emit: init_bt2
  
  script:
  def inputArgs = params.SE ? "-U $reads" : "-1 ${reads[0]} -2 ${reads[1]}"
  def dovetail = (params.mode == 'cnr' && !params.SE) ? "--dovetail" : ''
  """
  bowtie2 \
    -p $task.cpus \
    -x ${idx} \
    --no-mixed --no-unal --no-discordant \
    --local --very-sensitive-local \
    -X 1000 -k 4 --mm \
	$dovetail \
	$inputArgs 2> ${id}_bt2.log  | samtools view -bS -q 30 - > ${id}_init.bam
  samtools sort -@ $task.cpus ${id}_init.bam > ${id}_iSort.bam
  samtools index ${id}_iSort.bam
  """
}