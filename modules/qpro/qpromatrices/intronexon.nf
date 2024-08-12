/*
 * Generates Exon:Intron ratio for QC table
 * Scope: QPRO
 * Input: Stranded bigwigs
 * Emits: Exon:Intron ratio for each sample
 * Feeds: Collect QC Metrics script to output QC data-table
 */
 
process INTRONEXON {
  tag "Calculating Exon:Intro ratio for project: ${params.name}"
  label 'med_mem'


  input:
  path(bw)

  output:
  path("${params.name}_ExIn.log"), emit: intronExon_logs
  
  script:
  def myBW = bw.join(',')
  """
  calcInEx.R \
    -b $myBW \
	-t $task.cpus \
	-g $params.genome > ${params.name}_ExIn.log
  """
}
