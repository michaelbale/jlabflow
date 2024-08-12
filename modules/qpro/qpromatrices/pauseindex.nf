/*
 * Generates pause index matrix and outputs median PI for QC table
 * Scope: QPRO
 * Input: Stranded bigwigs
 * Emits: PI matrix annotated with gene name and transcript ID
 * Emits: Median PI for each sample
 * Feeds: Collect QC Metrics script to output QC data-table
 */
 
process PAUSEINDEX {
  tag "Calculating Pause Indices for project: ${params.name}"
  publishDir "$params.outdir/matrices/", mode: 'move', pattern: "${params.name}_PImat.txt"
  label 'med_mem'


  input:
  path(bw)

  output:
  path("${params.name}_PI.log"), emit: pause_logs
  path("${params.name}_PImat.txt"), emit: pi_mat
  
  script:
  def myBW = bw.join(',')
  """
  calcPI.R \
    -b $myBW \
	-p "${params.name}" \
	-t $task.cpus \
	-g $params.genome \
	-f $params.filter > ${params.name}_degradation.log
  """
}
