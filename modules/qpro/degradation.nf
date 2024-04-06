/*
 * Calculates RNA degradation metric from Smith, et al., & Sheffield (2021)
 * Scope: QPRO
 * Input: Cut-adapt (trim_galore!) output log with ID variable
 * Emits: One-line file with Degradation ratio
 * Feeds: Collect QC Metrics script to output QC data-table
 */
 
process DEGRADATION {
  tag "Calculating degradation ratio for project: ${params.name}"
  label 'small_mem'


  input:
  path(logs)

  output:
  path("${params.name}_degradation.log"), emit: degradation_logs

  script:
  def myLogs = logs.join(',')
  """
  calcRNADegradation.R -r $myLogs > ${params.name}_degradation.txt
  """
}
