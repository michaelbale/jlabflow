/*
 * Calculates RNA degradation metric from Smith, et al., & Sheffield (2021)
 * Scope: QPRO
 * Input: Cut-adapt (trim_galore!) output log with ID variable
 * Emits: One-line file with Degradation ratio
 * Feeds: Collect QC Metrics script to output QC data-table
 */
 





process DEGRADATION {
  tag "Calculating degradation ratio for ${id}"
  label 'small_mem'


  input:
  tuple val(id), path(log)

  output:
  path("${id}_degradation.log"), emit: degradation_logs

  script:
  """
  calcRNADegradation.R -r $log -o $id > ${id}_degradation.log
  """
}
