/*
 * Collect all log reports from various software into html (multiqc)
 * Input: ALL LOGS
 * Emits: multiqc reports
 */


process MULTIQC {
  tag "collating logs with Multiqc"
  label "small_mem"
  publishDir "$params.outdir", mode: 'move'
  
  input:
  path("*")
  
  output:
  path('multiqc_report.html')

  script:
  """
  multiqc .
  """
}
