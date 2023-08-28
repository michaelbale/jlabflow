process SUMMARIZEGENOMICLOCALIZATION {
  tag "Creating Peak localization with ${params.genome}"
  label 'small_mem'
  publishDir "${params.outdir}/callpeaks"

  input:
  path(peaks)
  path(annotations)

  output:
  path("${params.name}_peaksSummarized.pdf")

  script:
  def myPeaks = peaks.join(",")
  """
  summarizePeakLocalization.R -p $myPeaks -a $annotations  -t pdf -o "${params.name}_peaksSummarized"
  """
}