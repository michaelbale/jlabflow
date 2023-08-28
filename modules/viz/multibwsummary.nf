process MULTIBWSUMMARY{
  tag "Creating ctMat from $params.name with $regions"
  label 'big_mem'
  
  input:
  path(bigwigs)
  path(regions)
  
  output:
  path("${params.name}_bwSum.npz")
  
  script:
  """
  multiBigwigSummary BED-file \
    -b $bigwigs --BED $regions \
	-p $task.cpus --smartLabels \
	-out "${params.name}_bwSum.npz"
  """
}