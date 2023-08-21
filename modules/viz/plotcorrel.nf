process PLOTCORREL {
  tag 'Plotting Data Correlogram'
  label 'small_mem'
  publishDir "${params.outdir}/${dir}", mode: 'move'
  
  input:
  path(ctMat)
  val(dir)
  
  output:
  file("${params.name}_correl.svg")
  
  script:
  """
  plotCorrelation \
    -in $ctMat \
	-c spearman \
	-p heatmap \
	-T "${params.name} Project Correlogram" \
	--plotFileFormat svg \
	--removeOutliers \
	--plotNumbers \
	-o ${params.name}_correl.svg
  """
}