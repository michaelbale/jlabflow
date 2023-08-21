process PLOTPCA {
  tag 'Plotting Data Correlogram'
  label 'small_mem'
  publishDir "${params.outdir}/${dir}", mode: 'move'
  
  input:
  path(ctMat)
  val(dir)
  
  output:
  file("${params.name}_pca.svg")
  
  script:
  """
  plotPCA \
    -in $ctMat \
	-T "${params.name} Project Correlogram" \
	--plotFileFormat svg \
	-o ${params.name}_pca.svg
  """
}