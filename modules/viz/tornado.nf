process TORNADO {
  tag 'Creating Tornado plot with ${regions}'
  label 'big_mem'
  publishDir "${params.outdir}/${dir}", mode: 'move'
  
  input:
  path(bigwigs)
  path(regions)
  val(dir)
  
  output:
  file("${params.name}_tornadoPlot.svg")
  
  script:
  """
  computeMatrix reference-point \
    -p $task.cpus \
	-S $bigwigs -R $regions \
	-o tmp.npz -b 1500 -a 1500 \
	--missingDataAsZero --smartLabels
  plotHeatmap \
    -m tmp.npz \
	--colorMap RdBu \
	--zMin auto --zMax auto \
	-o ${params.name}_tornadoPlot.svg
  """
}