process GENRICH {
  tag "calling peaks using Genrich on ${id}"
  label 'med_mem'
  publishDir "$params.outdir/$dir", pattern: '*.narrowPeak', mode: 'copy'
  
  
  input:
  tuple val(id), path(bam)
  val(dir)
  
  output:
  tuple val(id), path(bam), path("${id}_genrich.narrowPeak"), emit: peak_tuples
  path("${id}_genrich.narrowPeak"), emit: peak_files
  
  script:
  def myOpts = bam.join(",")
  
  """
  Genrich -j -t ${myOpts} -o ${id}_genrich.narrowPeak
  """
}


