process COLLECTFRIP {
  tag 'summarizing frip scores for ${type} peaks'
  label 'small_mem'
  publishDir "${params.outdir}/", mode: 'move'
  
  input:
  val(type)
  val(dir)
  path(friplogs)
  
  output:
  file("summarized_${type}FRIPscores.txt")
  
  script:
  """
  echo -e "Sample:\tfRiP" > summarized_${type}FRIPscores.txt
  cat $friplogs >> summarized_${type}FRIPscores.txt
  """
}
  
