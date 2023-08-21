process NAMESORT {
  tag "Sorting bam ${id} by name for Genrich"
  label 'med_mem'
  
  input: 
  tuple val(id), path(bam)
  
  output:
  tuple val(id), path("${id}_nameSort.bam")
  
  script:
  """
  samtools sort -@ $task.cpus -n -o ${id}_nameSort.bam $bam
  """
}