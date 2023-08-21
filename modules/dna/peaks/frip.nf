process FRIP {
  tag "Calculating fRiP for ${id} with set ${peaks}"
  label 'med_mem'
  
  input:
  tuple val(id), path(bams), path(peaks)
  
  output:
  path("${id}_frip.txt"), emit: friplog
  
  script:
  """
  for i in $bams; do 
    samtools sort -@ $task.cpus -o \${i%_nameSort.bam}_tmpSort.bam \$i;
    samtools index \${i%_nameSort.bam}_tmpSort.bam;
  done
  calc_frip.py -i *tmpSort.bam -p $peaks -t $task.cpus -o ${id}_frip.txt
  """
}