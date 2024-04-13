/*
 * Create stranded bigwig file from nascent RNAseq
 * Scope: QPRO
 * Input: Final filtered bam file
 * Input: orientation
 * Emits: un-normalized bigwig
 * Feeds: workflow QPROSEQMATRICES
 * Feeds: workflow DATAVIZSTRANDED
 */
 
 
process BIGWIGSTRANDED {
  tag "creating $strand strand bigwig for $id"
  label 'big_mem'
  publishDir "$params.outdir/strandedbigwigs/${strand}Strand/", mode: 'copy', pattern: "${id}_${strand}.bw"
  
  input:
  tuple val(id), path(bam)
  val(strand)
  
  output:
  tuple path("${id}_${strand}.bw"), emit: stranded_bw
  
  script:
  def strandArg = (strand == 'plus') ? 'forward' : 'reverse'
  """
  samtools index -@ $task.cpus $bam
  bamCoverage -p $task.cpus \
    --bam $bam \
	-o ${id}_${strand}.bw \
	-bs 1 --normalizeUsing None \
        --skipNAs \
	--Offset 1 \
	--filterRNAstrand $strandArg
  """
}
