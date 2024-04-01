/*
 * Create rudimentary bigwig for visualizations in IGV (deeptools)
 * Input: final bam with sample ID
 * Emits: CPM or RPGC-normalized bigwig
 * Feeds: Heatmap plot generator(s)
 */

process BIGWIGRPGC{

	tag "Creating ${sampleID} bigwig in RPGC mode"
	publishDir "$params.outdir/bigwig_norm1x", mode: 'copy'
	label 'big_mem'

	input:
	tuple val(sampleID), file(finalBam) 
	
	output:
	path("${sampleID}_normTo1x.bw")

	script:
	def extendReads = params.SE ? '' : '--extendReads'
	"""
	samtools index $finalBam
	bamCoverage -p $task.cpus \
	  --bam ${finalBam} \
	  -o ${sampleID}_normTo1x.bw \
	  -bs 10 ${extendReads} --smoothLength 50 \
	  --normalizeUsing RPGC --effectiveGenomeSize ${params.rpgc} \
	  --skipNonCoveredRegions 
	"""
}
