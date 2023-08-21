/*
 * Remove Duplicated reads (picard MarkDuplicates)
 * Scope: DNA
 * Input: indexed bam file with grouping ID
 * Emits: Deduplication rate log
 * Emits: Dedeuplicated Bam
 * Feeds: Final Filter Step
 */

process RMDUPES {
   tag "Removing Dupes from ${id}"
   label 'med_mem'
   
   input:
   tuple val(id), path(bam), path(index)
   
   output:
   path("${id}_dups.log"), emit: dedup_log
   tuple val(id), file("${id}_rmDup.bam"), emit: dedup_bam
   
   script:
   """
   picard MarkDuplicates VERBOSITY=WARNING \
	INPUT=${bam} OUTPUT=${id}_rmDup.bam \
	METRICS_FILE=${id}_dups.log \
	REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT
   """
}
