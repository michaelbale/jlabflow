/*
 * Concatenate lanes from sequencer
 * Scope: GLOBAL
 * Input: Fastq files with grouping ID
 * Input: If PE or SE sequencing
 * Emits: Merged Fastq files
 * Emits: Grouping ID
 * Feeds: Main Workflow Selected
 */

process CATLANES {
    tag "Concatenating lanes for ${sampleID}"
    publishDir "${params.workDir}/${sampleID}", mode: 'copy', pattern: "*.gz"

    input:
    tuple val(sampleID), path(reads1), path(reads2) // reads2 can be an empty list for SE

    output:
    tuple val(sampleID), path("${sampleID}_*_combined.fq.gz"), optional: true, emit: combinedReads

    script:
    def r2Exists = (reads2 && !reads2.isEmpty()) ? true : false
    """
    zcat ${reads1.join(' ')} > ${sampleID}_R1_combined.fq
	gzip ${sampleID}_R1_combined.fq
    ${r2Exists ? "zcat ${reads2.join(' ')} > ${sampleID}_R2_combined.fq; gzip ${sampleID}_R2_combined.fq" : ''}
    """
}

