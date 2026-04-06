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
    tuple val(sampleID), path(reads1), path(reads2)

    output:
    tuple val(sampleID), path("${sampleID}_*_combined.fq.gz"), optional: true, emit: combinedReads

    script:
    // Convert Path objects to strings
    def r1files = reads1.collect { it.toString() }.join(' ')
    def r2files = reads2.collect { it.toString() }.join(' ')

    // Only include R2 commands if reads2 is non-empty
    def r2Cmd = reads2 ? """
        zcat ${r2files} > ${sampleID}_R2_combined.fq
        gzip ${sampleID}_R2_combined.fq
    """ : ''

    """
    zcat ${r1files} > ${sampleID}_R1_combined.fq
    gzip ${sampleID}_R1_combined.fq
    ${r2Cmd}
    """
}
