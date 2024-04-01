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
	    tag "Concatenating lanes into $params.workDir"
		publishDir "$params.workDir/$sampleID", mode: 'copy', pattern: "*.gz"
		label 'small_mem'
		
		input:
		tuple val(sampleID), path(reads)
		val(seqType)
		
		output:
		tuple val(sampleID), path("${sampleID}_*_init.fq.gz"), emit: reads_ch
		
		script:
		if( seqType == 'SE')
			"""
            echo $sampleID
            echo $reads
            echo $seqType
			cat ${reads} > ${sampleID}_R1_init.fq.gz
			"""
		else if( seqType == 'PE' )
			"""
			cat ${reads[0]} > ${sampleID}_R1_init.fq.gz
			cat ${reads[1]} > ${sampleID}_R2_init.fq.gz
			"""
	  }
