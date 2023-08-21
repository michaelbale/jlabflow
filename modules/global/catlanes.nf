getSampleID = {
	    (it =~ /(.+)_S\d+_L\d{3}/)[0][1]
	}
    
	Channel
	  .fromFilePairs(params.input, flat: true)
	  .map { prefix, r1, r2 -> tuple(getSampleID(prefix), r1, r2) }
	  .groupTuple()
	  .set {inFq_ch}
	

      process catLanes {
	    tag "Concatenating lanes into $params.workDir"
		publishDir "$params.workDir/$sampleID", mode: 'copy', pattern: "*.gz"
		label 'small_mem'
		
		input:
		tuple val(sampleID), path(R1), path(R2) from inFq_ch
				
		output:
		tuple val(sampleID), path("${sampleID}_*_init.fq.gz") into reads_ch
		
		script:
		"""
		zcat $R1 > ${sampleID}_R1_init.fq
		zcat $R2 > ${sampleID}_R2_init.fq
		
		gzip ${sampleID}_R1_init.fq
		gzip ${sampleID}_R2_init.fq
		"""
	  }