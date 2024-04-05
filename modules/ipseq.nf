params.genome = ''
params.genomes = []
params.bt2_index = params.genome ? params.genomes[ params.genome ].bt2Index ?: false : false
params.forbid = params.genome ? params.genomes[ params.genome ].forbid ?: false : false
params.rpgc = params.genome ? params.genomes[ params.genome ].rpgc ?: false : false
params.genesList = params.genome ? params.genomes[ params.genome ].genesList ?: false : false
params.annotations = params.genome ? params.genomes[ params.genome ].annotations?: false : false

include { IDXSTATS } from './global/idxstats'
include { TRIM } from './global/trim'
include { FASTQC } from './global/fastqc'

include { BOWTIE2MAP } from './dna/bt2map'
include { RMDUPES } from './dna/rmdupes'
include { FINALFILTER } from './dna/finalfilter'
include { BIGWIGRPGC } from './dna/bigwigbasic'
include { DATAVIZ } from './dataviz'

def pipelineInfo() {
log.info """\
		=================================================
            I P S E Q   P I P E L I N E v1
        =================================================
		Author: Michael J. Bale (mib4004@med.cornell.edu)
		
		Project ID: ${params.name}
        Genome: ${params.genome}
        Reads: ${params.input}
        Publish Directory: ${params.outdir}
        Forbidden List: ${params.forbid}
        BT2: ${params.bt2_index}
        """
         .stripIndent()
}
		 
workflow IPSEQ {
  
    take:
	  data
	  
	main:
	  pipelineInfo()
	  TRIM( data )
	  BOWTIE2MAP( params.bt2_index, TRIM.out.trimmed_reads )
	  FASTQC( TRIM.out.trimmed_reads )
	  RMDUPES( BOWTIE2MAP.out.init_bt2 )
	  IDXSTATS( RMDUPES.out.dedup_bam )
	  FINALFILTER( params.forbid, RMDUPES.out.dedup_bam )
	  BIGWIGRPGC( FINALFILTER.out )
  	  DATAVIZ( BIGWIGRPGC.out.collect() , params.genesList, 'ip')

	
	emit:
	  TRIM.out.trim_report
	    .mix( BOWTIE2MAP.out.bt2_logs )
	    .mix( FASTQC.out )
	    .mix( RMDUPES.out.dedup_log )
	    .mix( IDXSTATS.out )
	    .collect()

}
