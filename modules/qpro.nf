//Genome specific
params.genome = ''
params.genomes = []
params.bt2_index = params.genome ? params.genomes[ params.genome ].bt2Index ?: false : false
params.rpgc = params.genome ? params.genomes[ params.genome ].rpgc ?: false : false
params.genesList = params.genome ? params.genomes[ params.genome ].genesList ?: false : false
params.annotations = params.genome ? params.genomes[ params.genome ].annotations?: false : false
params.forbid = params.genome ? params.genomes[ params.genome ].forbid ?: false : false

include { IDXSTATS } from './global/idxstats'
include { TRIM } from './global/trim'
include { FASTQC } from './global/fastqc'
include { BOWTIE2MAP } from './dna/bt2map'
include { FINALFILTER } from './dna/finalfilter'
include { BIGWIGSTRANDED as BIGWIGSTRANDEDPLUS } from './qpro/bigwigstranded'
include { BIGWIGSTRANDED as BIGWIGSTRANDEDMINUS } from './qpro/bigwigstranded'
include { DEGRADATION } from './qpro/degradation'
include { UMIDUPES } from './qpro/umidupes'
//include { DATAVIZ as DATAVIZPLUS, DATAVIZ as DATAVIZMINUS } from './dataviz'


def pipelineInfo() {
log.info """\
		===========================================================
            N A S C E N T   R N A - S E Q    P I P E L I N E v1
        ===========================================================
		Author: Michael J. Bale (mib4004@med.cornell.edu)
		
		Project ID: ${params.name}
        Genome: ${params.genome}
        Reads: ${params.input}
        Publish Directory: ${params.outdir}
        Filter List: ${params.forbid}
        """
         .stripIndent()
}
workflow QPRO {
  
    take:
	  data
	  
	main:
	  pipelineInfo()
	  TRIM( data )
          TRIM.out.trimmed_reads.view()
	  DEGRADATION( TRIM.out.tupled_report )
	  BOWTIE2MAP( params.bt2_index, TRIM.out.trimmed_reads )
	  FASTQC( TRIM.out.trimmed_reads )
	  UMIDUPES( BOWTIE2MAP.out.init_bt2 )
	  IDXSTATS( UMIDUPES.out.dedup_bam )
	  FINALFILTER( params.forbid, UMIDUPES.out.dedup_bam )
	  BIGWIGSTRANDEDPLUS( FINALFILTER.out.final_bams, 'plus' )
	  BIGWIGSTRANDEDMINUS( FINALFILTER.out.final_bams, 'minus' )
	  //QPROMATRICES( //mixed bigwigs )
	  //DATAVIZPLUS( BIGWIGRPGC.out.collect() , params.plusGenes , 'qpro' )
	  //DATAVIZMINUS( BIGWIGRPGC.out.collect() , params.minusGenes , 'qpro' )
	  //COLLECTQCMETRICS(
        //DEGRADATION.out.degradation_logs.collect(),
	//	QPROMATRICES.out.pause_log.collect(),
	//	QPROMATRICES.out.intronExon_log.collect(),
	//	FINALFILTER.out.forbid_list_count.collect(),
	//	FINALFILTER.out.final_count.collect()		
	  //)
	
	emit:
	  TRIM.out.trim_report
	    .mix( BOWTIE2MAP.out.bt2_logs )
	    .mix( FASTQC.out )
	    .mix( UMIDUPES.out.dedup_log )
	    .mix( IDXSTATS.out )
	    .collect()

}
