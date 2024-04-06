#!/usr/bin/env nextflow 

nextflow.enable.dsl = 2

params.input = ''
params.genome = 'mm39'
params.name = 'JLABFLOW'

def helpUtils = evaluate(new File("${projectDir}/utils/help.groovy"))
if(params.help || params.mode == '') {
  helpUtils.showHelpMessage(params.mode)
  System.exit(0)
}

include { ATACSEQ } from './modules/atacseq'
include { IPSEQ } from './modules/ipseq'
include { QPRO } from './modules/qpro'
include { PREPAREINPUT } from './modules/prepareinput'
include { MULTIQC } from './modules/global/multiqc'


		
		 
workflow {
  
  
  main:
	
	PREPAREINPUT( params.input )
    if(params.mode == 'atac'){
	  if( params.SE ) {
	    error "ERROR: params mode == atac and SE == TRUE are incompatible; ATACSEQ should always be PE seq"
	  }
	  ATACSEQ( PREPAREINPUT.out )
	  logs = ATACSEQ.out
	} else if( params.mode == 'ip' || params.mode == 'cnr' || params.mode == 'cnt') {
	  if(params.mode == 'cnt' && params.SE) { 
	    error "ERROR: params mode == cnt and SE == TRUE are incompatible; cnt-IPSEQ should always be PE seq" 
	  }
      IPSEQ( PREPAREINPUT.out )
	  logs = IPSEQ.out	
	} else if (params.mode == 'rna') {
	  print("Not yet supported")
	} else if(params.mode == 'qpro') {
	  if( !params.SE ) { error "ERROR: qPRO must be SE - don't make me do more work plz (:" }
	  QPRO( PREPAREINPUT.out )
          logs = QPRO.out
	}
	
  MULTIQC( logs )
}

/* 
 * completion handler
 */
workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
}
