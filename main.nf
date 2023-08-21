#!/usr/bin/env nextflow 

nextflow.enable.dsl = 2

params.input = './testAtac/*{R1,R2}*'
params.genome = 'mm39'
params.name = 'oh please work'


include { ATACSEQ } from './modules/atacseq'
include { IPSEQ } from './modules/ipseq'
include { MULTIQC } from './modules/global/multiqc'



		 
		 
workflow {
  
  
  main:
    if(params.mode == 'atac'){
	  read_pairs = channel.fromFilePairs( params.input, checkIfExists: true ) 
	  ATACSEQ( read_pairs )
	  logs = ATACSEQ.out
	} else if( params.mode == 'ip') {
      read_pairs = channel.fromFilePairs( params.input, checkIfExists: true ) 
	  IPSEQ( read_pairs )
	  logs = IPSEQ.out	
	} else if (params.mode == 'rna') {
	  print("Not yet supported")
	}
  MULTIQC( logs )
}

/* 
 * completion handler
 */
workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
}
