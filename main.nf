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
//include { QPRO } from './modules/qpro'
include { CATLANES } from './modules/global/catlanes'
include { MULTIQC } from './modules/global/multiqc'


		
		 
workflow {
  
  
  main:
	if(params.catlanes){
	  Channel
	    .fromPath(params.input)
		.map { file ->
		  def sampleName = file.baseName.split('_')[0]
		  def isPE = !params.SE 
		  return tuple(sampleName, file, isPE)
		}
		.groupTuple(by: 0, size: 3) 
		.map { sampleName, files, isPE ->
		  def reads1 = files.findAll { it.baseName.contains('_R1_') }
		  def reads2 = isPE ? files.findAll { it.baseName.contains('_R2_') } : []
		  return tuple(sampleName, reads1, (isPE ? reads2 : []))
		}
		.set { initFq }
		reads = CATLANES( initFq )
    } else if (params.SE) {
      reads = Channel
        .fromPath(params.input)
        .map {
          file ->
            def sampleName = file.baseName.split('_')[0]
            return [sampleName, file]
        }
    } else {
        reads = Channel.fromFilePairs(params.input, checkIfExists: true)
    }	
	print(params.mode)
    if(params.mode == 'atac'){
	  if( params.SE ) {
	    error "ERROR: params mode == atac and SE == TRUE are incompatible; ATACSEQ should always be PE seq"
	  }
	  ATACSEQ( reads )
	  logs = ATACSEQ.out
	} else if( params.mode == 'ip' || params.mode == 'cnr' || params.mode == 'cnt') {
	  if(params.mode == 'cnt' && params.SE) { 
	    error "ERROR: params mode == cnt and SE == TRUE are incompatible; cnt-IPSEQ should always be PE seq" 
	  }
      IPSEQ( reads )
	  logs = IPSEQ.out	
	} else if (params.mode == 'rna') {
	  print("Not yet supported")
	} else if(params.mode == 'qpro') {
	  print('Not yet supported')
	}
	
  MULTIQC( logs )
}

/* 
 * completion handler
 */
workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
}
