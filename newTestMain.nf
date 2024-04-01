#!/usr/bin/env nextflow 

nextflow.enable.dsl = 2

params.input = ''
params.genome = 'mm39'
params.name = 'JLABFLOW'


include { CATLANES } from './modules/global/newcat'
include { TRIM } from './modules/global/trim'

//include { QPRO } from './modules/qpro'



print(params.input)
		 
		 
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
	if(params.mode == 'qpro'){
          TRIM( reads )
	}
}



/* 
 * completion handler
 */
workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
}
