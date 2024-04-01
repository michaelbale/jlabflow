#!/usr/bin/env nextflow 

nextflow.enable.dsl = 2

params.input = ''
params.genome = 'mm39'
params.name = 'JLABFLOW'


include { CATLANES } from './modules/global/catlanes'



print(params.input)
		 
		 
workflow {
    main:
    if (params.SE) {
        print('successfully parsed as SE')
        if (params.catlanes) {
            print('successfully parsed needing to catLanesSE')
            initFq = Channel
                .fromPath(params.input)
                .map { path -> 
                    def sampleName = path.baseName.split('_')[0]
					def r2 = []
                    return [sampleName, path] // Assuming path is what you meant to use
                }
				.groupTuple()
            reads = CATLANES(initFq, 'SE')  
        } else {
            print('no need to SE cat')
            reads = Channel
                .fromPath(params.input)
                .map {
                    file ->
                    def sampleName = file.baseName.split('_')[0]
                    return [sampleName, file]
                }
        }
    } else {
        print('successfully saw PE')
        if (params.catlanes) {
            print('successfully saw PE and catLanes')
            initFq = Channel
                .fromFilePairs(params.input, flat: true)
                .map{ prefix, r1, r2 -> 
				  def sampleName = prefix.split('_')[0]
				  return [sampleName, r1, r2]
				}
                .groupTuple()
            initFq.view()
            reads = CATLANES(initFq, 'PE')
        } else {
            reads = Channel.fromFilePairs(params.input, checkIfExists: true)
        }
    }
}


/* 
 * completion handler
 */
workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
}
