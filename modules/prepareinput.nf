/*
 * SUB-WORKFLOW: PREPAREINPUT
 * Modularization for preparing input data
 * Input: params.input
 * Output: reads prepped for workflows
 * Feeds: Base Workflow
 */

include { CATLANES } from './global/catlanes'


workflow PREPAREINPUT {

  take:
    input
  
  main:
    if(params.catlanes){
	  Channel
	    .fromPath(input)
		.map { file ->
		  def sampleName = file.baseName.split('_')[0]
		  def isPE = !params.SE 
		  return tuple(sampleName, file, isPE)
		}
		.groupTuple(by: 0) 
		.map { sampleName, files, isPE ->
		  def reads1 = files.findAll { it.baseName.contains('_R1_') }
		  def reads2 = isPE ? files.findAll { it.baseName.contains('_R2_') } : []
		  return tuple(sampleName, reads1, (isPE ? reads2 : []))
		}
		.set { initFq }
		reads = CATLANES( initFq )
    } else if (params.SE) {
      reads = Channel
        .fromPath(input)
        .map {
          file ->
            def sampleName = file.baseName.split('_')[0]
            return [sampleName, file]
        }
    } else {
        reads = Channel.fromFilePairs(input, checkIfExists: true)
    }
	
  emit:
    reads
}
