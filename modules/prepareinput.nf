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

    // Step 1: group input files by sampleID
    Channel.fromPath(input, checkIfExists: true)
        .map { file ->
            def sampleID = file.baseName.split('_')[0]
            return tuple(sampleID, file)
        }
        .groupTuple(by: 0)
        .map { sampleID, files ->
            files = files.sort { it.name }
            def reads1 = files.findAll { it.baseName.contains('_R1_') }
            def reads2 = files.findAll { it.baseName.contains('_R2_') }
            tuple(sampleID, reads1, reads2)
        }
        .set { groupedReads }

    // Step 2: conditionally run CATLANES for multi-lane data
    reads = params.catLanes
        ? CATLANES(groupedReads)       // multi-lane: concatenate
        : groupedReads                  // single-lane: pass through unchanged

    emit:
        reads
}
