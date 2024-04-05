
include { GENRICH } from './peaks/genrich'
include { GENRICH as GENRICHREP } from './peaks/genrich'
include { FRIP } from './peaks/frip'
include { FRIP as FRIPALL } from './peaks/frip'
include { NAMESORT } from './peaks/namesort'
include { COLLECTFRIP } from './peaks/processfrip'
include { SUMMARIZEGENOMICLOCALIZATION } from './peaks/summarizepeaks'


workflow CALLPEAKS {

  take:
  mode
  bams
  
  main:
  
  def getGroupID = {
    (it =~ /(.+)_\d+/)[0][1]
  }
  
  NAMESORT( bams )
  
  GENRICH( NAMESORT.out, 'Genrich_indPeaks' )
  FRIP( GENRICH.out.peak_tuples )
  SUMMARIZEGENOMICLOCALIZATION( GENRICH.out.peak_files.collect(), params.annotations )
  COLLECTFRIP( 'individual', FRIP.out.collect() )
  if(params.callConsensus){
    forGREP = NAMESORT.out.map{ arr -> tuple(getGroupID(arr[0]), arr[1]) }.groupTuple()
	forGREP.view()
   // GENRICHREP( forGREP, 'results' )
	// FRIPALL( GENRICHREP.out, true )
  }
  
}
