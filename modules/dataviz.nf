

include { MULTIBWSUMMARY } from './viz/multibwsummary'
include { PLOTCORREL } from './viz/plotcorrel'
include { PLOTPCA } from './viz/plotpca'
include { TORNADO } from './viz/tornado'
include { SUNSET } from './viz/sunset'


workflow DATAVIZ {
  
  take:
  bigwigs
  regions
  mode
  
  main:
  dir = mode == 'peaks' ? 'peakViz' : 'geneViz'
  MULTIBWSUMMARY( bigwigs, regions )
  PLOTPCA( MULTIBWSUMMARY.out, dir )
  PLOTCORREL( MULTIBWSUMMARY.out, dir)
  if(mode != 'peaks'){
    TORNADO( bigwigs, regions, dir )
  }
  if(mode != 'atac') {
    SUNSET( bigwigs, regions, dir )
  }
  
}