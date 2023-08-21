/*
 * SUB-WORKFLOW: DATAVIZ
 * Generates various visualizations from bigwig files with supplied region file(s)
 * Input: bigwigs
 * Input: region file
 * Input: mode; chooses which plots to make
 * Output: PCA Plot, Correl plot
 * Output: Tornado Plot, Sunset Plot
 * Feeds: Nothing
 */

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
