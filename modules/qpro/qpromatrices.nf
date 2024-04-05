/*
 * SUB-WORKFLOW: QPROMATICES
 * Generates various matrices from bigwig files for qpro analysis
 * Input: stranded bigwigs 
 * Output: Pause index matrix, median pausing index log
 * Output: Genebody matrix, TSS matrix
 * Output: Intron/Exon ratio 
 * Feeds: Collect qpro QC metrics
 */

include { GENEBODY } from './qpromatrices/genebody'
include { INTRONEXON } from './qpromatrices/intronexon'
include { PAUSEINDEX } from './qpromatrices/pauseindex'
include { TXSTARTSITE } from './qpromatrices/txstartsite'


workflow DATAVIZ {
  
  take:
  strandedBigwigs
  
  main:
  GENEBODY( strandedBigwigs )
  INTRONEXON( strandedBigwigs )
  PAUSEINDEX( strandedBigwigs )
  TXSTARTSITE( strandedBigwigs ) 
  
  emit:
  //stuff
  
  
}
