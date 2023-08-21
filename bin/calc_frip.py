#!/usr/bin/env python3
import argparse
import pysam
import sys
import deeptools.countReadsPerBin as crpb

def calcFrip(bamFiles, peakFile, cpus = 1):
  myRiP = [
    crpb.CountReadsPerBin(
      [i],
      bedFile = peakFile,
      numberOfProcessors = cpus
    ).run().sum(axis = 0)
    for i in bamFiles 
  ]
  
  myMapped = [
    pysam.AlignmentFile(i).mapped for i in bamFiles
  ]
  
  fRiP = dict(zip(bamFiles, [ i/j for i,j in zip(myRiP,myMapped) ]))
  
  return(fRiP)
  
def main(argv):

    parser = argparse.ArgumentParser(
    prog='calc_frip.py',
    description='''
      FOR USE WITH JLAB NEXTFLOW PIPELINE
      Reads bam file(s) to calculate the fraction of reads in peaks (fRiP)
      based on input peak file.
      Requires pysam and deeptools
    ''',
    epilog='''
      Email me at michaeljbale93@gmail.com with issues
    ''',
    add_help=True
    )

    parser.add_argument(
    '-i', '--bams', required = True, nargs = '+',
    help='''
      Bam files for reading
    '''
    )

    parser.add_argument(
    '-p', '--peakFile', required = True, nargs = 1,
    help='''
      Peak file for reading
    '''
    )

    parser.add_argument(
    '-o', '--outfile', required = True, nargs = 1,
    help='''
      File to write results to
    '''
    )

    parser.add_argument(
    '-t', '--threads', required = False, default = 1, type = int,
    help='''
      Number of processesors to use
    '''
    )
    
    args = vars(parser.parse_args(argv))
    
    frips = calcFrip(args['bams'], args['peakFile'], args['threads'])
    
    of = args['outfile'][0]
    
    with open(of, 'w') as f:
      for key, value in frips.items():
        f.write('%s: %s\n' % (key, value))
    f.close()


if __name__ == "__main__":
  main(sys.argv[1:])
  


        
      
      
