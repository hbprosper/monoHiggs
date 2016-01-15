#!/usr/bin/env python
# ---------------------------------------------------------------------------
import os, sys
from ROOT import gSystem
# ---------------------------------------------------------------------------
def main():
    gSystem.Load("libmonohiggs")
    from ROOT import monoHiggs
    
    inputFile = sys.argv[1] 
    prefix  = sys.argv[2] 
    pileup  = 50
    nevents = int(sys.argv[3])
    monoHiggs.hzz4l(inputFile, prefix, pileup, nevents)
    
# ---------------------------------------------------------------------------
try:
    if len(sys.argv) != 4:
        exit('''
    Usage: ./analyis.py [.root file] [prefix] [nevents]
        
        Note: the prefix must identify the final state (4mu, 4e, or 2e2mu)
        and the sample type (e.g., sig or bkg). It will be used as the
        prefix for the name of the root output file containing histograms.

        example: prefix = s_4mu''')

    main()
except KeyboardInterrupt:
    print "ciao!"
