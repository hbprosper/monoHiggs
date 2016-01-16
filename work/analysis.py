#!/usr/bin/env python
# ---------------------------------------------------------------------------
import os, sys
from ROOT import gSystem
# ---------------------------------------------------------------------------
def main():
    gSystem.Load("libmonohiggs")
    from ROOT import monoHiggs
    
    inputFile = sys.argv[1] 
    sample    = sys.argv[2]
    if len(sys.argv) > 3: 
        nevents = int(sys.argv[3])
    else:
        nevents = -1 # all events

    if len(sys.argv) > 4: 
        luminosity = float(sys.argv[4])
    else:
        luminosity = 30.0 # 1/fb

    if len(sys.argv) > 5: 
        xsection = float(sys.argv[5])
    else:
        xsection = 1.0    # fb

    if len(sys.argv) > 6: 
        pileup = int(sys.argv[6])
    else:
        pileup = 0        

        
    monoHiggs.hzz4l(inputFile, sample, nevents,
                    luminosity, xsection, pileup)
    
# ---------------------------------------------------------------------------
try:
    if len(sys.argv) < 3:
        exit('''
    Usage: ./analyis.py Delphes-Root-file sample
                        [nevents=-1]
                        [luminosity=30/fb]
                        [xsection=1fb]
                        [pileup=0]
        
        Note: the sample must identify the final state (4mu, 4e, or 2e2mu)
        and the sample type (e.g., sig or bkg). It will be used as the
        prefix for the name of the root output file containing histograms.

        example: sample = s_4mu''')

    main()
except KeyboardInterrupt:
    print "ciao!"
