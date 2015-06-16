#!/usr/bin/env python
# ---------------------------------------------------------------------------
import os, sys
from ROOT import gSystem
# ---------------------------------------------------------------------------
def main():
    gSystem.Load("libmonohiggs")
    from ROOT import monoHiggs
    
    inputFile = sys.argv[1] #"Higgs_hhxx_scalar_nohdecay_1GeV_13TeV.lhe_4leptons_CMS.root"
    prefix  = sys.argv[2] #"hzz4leptons"
    pileup  = 50
    nevents = int(sys.argv[3]) #10
    monoHiggs.analysis4mu(inputFile, prefix, pileup, nevents)
    
# ---------------------------------------------------------------------------
try:
    if len(sys.argv) != 4:
        exit("Usage: ./analyis.py [.root file] [prefix] [nevents]")

    main()
except KeyboardInterrupt:
    print "ciao!"
