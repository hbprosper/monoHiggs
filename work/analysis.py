#!/usr/bin/env python
# ---------------------------------------------------------------------------
import os, sys
from ROOT import gSystem
# ---------------------------------------------------------------------------
def main():
    gSystem.Load("libmonohiggs")
    from ROOT import monoHiggs
    
    inputFile = "Higgs_hhxx_scalar_nohdecay_1GeV_13TeV.lhe_4leptons_CMS.root"
    prefix  = "hzz4leptons"
    pileup  = 50
    nevents = 10
    monoHiggs.analysis4mu(inputFile, prefix, pileup)
    
# ---------------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print "ciao!"
