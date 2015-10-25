#!/usr/bin/env python

import sys,os
import ROOT
import math
from ROOT import TFile,TTree

rootfile = sys.argv[1]
print " ---> Opening file: ",rootfile
froot = TFile.Open(rootfile)
txAna = froot.Get("ggNtuplizer/EventTree")
if txAna == None :
    txAna = froot.Get("EventTree")
nEvts = txAna.GetEntries()
print "nEvtsInTree: %d" % nEvts

if len(sys.argv) > 2:
    nJobs = float(sys.argv[2])
    nEvtsPerJob = math.ceil(nEvts / nJobs)
    print "nEvtsPerJob: %d" % nEvtsPerJob
