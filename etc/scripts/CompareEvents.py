from sys import argv
import json 
from pprint import pprint 
from math import fabs
from math import cosh
import operator
import array
import sys
import ROOT

from optparse import OptionParser


parser = OptionParser()
parser.add_option("-d", "--doDifference",default=False, action="store_true", dest="doDiff")
parser.add_option("-i", "--doIntersection",default=False, action="store_true", dest="doIntersection")
parser.add_option("--compareAllVar", default=False, action="store_true", dest="compareAllV")


(options, args) = parser.parse_args()

#run D
#minrun = 203777
#maxrun = 208686

#run C
#minrun = 198049
#maxrun = 203742

#run A+B
#minrun = 190645
#maxrun = 196531

#all run 
minrun = 0
maxrun = 10000000

evtDump1 = {}
evtDump2 = {}

ROOT.gROOT.SetStyle("Plain")
ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetOptStat(111111)



#****************************************************


def getVariablesInTxt(inputFileName):

    evtDump = {}

    file = open(inputFileName, "r")
    
    for line in file:

        varList = {}

        try:
            
            for word in line.split():
           
                    variable, value = word.split(":")
                                       
                    if variable == "run" or variable == "lumi" or variable == "event" or variable == "type":
                        globals()[variable] = int(value)
                        
                    else:
                        varList[variable] = float(value)
                        if variable == "diphoMVA":
                            globals()[variable] = float(value)
                        

            if run < minrun or run > maxrun: continue
            if event < 0: continue   # filthy fix for globe events with neg. id
                    
            evtDump[ (run, lumi, event) ] = varList


        except Exception, e:
            print e
                    
    return evtDump



#*************** read input dumps and check for missing evts ***********************************************************


evtDump1 = getVariablesInTxt(args[0])
evtDump2 = getVariablesInTxt(args[1])
 
evtList1 = evtDump1.keys()
evtList2 = evtDump2.keys()

print "list1 = ", len(evtList1), " list2 = ", len(evtList2)

commonEvts = set(evtList1).intersection(set(evtList2))
missingEvts = list(set(evtList1) - set(evtList2))
surplus     = list(set(evtList2) - set(evtList1))

#print "intersection = ",len(commonEvts), " missing evts in list 2 = ", len(missingEvts), " surplus = ", len(surplus)


if len(evtList1) != 0:
    fract = float( len(missingEvts) )/float( len(evtList1) )    
    #print "Fraction of events missed by 2 in run ", minrun, "-", maxrun, ": ", fract

    
    
outputFile = open(args[2],"w")

if options.doDiff:
    for ievt in missingEvts:
        outputFile.write("run:%d lumi:%d event:%d \n" % ievt)
        
if options.doIntersection:
    for ievt in commonEvts:
        outputFile.write("run:%d lumi:%d event:%d \n" % ievt)

outputFile.close()



#h_Mass = ROOT.TH1F( "hMass", "diphoMass", 200, -4.0, 4.0 )

#for isevt in surplus:

    #allvari = evtDump2[isevt]
    #if allvari["ucat"] != -1:
        #print isevt, " cat:", allvari["tcat"], " untag:", allvari["ucat"]


for icat in range(8,10):

    globeVBF = []
    #for iallevt in evtDump1:
    for iallevt in missingEvts:
        allvarglobe = evtDump1[iallevt]
        if allvarglobe["cat"] == icat:
            globeVBF.append(iallevt)
            #print iallevt, " ", allvarglobe["cat"]
    print "number of missing events for globe cat ", icat, ":", len(globeVBF)
    

    evtinOne = []
    evtinTwo = []
    evt12 = []
    evtDiffOne = []
    evtDiffTwo = []
    
    for icevt in commonEvts:

        allvar1 = evtDump1[icevt]
        allvar2 = evtDump2[icevt]

        if allvar1["cat"] == icat:
            
            evtinOne.append(icevt)
            if allvar2["cat"] != icat:
                evtDiffOne.append(icevt)
                
                #print icevt, " cat:", allvar1["cat"], "-", allvar2["tcat"]
               # print icevt, " cat:", allvar1["cat"], "-", allvar2["cat"], " mass:", allvar1["mass"], "-", allvar2["mass"], " njets:", allvar1["numJets"], "-", allvar2["numJets"], " jet1 eta:", allvar1["jet1_eta"], "-", allvar2["jet1_eta"], " jet2 eta:", allvar1["jet2_eta"], "-", allvar2["jet2_eta"], " zeppenfeld:", allvar2["zeppenfeld"],  " dPhiGGJJ:", allvar2["dijdihdphi"], "combiMVA:", allvar1["combiMVA"], "-",  allvar2["combiMVA"], " dijetMVA:", allvar1["dijetMVA"], "-",  allvar2["dijetMVA"], " vbf_massJJ:", allvar1["vbf_massJJ"], "-",  allvar2["vbf_massJJ"], " dipho:", allvar1["diphoMVA"]

            
            
        #if allvar2["ucat"] == icat and allvar2["cat"] == -1:
        if allvar2["cat"] == icat:

            #print icevt, " cat:", allvar2["cat"]
            evtinTwo.append(icevt)
            if allvar1["cat"] != icat:
                evtDiffTwo.append(icevt)
                

    print "********** number of events in category:", icat, " for Globe = ", len(evtinOne), " for ggA = ", len(evtinTwo),"  diff = ", abs(float(len(evtinTwo))-float(len(evtinOne)))/float(len(evtinOne))*100, " %"
    print "************* evts in Globe and NOT in ggA = ", len(evtDiffOne), " in ggA and NOT in Globe = ", len(evtDiffTwo)

    #if allvar2["cat"] == 7 or allvar2["cat"] == 6 or allvar2["cat"] == 5:
     #   h_Mass.Fill( allvar2["mass"]-allvar1["mass"] )
    
    #if abs(sigmaM - sigmaO)/sigmaM > 0.003 and abs(allvar1["mass"] - allvar2["mass"])/allvar1["mass"] < 0.003:

       #print icevt, " cat:", allvar1["cat"], "-", allvar2["ucat"], " untag:",allvar2["ucat"], " pho1_eta:",allvar1["pho1_eta"],"-",allvar2["pho1_eta"], " pho2_eta:",allvar1["pho2_eta"],"-",allvar2["pho2_eta"], " mass:", allvar1["mass"], "-", allvar2["mass"]," diphomva:", allvar1["diphoMVA"], "-", allvar2["diphoMVA"], " njets:", allvar1["numJets"], "-", allvar2["numJets"], " combmva:", allvar1["combiMVA"], "-", allvar2["combiMVA"], " dijetmva:", allvar1["dijetMVA"], "-", allvar2["dijetMVA"]
        
     #   print icevt, " cat:", allvar1["cat"], "-", allvar2["ucat"],  " mass:", allvar1["mass"], "-", allvar2["mass"]," diphomva:", allvar1["diphoMVA"], "-", allvar2["diphoMVA"], " smom:", sigmaM, "-", allvar2["sigmamom"]

        


#print "********** number of events with diff = ", len(evtinTwo), "  = ", float(len(evtinTwo))/float(len(evtDump1)), " %"
#print "********** number of events affected by difference in variables : ", len(evtInDiff)
    
#d = ROOT.TCanvas ("c_mass")
#h_Mass.Draw()
#d.SetLogy()
#d.SaveAs("wrongMass.pdf")

      







        


   

    
        
    

