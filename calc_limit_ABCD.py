from cardFileWriter import cardFileWriter
from limit_helper import plotsignif , plotLimit , signal_bins_3fb
from math import exp,sqrt
import os,sys
import ROOT
import pickle
import array
import numpy as n
from Workspace.RA4Analysis.signalRegions import *

print "BLA"

regionToDPhi = {
  (5, 5) : {
    (250, 350) : {
      (500, -1) : 1.00
     }, 
    (350, 450) : {
      (500, -1) : 1.00
      }, 
    (450, -1) : {
      (500, -1) : 1.00
      }
    }, 
  (6, 7) : {
    (250, 350) : {
      (500, 750) : 1.00, (750, -1) : 1.00
     }, 
    (350, 450) : {
      (500, 750) : 1.00, (750, -1) : 1.00
      }, 
    (450, -1) : {
      (500, 1000) : 0.75, (1000, -1) : 0.75
      }
    }, 
  (8, -1) : {
    (250, 350) : {
      (500, 750) : 1.00, (750, -1) : 1.00
     }, 
    (350, 450) : {
      (500, -1) : 0.75
      }, 
    (450, -1) : {
      (500, -1) : 0.75
      }
    }

}

def dphiLimitToLabel(dphi):
  ndphi = int(100*dphi+0.5)
  result = None
  if ndphi==75:
    result = "1"
  elif ndphi==100:
    result = "2"
  assert result!=None
  return "D"+result

def njetBinToLabel(njBin):
  # simplified to lower boundary
  result = None
  if njBin[1]!=-1:
    result = "".join([str(x) for x in range(njBin[0],njBin[1]+1)])
  else:
    result = str(njBin[0])+"p"
  assert result!=None
  return "J"+result[0]
#  return "J"+result

def ltBinToLabel(ltBin):
  idxs = [ None, None ]
  for i in range(2):
    if ltBin[i]==250:
      idxs[i] = 1
    elif ltBin[i]==350:
      idxs[i] = 2
    elif ltBin[i]==450:
      idxs[i] = 3
    elif ltBin[i]==-1:
      idxs[i] = 4
    assert idxs[i]!=None
  return "L"+"".join([str(x) for x in range(idxs[0],idxs[1])])

def htBinToLabel(htBin):
  # simplified to lower limit
  idxs = [ None, None ]
  for i in range(2):
    if htBin[i]==500:
      idxs[i] = 1
    elif htBin[i]==750:
      idxs[i] = 2
    elif htBin[i]==1000:
      idxs[i] = 3
    elif htBin[i]==-1:
      idxs[i] = 4
    assert idxs[i]!=None
  return "H"+str(idxs[0])
#  return "H"+"".join([str(x) for x in range(idxs[0],idxs[1])])

def relErrForLimit(value,variance,sign=1):
  return 1.+sign*sqrt(variance)/value
  
ROOT.gROOT.LoadMacro("$CMSSW_BASE/src/Workspace/HEPHYPythonTools/scripts/root/tdrstyle.C")
ROOT.setTDRStyle()

path = os.environ["HOME"]+"/www/combine_tests/"
if not os.path.exists(path):
  os.makedirs(path)

path_table = os.environ["HOME"]+"/www/combine_tests/"
if not os.path.exists(path_table):
  os.makedirs(path_table)

text_path = "text_files"
if not os.path.exists(text_path):
  os.makedirs(text_path)

options = ['signif' , 'limit']
option = options[1]

lumi_bins = [1,2,3,4,5,6,7,8,9,10]
#lumi_bins = [1,2,3,4]
#lumi_bins = [3,10]
lumi_origin = 3


#bin_yields2 = pickle.load(file('/data/dspitzbart/lumi3.0yields_pkl_newOpt'))
#bin_yields2 = pickle.load(file('/data/dspitzbart/lumi3.0yields_pkl_final'))
#res = pickle.load(file('/data/easilar/PHYS14v3/withCSV/rCS_0b_3.0fbSlidingWcorrectionMuonChannel/singleLeptonic_Phys14V3__estimationResults_pkl_updated'))
#res = pickle.load(file(os.path.expandvars("$WORK/susy/DataDspitzbartResults2015/Prediction_SFtemplates_fullSR_lep_data_2.1/resultsFinal_withSystematics_pkl")))
res = pickle.load(file(os.path.expandvars("singleLeptonic_Spring15__estimationResults_pkl_kappa_corrected-150116.pkl")))

#pdg = 'pos'
#pdg = 'neg'
#pdg = 'both'

#
# consistency
#
njetBins = sorted(res.keys())
ltBins = [ ]
htBins = [ ]
for nj in njetBins:
  for lt in res[nj]:
    if not lt in ltBins:
      ltBins.append(lt)
    for ht in res[nj][lt]:
      if not ht in htBins:
        htBins.append(ht)
ltBins.sort()
htBins.sort()
print njetBins
print ltBins
print htBins
for njet in njetBins:
  print njetBinToLabel(njet)
for lt in ltBins:
  print ltBinToLabel(lt)
for ht in htBins:
  print htBinToLabel(ht)

signals = [
          {'color': ROOT.kBlue ,'name': 'S1500' ,'label': 'T5q^{4} 1.5/0.8/0.1'}, \
          {'color': ROOT.kRed  ,'name': 'S1200' ,'label': 'T5q^{4} 1.2/1.0/0.8'}, \
          {'color': ROOT.kBlack ,'name': 'S1000' ,'label': 'T5q^{4} 1.0/0.8/0.7'}, \
         ]

#signal = signals[2]


#
# prepare bins
#
procNames = [ "W", "tt", "other", "QCD" ]

sbBinNames = [ ]
sbBins = { }
mbBinNames = [ ]
mbBins = { }
for njet in njetBins[:]:
  for lt in ltBins[:]:
    if not lt in res[njet]:
      continue
    for ht in htBins[:]:
      if not ht in res[njet][lt]:
        continue
      dphiLimit = dphiLimitToLabel(regionToDPhi[njet][lt][ht])
      bNameBase = njetBinToLabel(njet) + ltBinToLabel(lt) + htBinToLabel(ht) + dphiLimit
      for r in [ "C", "S" ]:
        bName = bNameBase + r
        assert not bName in mbBinNames
        mbBinNames.append(bName)
        mbBins[bName] = ( njet, lt, ht )
        for sb in [ "W", "tt" ]:
          if r=="C":
            continue
          if sb=="W":
            sbName = "J3"  + ltBinToLabel(lt) + htBinToLabel(ht) + dphiLimit + r
            if not sbName in sbBinNames:
              sbBinNames.append(sbName)
              sbBins[sbName] = ( njet, lt, ht )
          elif sb=="tt":
            sbName = "J4"  + ltBinToLabel(lt) + htBinToLabel(ht) + dphiLimit + r
            if not sbName in sbBinNames:
              sbBinNames.append(sbName)
              sbBins[sbName] = ( njet, lt, ht )            

print mbBinNames
print sbBinNames                

signal_signif = []
for signal in signals[:1]:
  print signal


  c = cardFileWriter()
  c.defWidth=10
  c.precision=3
  c.maxUncNameWidth = 17
  c.maxUncStrWidth = 15
  #
  # define bins
  #
  for sbname in sbBinNames:
    c.addBin(sbname,procNames,sbname)
    sbres = res[sbBins[sbname][0]][sbBins[sbname][1]][sbBins[sbname][2]]
    # low vs. high dphi
    r = sbname[-2:]
    rDPhi = "low" if r=="C" else "high"
    rDPhi += "DPhi"
    # calculate missing numbers
    # yield (W) = observed - estimated tt - EWK(other)
    # wYield = sbres["y_crNJet_0b_"+rDPhi] - sbres["yTT_crNJet_0b_"+rDPhi] - \
    #    sbres["yRest_crNJet_0b_"+rDPhi+"_truth"]
    wYield = sbres["y_crNJet_0b_"+rDPhi] - sbres["yTT_crNJet_0b_"+rDPhi] # others are neglected in yield
    sbres["yW_crNJet_0b_"+rDPhi] = wYield
    # error on wYield
    # wVar = sbres["y_Var_crNJet_0b_"+rDPhi] + sbres["yTT_Var_crNJet_0b_"+rDPhi] + \
    #    sbres["yRest_Var_crNJet_0b_"+rDPhi+"_truth"]
    wVar = sbres["y_Var_crNJet_0b_"+rDPhi] + sbres["yTT_Var_crNJet_0b_"+rDPhi] # others are neglected in yield
    sbres["yW_Var_crNJet_0b_"+rDPhi] = wVar
    
    # W
    if sbname[:2]=="J3":
      # observation
      y_truth = sbres["yW_crNJet_0b_"+rDPhi+"_truth"] + \
          sbres["yTT_crNJet_0b_"+rDPhi+"_truth"] + \
          sbres["yRest_crNJet_0b_"+rDPhi+"_truth"]
      c.specifyObservation(sbname,int(y_truth+0.5))
      c.specifyExpectation(sbname,"signal",0.001)
#      c.specifyExpectation(sbname,"W",sbres["y_crNJet_0b_"+rDPhi]-sbres["yTT_crNJet_0b_"+rDPhi])
      c.specifyExpectation(sbname,"W",sbres["yW_crNJet_0b_"+rDPhi])
      c.specifyExpectation(sbname,"tt",sbres["yTT_crNJet_0b_"+rDPhi])
      c.specifyExpectation(sbname,"other",0.001) # others are neglected in yield
      c.specifyExpectation(sbname,"QCD",0.001)
    # tt
    elif sbname[:2]=="J4":
      # observation
      y_truth = sbres["yW_crNJet_1b_"+rDPhi+"_truth"] + \
          sbres["yTT_crNJet_1b_"+rDPhi+"_truth"] + \
          sbres["yRest_crNJet_1b_"+rDPhi+"_truth"]
      c.specifyObservation(sbname,int(y_truth+0.5))
      c.specifyExpectation(sbname,"signal",0.001)
      c.specifyExpectation(sbname,"W",0.001)
      c.specifyExpectation(sbname,"tt",sbres["y_crNJet_1b_"+rDPhi])
      c.specifyExpectation(sbname,"other",0.001) # others are neglected in yield
      c.specifyExpectation(sbname,"QCD",0.001)
    if sbname.endswith("C"):
      pass
#      # low DPhi
#      uncName = "yQCD" + sbname
#      c.addUncertainty(uncName,"lnN")
#      c.specifyUncertainty(uncName,sbname,"QCD",relErrForLimit(sbres["yQCD_crNJet_0b_lowDPhi"],sbres["yQCD_Var_crNJet_0b_lowDPhi"]))

      

  for mbname in mbBinNames:
    c.addBin(mbname,procNames,mbname)
    mbres = res[mbBins[mbname][0]][mbBins[mbname][1]][mbBins[mbname][2]]
    # low vs. high dphi
    r = mbname[-2:]
    if r=="C":
      rDPhi = "lowDPhi"
      # observation
      y_truth = mbres["yW_srNJet_0b_"+rDPhi+"_truth"] + \
          mbres["yTT_srNJet_0b_"+rDPhi+"_truth"] + \
          mbres["yRest_srNJet_0b_"+rDPhi+"_truth"] + \
          mbres["yQCD_srNJet_0b_"+rDPhi+"_truth"]
      c.specifyObservation(mbname,int(y_truth+0.5))
      # expectation
      c.specifyExpectation(mbname,"signal",0.01)
      c.specifyExpectation(mbname,"tt",mbres["yTT_srNJet_0b_"+rDPhi])
      c.specifyExpectation(mbname,"W",mbres["yW_srNJet_0b_"+rDPhi])
      c.specifyExpectation(mbname,"other",mbres["yRest_srNJet_0b_"+rDPhi+"_truth"])
      # !!!!! to be changed
      c.specifyExpectation(mbname,"QCD",0.001)
      # c.specifyExpectation(mbname,"QCD",mbres["yQCD_srNJet_0b"])
    else:
      rDPhi = "highDPhi"
      # observation
      y_truth = mbres["W_truth"] +  mbres["TT_truth"] + mbres["Rest_truth"]
      c.specifyObservation(mbname,int(y_truth+0.5))
      # expectation
      c.specifyExpectation(mbname,"signal",1.)
      c.specifyExpectation(mbname,"tt",mbres["TT_pred"])
      c.specifyExpectation(mbname,"W",mbres["W_pred"])
      c.specifyExpectation(mbname,"other",mbres["Rest_truth"])
      c.specifyExpectation(mbname,"QCD",0.001)

  for mbname in mbBinNames:
    mbnameBase = mbname[:-1]
    mbnameC = mbnameBase + "C"
    mbresC = res[mbBins[mbnameC][0]][mbBins[mbnameC][1]][mbBins[mbnameC][2]]
    mbnameS = mbnameBase + "S"
    mbresS = res[mbBins[mbnameS][0]][mbBins[mbnameS][1]][mbBins[mbnameS][2]]

    sbWnameBase = "J3" + mbnameBase[2:]
    # sbWnameC = sbWnameBase + "C"
    # sbWresC = res[sbBins[sbWnameC][0]][sbBins[sbWnameC][1]][sbBins[sbWnameC][2]]
    sbWnameS = sbWnameBase + "S"
    sbWresS = res[sbBins[sbWnameS][0]][sbBins[sbWnameS][1]][sbBins[sbWnameS][2]]

    sbttnameBase = "J4" + mbnameBase[2:]
    # sbttnameC = sbttnameBase + "C"
    # sbttresC = res[sbBins[sbttnameC][0]][sbBins[sbttnameC][1]][sbBins[sbttnameC][2]]
    sbttnameS = sbttnameBase + "S"
    sbttresS = res[sbBins[sbttnameS][0]][sbBins[sbttnameS][1]][sbBins[sbttnameS][2]]

    if mbname.endswith("S"):

      # correlation W regions: B and F / E and F
      uncName = "corrWBF" + mbname[:-1]
      c.addUncertainty(uncName,"lnN")
      c.specifyUncertainty(uncName,"J3"+mbname[2:-1]+"S","W",3.00)
      c.specifyUncertainty(uncName,mbname,"W",3.00)
      uncName = "corrWEF" + mbname[:-1]
      c.addUncertainty(uncName,"lnN")
      c.specifyUncertainty(uncName,mbname[:-1]+"C","W",3.00)
      c.specifyUncertainty(uncName,mbname,"W",3.00)
      # correlation tt regions: D and F / E and F
      uncName = "corrTTDF" + mbname[:-1]
      c.addUncertainty(uncName,"lnN")
      c.specifyUncertainty(uncName,"J4"+mbname[2:-1]+"S","tt",3.00)
      c.specifyUncertainty(uncName,mbname,"tt",3.00)
      uncName = "corrTTEF" + mbname[:-1]
      c.addUncertainty(uncName,"lnN")
      c.specifyUncertainty(uncName,mbname[:-1]+"C","tt",3.00)
      c.specifyUncertainty(uncName,mbname,"tt",3.00)
      # anticorrelated W/tt yields from fit
      uncName = "yWtt" + mbnameC
      c.addUncertainty(uncName,"lnN")
      c.specifyUncertainty(uncName,mbnameC,"W",relErrForLimit(mbresC["yW_srNJet_0b_lowDPhi"],mbresC["yW_Var_srNJet_0b_lowDPhi"]))
      c.specifyUncertainty(uncName,mbnameC,"tt",relErrForLimit(mbresC["yTT_srNJet_0b_lowDPhi"],mbresC["yTT_Var_srNJet_0b_lowDPhi"],-1))

    else:
      pass
#      # low DPhi
#      uncName = "yQCD" + mbname
#      c.addUncertainty(uncName,"lnN")
#      c.specifyUncertainty(uncName,mbname,"QCD",relErrForLimit(mbresC["yQCD_srNJet_0b_lowDPhi"],mbresC["yQCD_Var_srNJet_0b_lowDPhi"]))

  for sbname in sbBinNames:
    sbres = res[sbBins[sbname][0]][sbBins[sbname][1]][sbBins[sbname][2]]
    # statistical uncertainty SB CRs
    obs = c.observation[sbname]
    assert obs>0
    unc = 1. + 1./sqrt(obs)
    uncName = "stat" + sbname
    c.addUncertainty(uncName,"lnN")
    for mbname in mbBinNames:
      if not mbname.endswith("S"):
        continue
      if "J3"+mbname[2:-1]+"C"==sbname:
        c.specifyUncertainty(uncName,mbname,"W",unc)
      elif "J4"+mbname[2:-1]+"C"==sbname:
        c.specifyUncertainty(uncName,mbname,"tt",unc)
    # anticorrelated W/tt yields from fit
    if sbname.startswith("J3"):
      uncName = "yWtt" + sbname
      c.addUncertainty(uncName,"lnN")
      c.specifyUncertainty(uncName,sbname,"W",relErrForLimit(sbres["yW_crNJet_0b_highDPhi"],sbres["yW_Var_crNJet_0b_highDPhi"]))
      c.specifyUncertainty(uncName,sbname,"tt",relErrForLimit(sbres["yTT_crNJet_0b_highDPhi"],sbres["yTT_Var_crNJet_0b_highDPhi"],-1))

  #
  # global normalization
  #
  c.addUncertainty("lumi","lnN")
  c.addUncertainty("sigSyst","lnN")
  for bname in sbBinNames+mbBinNames:
    c.specifyUncertainty("lumi",bname,"signal",1.046)
    c.specifyUncertainty("sigSyst",bname,"signal",1.20)
    c.specifyUncertainty("lumi",bname,"other",1.046)
  #
  # QCD
  #
#  for bname in mbBinNames:
#    if bname.endswith("C"):
#      c.addUncertainty(uncName,"lnN")
#      c.specifyUncertainty(uncName,bname,"QCD",relErrForLimit()
#    c.specifyUncertainty("sigSyst",bname,"signal",1.20)
#    c.specifyUncertainty("lumi",bname,"other",1.046)

  c.writeToFile("calc_limit.txt")

