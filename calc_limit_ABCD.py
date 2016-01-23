from cardFileWriter import cardFileWriter
#from limit_helper import plotsignif , plotLimit , signal_bins_3fb
from math import exp,sqrt
import os,sys
import ROOT
import pickle
#import array
from Workspace.RA4Analysis.signalRegions import *

CorrSystSize = 1000.00

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

worstCaseSyst = {
  (5, 5) : {
    (250, 350) : {
      (500, -1) : 0.25
     }, 
    (350, 450) : {
      (500, -1) : 0.35
      }, 
    (450, -1) : {
      (500, -1) : 0.45
      }
    }, 
  (6, 7) : {
    (250, 350) : {
      (500, 750) : 1.00, (750, -1) : 0.30,
      (750, -1) : 1.00, (750, -1) : 0.30
     }, 
    (350, 450) : {
      (500, 750) : 1.00, (750, -1) : 0.40,
      (500, 750) : 1.00, (750, -1) : 0.55
      }, 
    (450, -1) : {
      (500, 1000) : 0.75, (1000, -1) : 0.30,
      (1000, -1) : 0.75, (1000, -1) : 0.80
      }
    }, 
  (8, -1) : {
    (250, 350) : {
      (500, 750) : 0.70, 
      (750, -1) : 0.60,
     }, 
    (350, 450) : {
      (500, -1) : 0.65
      }, 
    (450, -1) : {
      (500, -1) : 0.70
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
  result = 1.+sign*sqrt(variance)/value
  if result<0.:
    result = 0.01
  return result
  
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


#res = pickle.load(file(os.path.expandvars("singleLeptonic_Spring15__estimationResults_pkl_kappa_corrected-150116.pkl")))
sigres = pickle.load(file(os.path.expandvars("resultsFinal_withSystematics_andSignals_NewStructure_150120.pkl")))
bkgres = pickle.load(file(os.path.expandvars("resultsFinal_withSystematics_150121.pkl")))

#pdg = 'pos'
#pdg = 'neg'
#pdg = 'both'

#
# consistency
#
njetBins = sorted(bkgres.keys())
ltBins = [ ]
htBins = [ ]
for nj in njetBins:
  for lt in bkgres[nj]:
    if not lt in ltBins:
      ltBins.append(lt)
    for ht in bkgres[nj][lt]:
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
          {'color': ROOT.kBlue ,'name': 's1500' , 'mglu' : 1500, 'mlsp' : 100, 'label': 'T5q^{4} 1.5/0.8/0.1'}, \
          {'color': ROOT.kRed  ,'name': 's1200' , 'mglu' : 1200, 'mlsp' : 800, 'label': 'T5q^{4} 1.2/1.0/0.8'}, \
#          {'color': ROOT.kBlack ,'name': 's1000' , 'mglu' : 1000, 'mlsp' : 700, 'label': 'T5q^{4} 1.0/0.85/0.7'}, \
          {'color': ROOT.kBlack ,'name': 's1000' , 'mglu' : 1000, 'mlsp' : 100, 'label': 'T5q^{4} 1.0/0.55/0.1'}, \
         ]

#signal = signals[2]


#
# prepare bins
#
procNames = [ "W", "tt", "other", "QCD" ]

#nbins = 0
#for njet in njetBins[:]:
#  for lt in ltBins[:]:
#    if not lt in bkgres[njet]:
#      continue
#    for ht in htBins[:]:
#      if not ht in bkgres[njet][lt]:
#        continue
#      nbins += 1

if os.path.exists("results.log"):
  os.system("rm results.log; touch results.log")

sbBinNames = [ ]
sbBins = { }
mbBinNames = [ ]
mbBins = { }
for njet in njetBins[:1]:
  for lt in ltBins[:1]:
    if not lt in bkgres[njet]:
      continue
    for ht in htBins[:1]:
      if not ht in bkgres[njet][lt]:
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

for signal in signals[:]:
  print signal
  mglu = signal["mglu"]
  mlsp = signal["mlsp"]


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
    sbres = bkgres[sbBins[sbname][0]][sbBins[sbname][1]][sbBins[sbname][2]]
    sbsigres = sigres[sbBins[sbname][0]][sbBins[sbname][1]][sbBins[sbname][2]]
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
      if sbname.endswith("S"):  # need to include also CR yields
        c.specifyExpectation(sbname,"signal",sbsigres["signals"][mglu][mlsp]['yield_SB_W_SR'])
      else:
        c.specifyExpectation(sbname,"signal",0.)
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
      if sbname.endswith("S"):  # need to include also CR yields
        c.specifyExpectation(sbname,"signal",sbsigres["signals"][mglu][mlsp]['yield_SB_tt_SR'])
      else:
        c.specifyExpectation(sbname,"signal",0.)
      c.specifyExpectation(sbname,"W",0.001)
      c.specifyExpectation(sbname,"tt",sbres["y_crNJet_1b_"+rDPhi])
      c.specifyExpectation(sbname,"other",0.001) # others are neglected in yield
      c.specifyExpectation(sbname,"QCD",0.001)
    if sbname.endswith("C"):
      pass
#      # low DPhi
#      uncName = "yQCD" + sbname
#      c.addUncertainty(uncName,"lnN")
#      c.specifyUncertainty(uncName,sbname,"QCD",relErrForLimit(sbres["yQCD_crNJet_0b_lowDPhi"],\
#         sbres["yQCD_Var_crNJet_0b_lowDPhi"]))


  for mbname in mbBinNames:
    c.addBin(mbname,procNames,mbname)
    mbres = bkgres[mbBins[mbname][0]][mbBins[mbname][1]][mbBins[mbname][2]]
    mbsigres = sigres[mbBins[mbname][0]][mbBins[mbname][1]][mbBins[mbname][2]]
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
      c.specifyExpectation(mbname,"signal",mbsigres["signals"][mglu][mlsp]['yield_MB_CR']) # to be corrected!
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
      c.specifyExpectation(mbname,"signal",mbsigres["signals"][mglu][mlsp]['yield_MB_SR']) # to be corrected!
      c.specifyExpectation(mbname,"tt",mbres["TT_pred"])
      c.specifyExpectation(mbname,"W",mbres["W_pred"])
      c.specifyExpectation(mbname,"other",mbres["Rest_truth"])
      c.specifyExpectation(mbname,"QCD",0.001)

  #
  # global uncertainties
  #
  # c.addUncertainty("worst","lnN")
  c.addUncertainty("btag","lnN")
  c.addUncertainty("lumi","lnN")
  c.addUncertainty("sigSyst","lnN")
  for bname in sbBinNames+mbBinNames:
    c.specifyUncertainty("lumi",bname,"signal",1.046)
    c.specifyUncertainty("sigSyst",bname,"signal",1.20) # to be corrected!
    c.specifyUncertainty("lumi",bname,"other",1.046)

  for mbname in mbBinNames:
    mbnameBase = mbname[:-1]
    mbnameC = mbnameBase + "C"
    mbresC = bkgres[mbBins[mbnameC][0]][mbBins[mbnameC][1]][mbBins[mbnameC][2]]
    mbnameS = mbnameBase + "S"
    mbresS = bkgres[mbBins[mbnameS][0]][mbBins[mbnameS][1]][mbBins[mbnameS][2]]
    mbsigresS = sigres[mbBins[mbnameS][0]][mbBins[mbnameS][1]][mbBins[mbnameS][2]]

    sbWnameBase = "J3" + mbnameBase[2:]
    # sbWnameC = sbWnameBase + "C"
    # sbWresC = bkgres[sbBins[sbWnameC][0]][sbBins[sbWnameC][1]][sbBins[sbWnameC][2]]
    sbWnameS = sbWnameBase + "S"
    sbWresS = bkgres[sbBins[sbWnameS][0]][sbBins[sbWnameS][1]][sbBins[sbWnameS][2]]

    sbttnameBase = "J4" + mbnameBase[2:]
    # sbttnameC = sbttnameBase + "C"
    # sbttresC = bkgres[sbBins[sbttnameC][0]][sbBins[sbttnameC][1]][sbBins[sbttnameC][2]]
    sbttnameS = sbttnameBase + "S"
    sbttresS = bkgres[sbBins[sbttnameS][0]][sbBins[sbttnameS][1]][sbBins[sbttnameS][2]]

    if mbname.endswith("S"):

      # correlation W regions: B and F / E and F
      uncName = "corrWBF" + mbname[:-1]
      c.addUncertainty(uncName,"lnN")
      c.specifyUncertainty(uncName,"J3"+mbname[2:-1]+"S","W",CorrSystSize)
      c.specifyUncertainty(uncName,mbname,"W",CorrSystSize)
      uncName = "corrWEF" + mbname[:-1]
      c.addUncertainty(uncName,"lnN")
      c.specifyUncertainty(uncName,mbname[:-1]+"C","W",CorrSystSize)
      c.specifyUncertainty(uncName,mbname,"W",CorrSystSize)
      # correlation tt regions: D and F / E and F
      uncName = "corrTTDF" + mbname[:-1]
      c.addUncertainty(uncName,"lnN")
      c.specifyUncertainty(uncName,"J4"+mbname[2:-1]+"S","tt",CorrSystSize)
      c.specifyUncertainty(uncName,mbname,"tt",CorrSystSize)
      uncName = "corrTTEF" + mbname[:-1]
      c.addUncertainty(uncName,"lnN")
      c.specifyUncertainty(uncName,mbname[:-1]+"C","tt",CorrSystSize)
      c.specifyUncertainty(uncName,mbname,"tt",CorrSystSize)
      # anticorrelated W/tt yields from fit
      uncName = "yWtt" + mbnameC
      c.addUncertainty(uncName,"lnN")
      c.specifyUncertainty(uncName,mbnameC,"W",relErrForLimit(mbresC["yW_srNJet_0b_lowDPhi"],mbresC["yW_Var_srNJet_0b_lowDPhi"]))
      c.specifyUncertainty(uncName,mbnameC,"tt",relErrForLimit(mbresC["yTT_srNJet_0b_lowDPhi"],mbresC["yTT_Var_srNJet_0b_lowDPhi"],-1))
      # uncertainty on kappas
      uncName = "kappa" + mbnameS
      c.addUncertainty(uncName,"lnN")
      c.specifyUncertainty(uncName,mbnameS,"W",1.25) # to be corrected!
      c.specifyUncertainty(uncName,mbnameS,"tt",1.25) # to be corrected!
      # uncertainty on RCS_W (e+mu)/mu
      uncName = "rcsWemu" + mbnameS
      c.addUncertainty(uncName,"lnN")
      c.specifyUncertainty(uncName,mbnameS,"W",1.05) # to be corrected!
      # uncertainty on b tagging
      c.specifyUncertainty("btag",mbnameS,"W",1.05) # to be corrected!
      c.specifyUncertainty("btag",mbnameS,"W",1.05) # to be corrected!
      c.specifyUncertainty("btag",mbnameS,"other",1.05) # to be corrected!
      # stat. uncertainty on signal efficiency
      uncName = "statSeff" + mbnameS
      c.addUncertainty(uncName,"lnN")
      c.specifyUncertainty(uncName,mbnameS,"signal",1+mbsigresS["signals"][mglu][mlsp]["err_MB_SR"])
#      # WORST CASE SYST
#      uncName = "worst"+mbnameS
#      c.addUncertainty(uncName,"lnN")
#      c.specifyUncertainty(uncName,mbnameS,"W",1+worstCaseSyst[mbBins[mbname][0]][mbBins[mbname][1]][mbBins[mbname][2]])

    else:
      pass
#      # low DPhi
#      uncName = "yQCD" + mbname
#      c.addUncertainty(uncName,"lnN")
#      c.specifyUncertainty(uncName,mbname,"QCD",relErrForLimit(mbresC["yQCD_srNJet_0b_lowDPhi"],mbresC["yQCD_Var_srNJet_0b_lowDPhi"]))

  for sbname in sbBinNames:
    sbres = bkgres[sbBins[sbname][0]][sbBins[sbname][1]][sbBins[sbname][2]]
    # statistical uncertainty SB CRs
    obs = c.observation[sbname]
    assert obs>0
    unc = 1. + 1./sqrt(obs)
    uncName = "stat" + sbname[:-1] + "C"
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
      c.specifyUncertainty(uncName,sbname,"W",relErrForLimit(sbres["yW_crNJet_0b_highDPhi"], 
                                                             sbres["yW_Var_crNJet_0b_highDPhi"]))
      c.specifyUncertainty(uncName,sbname,"tt",relErrForLimit(sbres["yTT_crNJet_0b_highDPhi"], \
                                                                sbres["yTT_Var_crNJet_0b_highDPhi"],-1))
      if uncName=="yWttJ3L3H1D2S":
        print sbres["yW_crNJet_0b_highDPhi"],sbres["yW_Var_crNJet_0b_highDPhi"], \
            relErrForLimit(sbres["yW_crNJet_0b_highDPhi"],sbres["yW_Var_crNJet_0b_highDPhi"])
        print sbres["yTT_crNJet_0b_highDPhi"],sbres["yTT_Var_crNJet_0b_highDPhi"], \
            relErrForLimit(sbres["yTT_crNJet_0b_highDPhi"],sbres["yTT_Var_crNJet_0b_highDPhi"],-1)

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
  #
  # comments
  #
  txt = open("calc_limit.txt","a")
  txt.write("#\n")
  txt.write("# List of uncertainties\n")
  txt.write("#\n")
  txt.write("# corrWBFJxLyHzDu ... correlation W: SB/highDPhi with MB/highDPhi\n")
  txt.write("# corrWEFJxLyHzDu ... correlation W: MB/lowDPhi with MB/highDPhi\n")
  txt.write("# corrTTDFJxLyHzDu .. correlation tt: SB/highDPhi with MB/highDPhi\n")
  txt.write("# corrTTEFJxLyHzDu .. correlation tt: MB/lowDPhi with MB/highDPhi\n")
  txt.write("# yWttJxLyHzDuC ..... anti-correlated W/tt fraction fit systematics in MB CR\n")
  txt.write("# yQCDJxLyHzDuC ..... uncertainty QCD estimate in MB CR\n")
  txt.write("# statJ[34]LyHzDuC .. stat. uncertainty from yield in SB lowDPhi \n")
  txt.write("# yWttJ[34]LyHzDuS .. anti-correlated W/tt fraction fit systematics in W SB highDPhi\n")
  txt.write("# yWttJ[34]LyHzDuC ?? anti-correlated W/tt fraction fit systematics in W SB lowDPhi\n")
  txt.write("# lumi .............. luminosity\n")
  txt.write("# sigSyst ........... approximated total signal systematics")
  txt.close()

  stdout = sys.stdout
  sys.stdout = open("results.log","a")
  print 'Result ',mbBinNames[0]," , ",signal["name"]," : ",c.calcLimit(options="--run blind")
  sys.stdout.close()
  sys.stdout = stdout
