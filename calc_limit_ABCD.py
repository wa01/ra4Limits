from cardFileWriter import cardFileWriter
from limit_helper import plotsignif , plotLimit , signal_bins_3fb
from math import exp
import os,sys
import ROOT
import pickle
import array
import numpy as n
from Workspace.RA4Analysis.signalRegions import *

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
procNames = [ "W", "tt", "other" ]

sbBinNames = [ ]
sbBins = { }
mbBinNames = [ ]
mbBins = { }
for njet in njetBins[:1]:
  for lt in ltBins[:1]:
    if not lt in res[njet]:
      continue
    for ht in htBins[:1]:
      if not ht in res[njet][lt]:
        continue
      dphiLimit = dphiLimitToLabel(regionToDPhi[njet][lt][ht])
      bNameBase = njetBinToLabel(njet) + ltBinToLabel(lt) + htBinToLabel(ht) + dphiLimit
      for r in [ "CR", "SR" ]:
        bName = bNameBase + r
        assert not bName in mbBinNames
        mbBinNames.append(bName)
        mbBins[bName] = ( njet, lt, ht )
        for sb in [ "W", "tt" ]:
          if r=="CR":
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
  c.defWidth=12
  c.precision=6
  #
  # define bins
  #
  for sbname in sbBinNames:
    c.addBin(sbname,procNames,sbname)
    sbres = res[sbBins[sbname][0]][sbBins[sbname][1]][sbBins[sbname][2]]
    # low vs. high dphi
    r = sbname[-2:]
    rDPhi = "low" if r=="CR" else "high"
    rDPhi += "DPhi"
    # observation
    y_truth = sbres["yW_crNJet_0b_"+rDPhi+"_truth"] + \
        sbres["ytt_crNJet_0b_"+rDPhi+"_truth"] + \
        sbres["yRest_crNJet_0b_"+rDPhi+"_truth"]
    c.specifyObservation(sbname,int(y_truth+0.5))
    # W
    if sbname[:2]=="J3":
      pass
    # tt
    elif sbname[:2]=="J4":
      c.specifyExpectation(sbname,"W",sbres["yW_crNJet_1b_"+rDPhi+"_truth"])
      c.specifyExpectation(sbname,"W",sbres["yW_crNJet_1b_"+rDPhi+"_truth"])
      

  for mbname in mbBinNames:
    c.addBin(mbname,procNames,mbname)

sys.exit(0)
#          bNameBase = "J"+str(njet(0)
#        c.addBin

for signal in signals[:1]:
  for njet in njetBins[:1]:
    for lt in ltBins[:1]:
      for ht in htBins[:1]:


        found_bin.append(\
        { #'closure': res[njet][lt][ht]['tot_clos']  ,
            'Berror':res[njet][lt][ht]['tot_pred_err'] ,'B': res[njet][lt][ht]['tot_pred'], 'nJet': (6, 7), 'S1000': res[njet][lt][ht]['T5q^{4} 1.0/0.8/0.7_yield'], 'HT': ht, 'LT': lt, 'S1500': res[njet][lt][ht]['T5q^{4} 1.5/0.8/0.1_yield'], 'deltaPhi': 0.75, 'S1200': res[njet][lt][ht]['T5q^{4} 1.2/1.0/0.8_yield']},\
                        )


  print found_bin

  x_s = n.zeros(11, dtype=float)
  y_s = n.zeros(11, dtype=float)
  y_1min = n.zeros(11, dtype=float)
  y_2min = n.zeros(11, dtype=float)
  y_1max = n.zeros(11, dtype=float)
  y_2max = n.zeros(11, dtype=float)
  for lum in lumi_bins:
    print "lum:" , lum
    print "lumi:" , lum*1000
    #c.addUncertainty('Lumi', 'lnN')
    #c.specifyFlatUncertainty('Lumi', 1.20)
    c.addUncertainty('JES', 'lnN')
    #c.addUncertainty('closure', 'lnN')
    c.addUncertainty('predUnc', 'lnN')
    c.addUncertainty('SigAcc', 'lnN')
    for i , bin in enumerate(found_bin):
      #print i ,bin['HT'],bin['LT'] ,bin['nJet'],bin['B'] , bin['S1500']

      bkg_Y = { 'name': 'bin_'+str(i), 'value': float(bin['B'])/float(lumi_origin), 'label': 'ttJets+WJets'}

      signal.update({'value': float(bin[signal['name']])/float(lumi_origin)})
      #print signal

      c.addBin(bkg_Y['name'], ['bkg'], bkg_Y['name'])


      #print "Bin"+str(i), y
      c.specifyObservation(bkg_Y['name'], int(signal['value']*lum))
      c.specifyExpectation(bkg_Y['name'], 'bkg', bkg_Y['value']*lum)
      c.specifyExpectation(bkg_Y['name'], 'signal', signal['value']*lum)

      c.specifyUncertainty('JES', bkg_Y['name'], 'bkg', 1.2)
      #c.specifyUncertainty('closure', bkg_Y['name'], 'bkg', 1+bin['closure'])
      c.specifyUncertainty('predUnc', bkg_Y['name'], 'bkg', 1+bin['Berror'])
      c.specifyUncertainty('SigAcc', bkg_Y['name'], 'signal', 1.2)
    #####End of Bin loop######


    if option == 'limit' :

      limit = c.calcLimit()
      print 'limit:' ,  limit
      limit_med = limit['0.500']
      print "limit median :" , limit_med
      y_2min[int(lum)] = limit['0.025']
      y_1min[int(lum)] = limit['0.160']
      y_s[int(lum)]    = limit_med
      y_1max[int(lum)] = limit['0.840']
      y_2max[int(lum)] = limit['0.975']
      x_s[int(lum)] = int(lum)


  signal_signif.append({'label':signal['label'] ,'name':signal['name'],'color':signal['color'] , 'y_m':y_s[1:] , 'x': x_s[1:] , 'y1_min':y_1min[1:],'y2_min': y_2min[1:] ,'y1_max': y_1max[1:] ,'y2_max':y_2max[1:] })

print signal_signif

for signal in signal_signif:
   plotLimit(signal,path,option+pdg+'_fullbin_'+str(lumi_origin),lumi_origin)


