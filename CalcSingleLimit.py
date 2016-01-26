from cardFileWriter import cardFileWriter
#from limit_helper import plotsignif , plotLimit , signal_bins_3fb
from math import exp,sqrt,isnan
import os,sys
import ROOT

def relErrForLimit(value,variance,sign=1):
    result = 1.+sign*sqrt(variance)/value
    if result<0.:
        result = 0.01
    return result

def relErrorsOnFractions(yields,varYields):
    ny = len(yields)
    assert ny==len(varYields)
    sy = sum(yields)
    a = ROOT.TMatrixD(ny,ny)
    c = ROOT.TMatrixDSym(ny)
    for i in range(ny):
        rai = ROOT.TMatrixDRow(a,i)
        rci = ROOT.TMatrixDRow(c,i)
        for j in range(ny):
            if j==i:
                rci[j] = varYields[i]
            else:
                rci[j] = 0.
            raj = ROOT.TMatrixDRow(a,i)
            v = -yields[i]
            if j==i:
                v += sy
            v /= sy*sy
            rai[j] = v
    d = c.Similarity(a)
    result = [ ]
    for i in range(ny):
        result.append(sqrt(d[i][i])/(yields[i]/sy))
    return result
#    for i in range(ny):
#        line = ""
#        for j in range(ny):
#            if j==i:
#                line += "{0:8.4f}".format(sqrt(d[i][j]))
#            else:
#                line += "{0:8.4f}".format(d[i][j]/sqrt(d[i][i]*d[j][j]))
#        print line
            
    
class CalcSingleLimit:

    def __init__(self,bkgres,sbBinNames,sbBins,mbBinNames,mbBins,sigres,signal):

        self.name = "calc_single_limit"
        self.runLimit = False
        self.runBlind = False
        self.useBins = range(len(mbBinNames))
        self.corrSystSize = 100.0
        self.procNames = [ "W", "tt", "other", "QCD" ]

        self.bkgres = bkgres
        self.sbBinNames = sbBinNames
        self.sbBins = sbBins
        self.mbBinNames = mbBinNames
        self.mbBins = mbBins

        self.sigres = sigres
        self.signal = signal
        self.mglu = signal["mglu"]
        self.mlsp = signal["mlsp"]

        self.c = cardFileWriter()
        self.c.defWidth=10
        self.c.precision=3
        self.c.maxUncNameWidth = 17
        self.c.maxUncStrWidth = 15

    def subDict(self,d,bins):
        return d[bins[0]][bins[1]][bins[2]]

    def sigSubDict(self,d):
        return d["signals"][self.mglu][self.mlsp]

    def limitSinglePoint(self):

      mbBinNames = [ ]
      for i,n in enumerate(self.mbBinNames):
          if i in self.useBins:
              mbBinNames.append(n)
      sbBinNames = [ ]
      for n in self.sbBinNames:
          for m in mbBinNames:
              if n[2:]==m[2:]:
                  sbBinNames.append(n)
                  break
      print "Using mb bins ",mbBinNames
      print "Using sb bins ",sbBinNames
      #
      # bin definition; observed and expected counts
      #
      for sbname in sbBinNames:
        sbnameS = sbname + "S"
        self.c.addBin(sbnameS,self.procNames,sbnameS)
        sbres = self.subDict(self.bkgres,self.sbBins[sbname])
        self.sbsigres = self.subDict(self.sigres,self.sbBins[sbname])

        #
        # debug
        #
        #if sbname.startswith("J4"):
        #    numerator = sbres["y_crNJet_1b_highDPhi"] - sbres["yQCD_crNJet_1b_highDPhi"]
        #    denominator = sbres["y_crNJet_1b_lowDPhi"] - sbres["yQCD_crNJet_1b_lowDPhi"]
        #    print "ttbar SB :",numerator,denominator,numerator/denominator,sbres["rCS_crLowNJet_1b"]
        #    print sbres["rCS_crLowNJet_1b"]["rCS"],sbres["rCS_crLowNJet_1b_kappa"]["rCS"], \
        #        sbres["rCS_crLowNJet_1b_onlyTT"]["rCS"],sbres["rCS_srNJet_0b_onlyTT"]["rCS"]
        #


        r = "S"
        rDPhi = "low" if r=="C" else "high"
        rDPhi += "DPhi"
        # calculate missing numbers for both low and high dPhi
        # yield (W) = observed - estimated tt # others are neglected in yield
        wYield = sbres["y_crNJet_0b_lowDPhi"] - sbres["yTT_crNJet_0b_lowDPhi"] - \
            sbres["yRest_crNJet_0b_lowDPhi_truth"]
        sbres["yW_crNJet_0b_lowDPhi"] = wYield
        wYield = sbres["y_crNJet_0b_highDPhi"] - sbres["yTT_crNJet_0b_highDPhi"] - \
            sbres["yRest_crNJet_0b_highDPhi_truth"]
        sbres["yW_crNJet_0b_highDPhi"] = wYield
        # error on wYield
        # yield (W) = observed - estimated tt # others are neglected in yield
        wVar = sbres["y_Var_crNJet_0b_lowDPhi"] + sbres["yTT_Var_crNJet_0b_lowDPhi"] + sbres["yRest_Var_crNJet_0b_lowDPhi_truth"]
        sbres["yW_Var_crNJet_0b_lowDPhi"] = wVar
        wVar = sbres["y_Var_crNJet_0b_highDPhi"] + sbres["yTT_Var_crNJet_0b_highDPhi"] + sbres["yRest_Var_crNJet_0b_highDPhi_truth"]
        sbres["yW_Var_crNJet_0b_highDPhi"] = wVar
        #
        # define W sideband
        #
        if sbnameS[:2]=="J3":
          # observation
          # currently: derive from truth (!*! assume no QCD in SB/SR)
          #y_truth = sbres["yW_crNJet_0b_"+rDPhi+"_truth"] + \
          #    sbres["yTT_crNJet_0b_"+rDPhi+"_truth"] + \
          #    sbres["yRest_crNJet_0b_"+rDPhi+"_truth"]
          self.c.specifyObservation(sbnameS,int(sbres["y_crNJet_0b_highDPhi"]+0.5))
          self.c.specifyExpectation(sbnameS,"signal",self.sigSubDict(self.sbsigres)['yield_SB_W_SR'])
          self.c.specifyExpectation(sbnameS,"W",sbres["yW_crNJet_0b_"+rDPhi])
          self.c.specifyExpectation(sbnameS,"tt",sbres["yTT_crNJet_0b_"+rDPhi])
          self.c.specifyExpectation(sbnameS,"other",sbres["yRest_crNJet_0b_"+rDPhi+"_truth"])
          self.c.specifyExpectation(sbnameS,"QCD",0.001) # QCD is neglected in yield
        #
        # define tt sideband
        # 
        elif sbnameS[:2]=="J4":
          # observation
          #y_truth = sbres["yW_crNJet_1b_"+rDPhi+"_truth"] + \
          #    sbres["yTT_crNJet_1b_"+rDPhi+"_truth"] + \
          #    sbres["yRest_crNJet_1b_"+rDPhi+"_truth"]
          self.c.specifyObservation(sbnameS,int(sbres["y_crNJet_1b_highDPhi"]+0.5))
          self.c.specifyExpectation(sbnameS,"signal",self.sigSubDict(self.sbsigres)['yield_SB_tt_SR'])
          self.c.specifyExpectation(sbnameS,"W",0.001)
          self.c.specifyExpectation(sbnameS,"tt",sbres["y_crNJet_1b_"+rDPhi])
          self.c.specifyExpectation(sbnameS,"other",0.001) 
#          self.c.specifyExpectation(sbnameS,"other",sbres["yRest_crNJet_1b_"+rDPhi+"_truth"]) 
          self.c.specifyExpectation(sbnameS,"QCD",sbres["yQCD_crNJet_1b_highDPhi"])

      for mbname in mbBinNames:
        mbres = self.subDict(self.bkgres,self.mbBins[mbname])
        self.mbsigres = self.subDict(self.sigres,self.mbBins[mbname])
        #
        # low dPhi (CR)
        #
        mbnameC = mbname + "C"
        self.c.addBin(mbnameC,self.procNames,mbnameC)
        rDPhi = "lowDPhi"
        # observation
        #y_truth = mbres["yW_srNJet_0b_"+rDPhi+"_truth"] + \
        #    mbres["yTT_srNJet_0b_"+rDPhi+"_truth"] + \
        #    mbres["yRest_srNJet_0b_"+rDPhi+"_truth"] + \
        #    mbres["yQCD_srNJet_0b_"+rDPhi+"_truth"]
        self.c.specifyObservation(mbnameC,int(mbres["y_srNJet_0b_lowDPhi"]+0.5))
        # expectation
        self.c.specifyExpectation(mbnameC,"signal",self.sigSubDict(self.mbsigres)['yield_MB_CR'])
        self.c.specifyExpectation(mbnameC,"tt",mbres["yTT_srNJet_0b_"+rDPhi])
        self.c.specifyExpectation(mbnameC,"W",mbres["yW_srNJet_0b_"+rDPhi])
        self.c.specifyExpectation(mbnameC,"other",mbres["yRest_srNJet_0b_"+rDPhi+"_truth"])
        self.c.specifyExpectation(mbnameC,"QCD",mbres["yQCD_srNJet_0b_"+rDPhi])
        #
        # high dPhi
        #
        mbnameS = mbname + "S"
        self.c.addBin(mbnameS,self.procNames,mbnameS)
        rDPhi = "highDPhi"
        # observation
        #y_truth = mbres["W_truth"] +  mbres["TT_truth"] + mbres["Rest_truth"]
        self.c.specifyObservation(mbnameS,int(mbres["y_srNJet_0b_highDPhi"]+0.5))
        # expectation
        self.c.specifyExpectation(mbnameS,"signal",self.sigSubDict(self.mbsigres)['yield_MB_SR'])
        self.c.specifyExpectation(mbnameS,"tt",mbres["TT_pred_final"])
        self.c.specifyExpectation(mbnameS,"W",mbres["W_pred_final"])
        self.c.specifyExpectation(mbnameS,"other",mbres["Rest_truth"])
        self.c.specifyExpectation(mbnameS,"QCD",0.001)

      #
      # global uncertainties
      #
      # self.c.addUncertainty("worst","lnN")
      self.c.addUncertainty("btag","lnN")
      self.c.addUncertainty("lumi","lnN")
      self.c.addUncertainty("sigSyst","lnN")
      self.c.addUncertainty("xsecOther","lnN")
      for bname in sbBinNames:
        sbname = bname + "S"
        self.c.specifyUncertainty("lumi",sbname,"signal",1.046)
        self.c.specifyUncertainty("sigSyst",sbname,"signal",1.20) # to be corrected!
        self.c.specifyUncertainty("lumi",sbname,"other",1.046)
        self.c.specifyUncertainty("xsecOther",sbname,"other",1.50)
      for bname in mbBinNames:
        for r in [ "C", "S" ]:
          mbname = bname + r
          self.c.specifyUncertainty("lumi",mbname,"signal",1.046)
          self.c.specifyUncertainty("sigSyst",mbname,"signal",1.20) # to be corrected!
          self.c.specifyUncertainty("lumi",mbname,"other",1.046)
          self.c.specifyUncertainty("xsecOther",sbname,"other",1.50)
      #
      # correlations between MB/SR and MB/CR or SB/SR
      #
      for mbname in mbBinNames:
        bname = mbname[2:]
        mbnameC = mbname + "C"
        mbnameS = mbname + "S"
        mbres = self.subDict(self.bkgres,self.mbBins[mbname])
        self.mbsigres = self.subDict(self.sigres,self.mbBins[mbname])

        sbWname = "J3" + bname
        sbWnameS = sbWname + "S"
        sbWresS = self.subDict(self.bkgres,self.sbBins[sbWname])

        sbttname = "J4" + bname
        sbttnameS = sbttname + "S"
        sbttresS = self.subDict(self.bkgres,self.sbBins[sbttname])
        #
        # correlation W regions: B and F / E and F
        #
        uncName = "corrWBF" + mbname
        self.c.addUncertainty(uncName,"lnU")
        self.c.specifyUncertainty(uncName,"J3"+bname+"S","W",self.corrSystSize)
        self.c.specifyUncertainty(uncName,mbnameS,"W",self.corrSystSize)
        uncName = "corrWEF" + mbname
        self.c.addUncertainty(uncName,"lnU")
        self.c.specifyUncertainty(uncName,mbnameC,"W",self.corrSystSize)
        self.c.specifyUncertainty(uncName,mbnameS,"W",self.corrSystSize)
        #
        # correlation tt regions: D and F / E and F
        #
        uncName = "corrTTDF" + mbname
        self.c.addUncertainty(uncName,"lnU")
        self.c.specifyUncertainty(uncName,"J4"+bname+"S","tt",self.corrSystSize)
        self.c.specifyUncertainty(uncName,mbnameS,"tt",self.corrSystSize)
        uncName = "corrTTEF" + mbname
        self.c.addUncertainty(uncName,"lnU")
        self.c.specifyUncertainty(uncName,mbnameC,"tt",self.corrSystSize)
        self.c.specifyUncertainty(uncName,mbnameS,"tt",self.corrSystSize)
      #
      # (anti-)correlations from fitted W/tt yields
      #
      for mbname in mbBinNames:
        bname = mbname[2:]
        mbnameC = mbname + "C"
        mbnameS = mbname + "S"
        mbres = self.subDict(self.bkgres,self.mbBins[mbname])
        self.mbsigres = self.subDict(self.sigres,self.mbBins[mbname])

        sbWname = "J3" + bname
        sbWnameS = sbWname + "S"
        uncName = "yWtt" + sbWnameS
        if not uncName in self.c.uncertainties:
            #
            # anticorrelated W/tt yields from fit (translated from low dphi region)
            #   use errors on fraction (total normalization is included in Poisson error of bin)
            #
            sbWres = self.subDict(self.bkgres,self.sbBins[sbWname])
            ys = [ sbWres["yW_crNJet_0b_lowDPhi"], sbWres["yTT_crNJet_0b_lowDPhi"], sbWres["yRest_crNJet_0b_lowDPhi_truth"] ]
            vys = [ sbWres["yW_Var_crNJet_0b_lowDPhi"], sbWres["yTT_Var_crNJet_0b_lowDPhi"], sbWres["yRest_Var_crNJet_0b_lowDPhi_truth"] ]
            fracErrs = relErrorsOnFractions(ys,vys)
            self.c.addUncertainty(uncName,"lnN")
            self.c.specifyUncertainty(uncName,sbWnameS,"W",1.+fracErrs[0])
            self.c.specifyUncertainty(uncName,sbWnameS,"tt",1.-fracErrs[1])
        # !*! need to change from error on yield to error on fraction since total normalization fluctuation
        #     already accounted for!!!
        # !*! should be correlated between MB SR and CR
        #
        ys = [ mbres["yW_srNJet_0b_lowDPhi"], mbres["yTT_srNJet_0b_lowDPhi"], \
                   mbres["yRest_srNJet_0b_lowDPhi_truth"],  mbres["yQCD_srNJet_0b_lowDPhi"] ]
        # temporary fix for QCD variance 
        vQCD = mbres["yQCD_Var_srNJet_0b_lowDPhi"]
        if isnan(vQCD):
            print "Replacing nan for yQCD_Var_srNJet_0b_lowDPhi in ",mbnameC
            vQCD =  mbres["yQCD_srNJet_0b_lowDPhi"]**2
        vys = [ mbres["yW_Var_srNJet_0b_lowDPhi"], mbres["yTT_Var_srNJet_0b_lowDPhi"], \
                   mbres["yRest_Var_srNJet_0b_lowDPhi_truth"],  vQCD ]
        fracErrs = relErrorsOnFractions(ys,vys)
        uncName = "yWtt" + mbnameC
        self.c.addUncertainty(uncName,"lnN")
        self.c.specifyUncertainty(uncName,mbnameC,"W",1.+fracErrs[0])
        self.c.specifyUncertainty(uncName,mbnameC,"tt",1.-fracErrs[1])
        self.c.specifyUncertainty(uncName,mbnameS,"W",1.+fracErrs[0])
        self.c.specifyUncertainty(uncName,mbnameS,"tt",1.-fracErrs[1])

      #
      # other systematics on (total) prediction in MB/SR
      #
      for mbname in mbBinNames:
        bname = mbname[2:]
        mbnameC = mbname + "C"
        mbnameS = mbname + "S"
        mbres = self.subDict(self.bkgres,self.mbBins[mbname])
        self.mbsigres = self.subDict(self.sigres,self.mbBins[mbname])

        sbWname = "J3" + bname
        # sbWnameC = sbWnameBase + "C"
        # sbWresC = self.subDict(self.bkgres,self.sbBins[sbWnameC])
        sbWnameS = sbWname + "S"
        sbWresS = self.subDict(self.bkgres,self.sbBins[sbWname])

        sbttname = "J4" + bname
        # sbttnameC = sbttnameBase + "C"
        # sbttresC = self.subDict(self.bkgres,self.sbBins[sbttname])
        sbttnameS = sbttname + "S"
        sbttresS = self.subDict(self.bkgres,self.sbBins[sbttname])
        # uncertainty on RCS_W (e+mu)/mu
        uncName = "rcsWemu" + mbnameS
        self.c.addUncertainty(uncName,"lnN")
        self.c.specifyUncertainty(uncName,mbnameS,"W",1.+mbres["systematics"]["ratio_mu_elemu"])
        # uncertainty on b tagging
        if not "btag" in self.c.uncertainties:
            self.c.addUncertainty("btag","lnN")
        self.c.specifyUncertainty("btag",mbnameS,"W",1.+mbres["systematics"]["btagSF"])
        self.c.specifyUncertainty("btag",mbnameS,"tt",1.+mbres["systematics"]["btagSF"])
        self.c.specifyUncertainty("btag",mbnameS,"other",1.+mbres["systematics"]["btagSF"])
        # uncertainty on top pt (!*! should rescale total rel. uncertainty for application on ttbar only)
        if not "topPt" in self.c.uncertainties:
            self.c.addUncertainty("topPt","lnN")
        self.c.specifyUncertainty("topPt",mbnameS,"W",1.+mbres["systematics"]["topPt"])
        self.c.specifyUncertainty("topPt",mbnameS,"tt",1.+mbres["systematics"]["topPt"])
        self.c.specifyUncertainty("topPt",mbnameS,"other",1.+mbres["systematics"]["topPt"])
        # uncertainty on lepton SFs
        if not "leptonSF" in self.c.uncertainties:
            self.c.addUncertainty("leptonSF","lnN")
        self.c.specifyUncertainty("leptonSF",mbnameS,"signal",1.+mbres["systematics"]["lepSF"])
        self.c.specifyUncertainty("leptonSF",mbnameS,"W",1.+mbres["systematics"]["lepSF"])
        self.c.specifyUncertainty("leptonSF",mbnameS,"tt",1.+mbres["systematics"]["lepSF"])
        self.c.specifyUncertainty("leptonSF",mbnameS,"other",1.+mbres["systematics"]["lepSF"])
        # stat uncertainty on kappaW, kappaTT and kappa_b
        self.c.addUncertainty("kappaW"+mbnameS,"lnN")
        self.c.specifyUncertainty("kappaW"+mbnameS,mbnameS,"W",1.+mbres["systematics"]["kappa_W"])
        self.c.addUncertainty("kappaTT"+mbnameS,"lnN")
        self.c.specifyUncertainty("kappaTT"+mbnameS,mbnameS,"tt",1.+mbres["systematics"]["kappa_TT"])
        self.c.addUncertainty("kappab"+mbnameS,"lnN")
        self.c.specifyUncertainty("kappab"+mbnameS,mbnameS,"tt",1.+mbres["systematics"]["kappa_b"])
        # Rcs systematics W and tt ("linear fit")
        self.c.addUncertainty("rcsW"+mbnameS,"lnN")
        self.c.specifyUncertainty("rcsW"+mbnameS,mbnameS,"W",1.+mbres["systematics"]["rcs_W"])
        self.c.addUncertainty("rcsTT"+mbnameS,"lnN")
        self.c.specifyUncertainty("rcsTT"+mbnameS,mbnameS,"tt",1.+mbres["systematics"]["rcs_tt"])
        # QCD systematics
        self.c.addUncertainty("QCD"+mbnameS,"lnN")
        self.c.specifyUncertainty("QCD"+mbnameS,mbnameS,"W",1.+mbres["systematics"]["QCD"])
        self.c.specifyUncertainty("QCD"+mbnameS,mbnameS,"tt",1.+mbres["systematics"]["QCD"])
        self.c.specifyUncertainty("QCD"+mbnameS,mbnameS,"other",1.+mbres["systematics"]["QCD"])
        # dilepton
        self.c.addUncertainty("DiLep"+mbnameS,"lnN")
        self.c.specifyUncertainty("DiLep"+mbnameS,mbnameS,"W",1.+mbres["systematics"]["dilep"])
        self.c.specifyUncertainty("DiLep"+mbnameS,mbnameS,"tt",1.+mbres["systematics"]["dilep"])
        self.c.specifyUncertainty("DiLep"+mbnameS,mbnameS,"other",1.+mbres["systematics"]["dilep"])
        # PU systematics
        if not "PU" in self.c.uncertainties:
            self.c.addUncertainty("PU","lnN")
        self.c.specifyUncertainty("PU",mbnameS,"W",1.+mbres["systematics"]["pileup"])
        self.c.specifyUncertainty("PU",mbnameS,"tt",1.+mbres["systematics"]["pileup"])
        self.c.specifyUncertainty("PU",mbnameS,"other",1.+mbres["systematics"]["pileup"])
        # Cross sections & W polarization
        if not "xsecW" in self.c.uncertainties:
            self.c.addUncertainty("xsecW","lnN")
        self.c.specifyUncertainty("xsecW",mbnameS,"W",1.+mbres["systematics"]["Wxsec"])
        self.c.specifyUncertainty("xsecW",mbnameS,"tt",1.+mbres["systematics"]["Wxsec"])
        self.c.specifyUncertainty("xsecW",mbnameS,"other",1.+mbres["systematics"]["Wxsec"])
        if not "xsecTT" in self.c.uncertainties:
            self.c.addUncertainty("xsecTT","lnN")
        self.c.specifyUncertainty("xsecTT",mbnameS,"W",1.+mbres["systematics"]["TTxsec"])
        self.c.specifyUncertainty("xsecTT",mbnameS,"tt",1.+mbres["systematics"]["TTxsec"])
        self.c.specifyUncertainty("xsecTT",mbnameS,"other",1.+mbres["systematics"]["TTxsec"])
        if not "WPol" in self.c.uncertainties:
            self.c.addUncertainty("WPol","lnN")
        self.c.specifyUncertainty("WPol",mbnameS,"W",1.+mbres["systematics"]["Wpol"])
        self.c.specifyUncertainty("WPol",mbnameS,"tt",1.+mbres["systematics"]["Wpol"])
        self.c.specifyUncertainty("WPol",mbnameS,"other",1.+mbres["systematics"]["Wpol"])
        # stat. uncertainty on signal efficiency
        uncName = "statSeff" + mbnameS
        self.c.addUncertainty(uncName,"lnN")
        self.c.specifyUncertainty(uncName,mbnameS,"signal",1+self.sigSubDict(self.mbsigres)["err_MB_SR"])
      #      # WORST CASE SYST
      #      uncName = "worst"+mbnameS
      #      self.c.addUncertainty(uncName,"lnN")
      #      self.c.specifyUncertainty(uncName,mbnameS,"W",1+worstCaseSyst[self.mbBins[mbname][0]][self.mbBins[mbname][1]][self.mbBins[mbname][2]])

      #
      # SB related systematics propagated to MB
      #
      for sbname in sbBinNames:
        sbres = self.subDict(self.bkgres,self.sbBins[sbname])
        sbnameC = sbname + "C"
        sbnameS = sbname + "S"
        # statistical uncertainty SB CRs
        uncName = "stat" + sbnameC
        self.c.addUncertainty(uncName,"lnN")
        for mbname in mbBinNames:
          bname = mbname[2:]
          mbnameS = mbname + "S"
          if "J3"+bname==sbname:
              wYield = sbres["yW_crNJet_0b_lowDPhi"]
              wVar = sbres["yW_Var_crNJet_0b_lowDPhi"]
              self.c.specifyUncertainty(uncName,mbnameS,"W",1.+sqrt(wVar)/wYield)
          elif "J4"+bname==sbname:
              ttYield = sbres["y_crNJet_1b_lowDPhi"]
              ttVar = sbres["y_crNJet_1b_lowDPhi"]
              self.c.specifyUncertainty(uncName,mbnameS,"tt",1.+sqrt(ttVar)/ttYield)


      #
      # QCD uncertainties in CR / SB (!*! double counting?)
      #
      useAddQCD = False
      if useAddQCD:
          for sbname in sbBinNames:
            sbnameS = sbname + "S"
            sbres = self.subDict(self.bkgres,self.sbBins[sbname])
            self.c.addUncertainty("qcd"+sbnameS,"lnN")
            if sbname.startswith("J3"):
                self.c.specifyUncertainty("qcd"+sbnameS,sbnameS,"QCD",2.0)
            elif sbname.startswith("J4"):
                vQCD = sbres["yQCD_Var_crNJet_1b_highDPhi"]
                if isnan(vQCD):
                    print "Replacing nan for sbres yQCD_Var_crNJet_1b_highDPhi in ",sbnameS
                    vQCD =  sbres["yQCD_crNJet_1b_highDPhi"]**2
                print "QCD",sbname,sbres["yQCD_crNJet_1b_highDPhi"],sbres["yQCD_Var_crNJet_1b_highDPhi"]
                self.c.specifyUncertainty("qcd"+sbnameS,sbnameS,"QCD", \
                                              relErrForLimit(sbres["yQCD_crNJet_1b_highDPhi"],vQCD))
          for mbname in mbBinNames:
            mbnameS = mbname + "C"
            mbres = self.subDict(self.bkgres,self.mbBins[mbname])
            self.c.addUncertainty("qcd"+mbnameC,"lnN")
            # temporary fix for QCD variance 
            vQCD = mbres["yQCD_Var_srNJet_0b_lowDPhi"]
            if isnan(vQCD):
                vQCD =  mbres["yQCD_srNJet_0b_lowDPhi"]**2
            self.c.specifyUncertainty("qcd"+mbnameC,mbnameC,"QCD", \
                                          relErrForLimit(sbres["yQCD_srNJet_0b_lowDPhi"],vQCD))

          
      txtname = self.name + ".txt"
      self.c.writeToFile(txtname)
      #
      # comments
      #
      txt = open(txtname,"a")
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
      txt.write("# sigSyst ........... approximated total signal systematics\n")
      txt.close()

      if self.runLimit:
          stdout = sys.stdout
          sys.stdout = open(self.name+".log","a")
          opts = ""
          if self.runBlind:
              opts = "--run blind"
          print 'Result ',mbBinNames[0]," , ",self.signal["name"],self.signal["mglu"],self.signal["mlsp"]," : ",self.c.calcLimit(options=opts)
          sys.stdout.close()
          sys.stdout = stdout

