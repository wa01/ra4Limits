from cardFileWriter import cardFileWriter
#from limit_helper import plotsignif , plotLimit , signal_bins_3fb
from math import exp,sqrt
import os,sys

def relErrForLimit(value,variance,sign=1):
    result = 1.+sign*sqrt(variance)/value
    if result<0.:
        result = 0.01
    return result

class CalcSingleLimit:

    def __init__(self,bkgres,sbBinNames,sbBins,mbBinNames,mbBins,sigres,signal):

        self.corrSystSize = 1000.00
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

      #
      # define bins
      #
      for sbname in self.sbBinNames:
        sbnameS = sbname + "S"
        self.c.addBin(sbnameS,self.procNames,sbnameS)
        sbres = self.subDict(self.bkgres,self.sbBins[sbname])
        self.sbsigres = self.subDict(self.sigres,self.sbBins[sbname])
        r = "S"
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
        if sbnameS[:2]=="J3":
          # observation
          y_truth = sbres["yW_crNJet_0b_"+rDPhi+"_truth"] + \
              sbres["yTT_crNJet_0b_"+rDPhi+"_truth"] + \
              sbres["yRest_crNJet_0b_"+rDPhi+"_truth"]
          self.c.specifyObservation(sbnameS,int(y_truth+0.5))
          self.c.specifyExpectation(sbnameS,"signal",self.sigSubDict(self.sbsigres)['yield_SB_W_SR'])
        #      self.c.specifyExpectation(sbnameS,"W",sbres["y_crNJet_0b_"+rDPhi]-sbres["yTT_crNJet_0b_"+rDPhi])
          self.c.specifyExpectation(sbnameS,"W",sbres["yW_crNJet_0b_"+rDPhi])
          self.c.specifyExpectation(sbnameS,"tt",sbres["yTT_crNJet_0b_"+rDPhi])
          self.c.specifyExpectation(sbnameS,"other",0.001) # others are neglected in yield
          self.c.specifyExpectation(sbnameS,"QCD",0.001)
        # tt
        elif sbnameS[:2]=="J4":
          # observation
          y_truth = sbres["yW_crNJet_1b_"+rDPhi+"_truth"] + \
              sbres["yTT_crNJet_1b_"+rDPhi+"_truth"] + \
              sbres["yRest_crNJet_1b_"+rDPhi+"_truth"]
          self.c.specifyObservation(sbnameS,int(y_truth+0.5))
          self.c.specifyExpectation(sbnameS,"signal",self.sigSubDict(self.sbsigres)['yield_SB_tt_SR'])
          self.c.specifyExpectation(sbnameS,"W",0.001)
          self.c.specifyExpectation(sbnameS,"tt",sbres["y_crNJet_1b_"+rDPhi])
          self.c.specifyExpectation(sbnameS,"other",0.001) # others are neglected in yield
          self.c.specifyExpectation(sbnameS,"QCD",0.001)


      for mbname in self.mbBinNames:
        mbres = self.subDict(self.bkgres,self.mbBins[mbname])
        self.mbsigres = self.subDict(self.sigres,self.mbBins[mbname])
        #
        # low dPhi
        #
        mbnameC = mbname + "C"
        self.c.addBin(mbnameC,self.procNames,mbnameC)
        rDPhi = "lowDPhi"
        # observation
        y_truth = mbres["yW_srNJet_0b_"+rDPhi+"_truth"] + \
            mbres["yTT_srNJet_0b_"+rDPhi+"_truth"] + \
            mbres["yRest_srNJet_0b_"+rDPhi+"_truth"] + \
            mbres["yQCD_srNJet_0b_"+rDPhi+"_truth"]
        self.c.specifyObservation(mbnameC,int(y_truth+0.5))
        # expectation
        self.c.specifyExpectation(mbnameC,"signal",self.sigSubDict(self.mbsigres)['yield_MB_CR']) # to be corrected!
        self.c.specifyExpectation(mbnameC,"tt",mbres["yTT_srNJet_0b_"+rDPhi])
        self.c.specifyExpectation(mbnameC,"W",mbres["yW_srNJet_0b_"+rDPhi])
        self.c.specifyExpectation(mbnameC,"other",mbres["yRest_srNJet_0b_"+rDPhi+"_truth"])
        # !!!!! to be changed
        self.c.specifyExpectation(mbnameC,"QCD",0.001)
        # self.c.specifyExpectation(mbnameC,"QCD",mbres["yQCD_srNJet_0b"])
        #
        # high dPhi
        #
        mbnameS = mbname + "S"
        self.c.addBin(mbnameS,self.procNames,mbnameS)
        rDPhi = "highDPhi"
        # observation
        y_truth = mbres["W_truth"] +  mbres["TT_truth"] + mbres["Rest_truth"]
        self.c.specifyObservation(mbnameS,int(y_truth+0.5))
        # expectation
        self.c.specifyExpectation(mbnameS,"signal",self.sigSubDict(self.mbsigres)['yield_MB_SR']) # to be corrected!
        self.c.specifyExpectation(mbnameS,"tt",mbres["TT_pred"])
        self.c.specifyExpectation(mbnameS,"W",mbres["W_pred"])
        self.c.specifyExpectation(mbnameS,"other",mbres["Rest_truth"])
        self.c.specifyExpectation(mbnameS,"QCD",0.001)

      #
      # global uncertainties
      #
      # self.c.addUncertainty("worst","lnN")
      self.c.addUncertainty("btag","lnN")
      self.c.addUncertainty("lumi","lnN")
      self.c.addUncertainty("sigSyst","lnN")
      for bname in self.sbBinNames:
        sbname = bname + "S"
        self.c.specifyUncertainty("lumi",sbname,"signal",1.046)
        self.c.specifyUncertainty("sigSyst",sbname,"signal",1.20) # to be corrected!
        self.c.specifyUncertainty("lumi",sbname,"other",1.046)
      for bname in self.mbBinNames:
        for r in [ "C", "S" ]:
          mbname = bname + r
          print mbname
          self.c.specifyUncertainty("lumi",mbname,"signal",1.046)
          self.c.specifyUncertainty("sigSyst",mbname,"signal",1.20) # to be corrected!
          self.c.specifyUncertainty("lumi",mbname,"other",1.046)

      for mbname in self.mbBinNames:
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
        #
        # signal regions
        #
        # correlation W regions: B and F / E and F
        uncName = "corrWBF" + mbname
        self.c.addUncertainty(uncName,"lnN")
        self.c.specifyUncertainty(uncName,"J3"+bname+"S","W",self.corrSystSize)
        self.c.specifyUncertainty(uncName,mbnameS,"W",self.corrSystSize)
        uncName = "corrWEF" + mbname
        self.c.addUncertainty(uncName,"lnN")
        self.c.specifyUncertainty(uncName,mbnameC,"W",self.corrSystSize)
        self.c.specifyUncertainty(uncName,mbnameS,"W",self.corrSystSize)
        # correlation tt regions: D and F / E and F
        uncName = "corrTTDF" + mbname
        self.c.addUncertainty(uncName,"lnN")
        self.c.specifyUncertainty(uncName,"J4"+bname+"S","tt",self.corrSystSize)
        self.c.specifyUncertainty(uncName,mbnameS,"tt",self.corrSystSize)
        uncName = "corrTTEF" + mbname
        self.c.addUncertainty(uncName,"lnN")
        self.c.specifyUncertainty(uncName,mbnameC,"tt",self.corrSystSize)
        self.c.specifyUncertainty(uncName,mbnameS,"tt",self.corrSystSize)
        # anticorrelated W/tt yields from fit
        uncName = "yWtt" + mbnameC
        self.c.addUncertainty(uncName,"lnN")
        self.c.specifyUncertainty(uncName,mbnameC,"W",relErrForLimit(mbres["yW_srNJet_0b_lowDPhi"],mbres["yW_Var_srNJet_0b_lowDPhi"]))
        self.c.specifyUncertainty(uncName,mbnameC,"tt",relErrForLimit(mbres["yTT_srNJet_0b_lowDPhi"],mbres["yTT_Var_srNJet_0b_lowDPhi"],-1))
        # uncertainty on kappas
        uncName = "kappa" + mbnameS
        self.c.addUncertainty(uncName,"lnN")
        self.c.specifyUncertainty(uncName,mbnameS,"W",1.25) # to be corrected!
        self.c.specifyUncertainty(uncName,mbnameS,"tt",1.25) # to be corrected!
        # uncertainty on RCS_W (e+mu)/mu
        uncName = "rcsWemu" + mbnameS
        self.c.addUncertainty(uncName,"lnN")
        self.c.specifyUncertainty(uncName,mbnameS,"W",1.05) # to be corrected!
        # uncertainty on b tagging
        self.c.specifyUncertainty("btag",mbnameS,"W",1.05) # to be corrected!
        self.c.specifyUncertainty("btag",mbnameS,"W",1.05) # to be corrected!
        self.c.specifyUncertainty("btag",mbnameS,"other",1.05) # to be corrected!
        # stat. uncertainty on signal efficiency
        uncName = "statSeff" + mbnameS
        self.c.addUncertainty(uncName,"lnN")
        self.c.specifyUncertainty(uncName,mbnameS,"signal",1+self.sigSubDict(self.mbsigres)["err_MB_SR"])
      #      # WORST CASE SYST
      #      uncName = "worst"+mbnameS
      #      self.c.addUncertainty(uncName,"lnN")
      #      self.c.specifyUncertainty(uncName,mbnameS,"W",1+worstCaseSyst[self.mbBins[mbname][0]][self.mbBins[mbname][1]][self.mbBins[mbname][2]])


      for sbname in self.sbBinNames:
        sbres = self.subDict(self.bkgres,self.sbBins[sbname])
        sbnameC = sbname + "C"
        sbnameS = sbname + "S"
        # statistical uncertainty SB CRs
        obs = self.c.observation[sbnameS] # to be corrected !!!
        assert obs>0
        unc = 1. + 1./sqrt(obs)
        uncName = "stat" + sbnameC
        self.c.addUncertainty(uncName,"lnN")
        for mbname in self.mbBinNames:
          bname = mbname[2:]
          mbnameS = mbname + "S"
          if "J3"+bname+"C"==sbname:
            self.c.specifyUncertainty(uncName,mbname,"W",unc)
          elif "J4"+bname+"C"==sbname:
            self.c.specifyUncertainty(uncName,mbname,"tt",unc)
        # anticorrelated W/tt yields from fit
        if sbname.startswith("J3"):
          uncName = "yWtt" + sbnameS
          self.c.addUncertainty(uncName,"lnN")
          self.c.specifyUncertainty(uncName,sbnameS,"W",relErrForLimit(sbres["yW_crNJet_0b_highDPhi"], 
                                                                 sbres["yW_Var_crNJet_0b_highDPhi"]))
          self.c.specifyUncertainty(uncName,sbnameS,"tt",relErrForLimit(sbres["yTT_crNJet_0b_highDPhi"], \
                                                                    sbres["yTT_Var_crNJet_0b_highDPhi"],-1))
          if uncName=="yWttJ3L3H1D2S":
            print sbres["yW_crNJet_0b_highDPhi"],sbres["yW_Var_crNJet_0b_highDPhi"], \
                relErrForLimit(sbres["yW_crNJet_0b_highDPhi"],sbres["yW_Var_crNJet_0b_highDPhi"])
            print sbres["yTT_crNJet_0b_highDPhi"],sbres["yTT_Var_crNJet_0b_highDPhi"], \
                relErrForLimit(sbres["yTT_crNJet_0b_highDPhi"],sbres["yTT_Var_crNJet_0b_highDPhi"],-1)

      #
      # QCD
      #
      #  for bname in self.mbBinNames:
      #    if bname.endswith("C"):
      #      self.c.addUncertainty(uncName,"lnN")
      #      self.c.specifyUncertainty(uncName,bname,"QCD",relErrForLimit()
      #    self.c.specifyUncertainty("sigSyst",bname,"signal",1.20)
      #    self.c.specifyUncertainty("lumi",bname,"other",1.046)

      self.c.writeToFile("calc_limit.txt")
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
      return

      stdout = sys.stdout
      sys.stdout = open("results.log","a")
      print 'Result ',self.mbBinNames[0]," , ",signal["name"]," : ",self.c.calcLimit(options="--run blind")
      sys.stdout.close()
      sys.stdout = stdout

