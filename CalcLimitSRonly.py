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
      print "Using mb bins ",mbBinNames
      #
      # bin definition; observed and expected counts
      #
      for mbname in mbBinNames:
        mbres = self.subDict(self.bkgres,self.mbBins[mbname])
        self.mbsigres = self.subDict(self.sigres,self.mbBins[mbname])
        #
        # high dPhi
        #
        mbnameS = mbname + "S"
        self.c.addBin(mbnameS,self.procNames,mbnameS)
        rDPhi = "highDPhi"
        # observation
        y_truth = mbres["W_truth"] +  mbres["TT_truth"] + mbres["Rest_truth"]
        self.c.specifyObservation(mbnameS,int(y_truth+0.5)) # !*! to be corrected
        # expectation
        self.c.specifyExpectation(mbnameS,"signal",self.sigSubDict(self.mbsigres)['yield_MB_SR'])
        self.c.specifyExpectation(mbnameS,"tt",mbres["TT_pred"])
        self.c.specifyExpectation(mbnameS,"W",mbres["W_pred"])
        self.c.specifyExpectation(mbnameS,"other",mbres["Rest_truth"])
        self.c.specifyExpectation(mbnameS,"QCD",0.001)

      #
      # global uncertainties
      #
      # self.c.addUncertainty("worst","lnN")
#      self.c.addUncertainty("btag","lnN")
      self.c.addUncertainty("lumi","lnN")
      self.c.addUncertainty("sigSyst","lnN")
      self.c.addUncertainty("xsecOther","lnN")
      for bname in mbBinNames:
        for r in [ "S" ]:
          mbname = bname + r
          self.c.specifyUncertainty("lumi",mbname,"signal",1.046)
          self.c.specifyUncertainty("sigSyst",mbname,"signal",1.20) # to be corrected!
          self.c.specifyUncertainty("lumi",mbname,"other",1.046)
          self.c.specifyUncertainty("xsecOther",mbname,"other",1.50)
      #
      # other systematics on (total) prediction in MB/SR
      #
      for mbname in mbBinNames:
        bname = mbname[2:]
        mbnameS = mbname + "S"
        mbres = self.subDict(self.bkgres,self.mbBins[mbname])
        self.mbsigres = self.subDict(self.sigres,self.mbBins[mbname])

        # total uncertainty
#        uncName = "tot" + mbnameS
#        self.c.addUncertainty(uncName,"lnN")
#        self.c.specifyUncertainty(uncName,mbnameS,"W",relErrForLimit(mbres["W_pred"],mbres["W_pred_err"]**2))
#        self.c.specifyUncertainty(uncName,mbnameS,"tt",relErrForLimit(mbres["TT_pred"],mbres["TT_pred_err"]**2))
#        self.c.specifyUncertainty(uncName,mbnameS,"other",relErrForLimit(mbres["Rest_truth"],mbres["Rest_truth_err"]**2))
        uncName = "W" + mbnameS
        self.c.addUncertainty(uncName,"lnN")
        self.c.specifyUncertainty(uncName,mbnameS,"W",relErrForLimit(mbres["W_pred"],mbres["W_pred_err"]**2))
        uncName = "tt" + mbnameS
        self.c.addUncertainty(uncName,"lnN")
        self.c.specifyUncertainty(uncName,mbnameS,"tt",relErrForLimit(mbres["TT_pred"],mbres["TT_pred_err"]**2))
        uncName = "other" + mbnameS
        self.c.addUncertainty(uncName,"lnN")
        self.c.specifyUncertainty(uncName,mbnameS,"other",relErrForLimit(mbres["Rest_truth"],mbres["Rest_truth_err"]**2))
        # stat. uncertainty on signal efficiency
        uncName = "statSeff" + mbnameS
        self.c.addUncertainty(uncName,"lnN")
        self.c.specifyUncertainty(uncName,mbnameS,"signal",1+self.sigSubDict(self.mbsigres)["err_MB_SR"])
      #      # WORST CASE SYST
      #      uncName = "worst"+mbnameS
      #      self.c.addUncertainty(uncName,"lnN")
      #      self.c.specifyUncertainty(uncName,mbnameS,"W",1+worstCaseSyst[self.mbBins[mbname][0]][self.mbBins[mbname][1]][self.mbBins[mbname][2]])



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
          print 'Result ',mbBinNames[0]," , ",self.signal["name"],self.signal["mglu"],self.signal["mlsp"]," : ",self.c.calcLimit(options="--run blind")
          sys.stdout.close()
          sys.stdout = stdout

