#!/usr/bin/env python
import re
from sys import argv, stdout, stderr, exit
from optparse import OptionParser

# import ROOT with a fix to get batch mode (http://root.cern.ch/phpBB3/viewtopic.php?t=3198)
hasHelp = False
#for X in ("-h", "-?", "--help"):
#    if X in argv:
#        hasHelp = True
#        argv.remove(X)
#argv.append( '-b-' )
import ROOT
#ROOT.gROOT.SetBatch(True)
ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
#argv.remove( '-b-' )
if hasHelp: argv.append("-h")

parser = OptionParser(usage="usage: %prog [options] in.root  \nrun with --help to get list of options")
parser.add_option("--vtol", "--val-tolerance", dest="vtol", default=0.30, type="float", help="Report nuisances whose value changes by more than this amount of sigmas")
parser.add_option("--stol", "--sig-tolerance", dest="stol", default=0.10, type="float", help="Report nuisances whose sigma changes by more than this amount")
parser.add_option("--vtol2", "--val-tolerance2", dest="vtol2", default=2.0, type="float", help="Report severely nuisances whose value changes by more than this amount of sigmas")
parser.add_option("--stol2", "--sig-tolerance2", dest="stol2", default=0.50, type="float", help="Report severely nuisances whose sigma changes by more than this amount")
parser.add_option("-a", "--all",      dest="all",    default=False,  action="store_true", help="Print all nuisances, even the ones which are unchanged w.r.t. pre-fit values.")
parser.add_option("-A", "--absolute", dest="abs",    default=False,  action="store_true", help="Report also absolute values of nuisance values and errors, not only the ones normalized to the input sigma")
parser.add_option("-p", "--poi",      dest="poi",    default="r",    type="string",  help="Name of signal strength parameter (default is 'r' as per text2workspace.py)")
parser.add_option("-f", "--format",   dest="format", default="text", type="string",  help="Output format ('text', 'latex', 'twiki'")
parser.add_option("-g", "--histogram", dest="plotfile", default=None, type="string", help="If true, plot the pulls of the nuisances to the given file.")

(options, args) = parser.parse_args()
if len(args) == 0:
    parser.print_usage()
    exit(1)

file = ROOT.TFile(args[0])
if file == None: raise RuntimeError, "Cannot open file %s" % args[0]
fit_b  = file.Get("fit_b")
fpf_b = fit_b.floatParsFinal()
prefit = file.Get("nuisances_prefit")
if fit_b == None or fit_b.ClassName()   != "RooFitResult": raise RuntimeError, "File %s does not contain the output of the background fit 'fit_b'" % args[0]
if prefit == None or prefit.ClassName() != "RooArgSet":    raise RuntimeError, "File %s does not contain the prefit nuisances 'nuisances_prefit'"  % args[0]

results = {
#    "DiLep" : { },
    "QCD" : { },
    "corrTTDF" : { },
    "corrTTEF" : { },
    "corrWBF" : { },
    "corrWEF" : { },
#    "jec" : { },
    "kappaTT" : { },
    "kappaW" : { },
    "kappab" : { },
#    "rcs" : { },
    "rcsWemu" : { },
    "statSeff" : { },
    "yWtt" : { },
#    "isr" : { },
    "stat" : { },
    "other" : { }
}


for i in range(fpf_b.getSize()):
    nuis_b = fpf_b.at(i)
    name   = nuis_b.GetName();
#    nuis_b = fpf_b.find(name)
    nuis_p = prefit.find(name)

    cat = "other"
    for k in sorted(results.keys(),reverse=True):
        if k!="other" and name.startswith(k):
            cat = k
            break
    assert not name in results[cat]
    results[cat][name] = { }

    mean_p, sigma_p = 0,0
    if nuis_p == None:
        if not options.abs: continue
        print "No prefit for ",name
        continue
    else:
        print name,nuis_p.getMin(),nuis_p.getMax()
        mean_p, sigma_p = (nuis_p.getVal(), nuis_p.getError())
        results[cat][name]["mean_p"] = mean_p
        results[cat][name]["sigma_p"] = sigma_p
    for fit_name, nuis_x in [('b', nuis_b)]:
        if nuis_x == None:
            pass
        else:
            results[cat][name]["nuis_val"] = nuis_x.getVal()
            results[cat][name]["nuis_err"] = nuis_x.getError()

ROOT.gStyle.SetOptTitle(0)
ROOT.gROOT.cd()
canvases = [ ]
gobjects = { }
for c in results:
    names = sorted(results[c].keys())
    nnames = len(names)
    if nnames==0:
        continue
    gobjects[c] = [ ]
    gobjects[c].append(ROOT.TH1F("h"+c,"h"+c,nnames,-0.5,nnames-0.5))
    gobjects[c][-1].SetStats(0)
    xaxis = gobjects[c][-1].GetXaxis()
    for i,n in enumerate(names):
        xaxis.SetBinLabel(i+1,n)
    for postfix in [ "pre", "post" ]:
        g = ROOT.TGraphErrors()
        g.SetName("g"+cat+postfix)
        gobjects[c].append(g)
        g.SetLineWidth(2)
    gobjects[c][1].SetLineColor(ROOT.kBlue)
    gobjects[c][2].SetLineColor(ROOT.kRed)
    ymin = None
    ymax = None
    for i,n in enumerate(names):
        y = results[c][n]["mean_p"]
        ey = results[c][n]["sigma_p"]
        if ymin==None or (y-ey)<ymin:
            ymin = y - ey
        if ymax==None or (y+ey)>ymax:
            ymax = y + ey
        gobjects[c][1].SetPoint(i,float(i)-0.20,y)
        gobjects[c][1].SetPointError(i,0.20,ey)

        y = results[c][n]["nuis_val"]
        ey = results[c][n]["nuis_err"]
        if ymin==None or (y-ey)<ymin:
            ymin = y - ey
        if ymax==None or (y+ey)>ymax:
            ymax = y + ey
        gobjects[c][2].SetPoint(i,float(i)+0.20,y)
        gobjects[c][2].SetPointError(i,0.20,ey)

    dy = ymax - ymin
    gobjects[c][0].SetMinimum(ymin-0.05*dy)
    gobjects[c][0].SetMaximum(ymax+0.05*dy)
    canvases.append(ROOT.TCanvas("c"+c,"c"+c,800,600))
    gobjects[c][0].Draw()
    gobjects[c][1].Draw("e0 z same")
    gobjects[c][2].Draw("e0 z same")
    leg = ROOT.TLegend(0.68,0.78,0.88,0.88)
    leg.SetFillColor(ROOT.kWhite)
    leg.SetBorderSize(0)
    leg.AddEntry(gobjects[c][1],"Prefit","l")
    leg.AddEntry(gobjects[c][2],"Post bkg fit","l")
    leg.Draw()
    canvases[-1].Update()
    canvases[-1].SaveAs("nuisances-"+c+".pdf")

raw_input("Enter")

