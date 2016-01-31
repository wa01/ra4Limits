import ROOT
from array import array

filename = ""

def toGraph2D(name,title,length,x,y,z):
    result = ROOT.TGraph2D(name,title)
    for i in range(length):
        result.SetPoint(i,x[i],y[i],z[i])
    return result

def toGraph(name,title,length,x,y):
    result = ROOT.TGraph(name,title)
    for i in range(length):
        result.SetPoint(i,x[i],y[i])
    return result

def SetupColors():
    num = 5
    bands = 255
    colors = [ ]
    stops = [0.00, 0.34, 0.61, 0.84, 1.00]
    red = [0.50, 0.50, 1.00, 1.00, 1.00]
    green = [0.50, 1.00, 1.00, 0.60, 0.50]
    blue = [1.00, 1.00, 0.50, 0.40, 0.50]
    arr_stops = array('d', stops)
    arr_red = array('d', red)
    arr_green = array('d', green)
    arr_blues = array('d', blues)
    # num = 6
    # red[num] =   {1.,0.,0.,0.,1.,1.}
    # green[num] = {0.,0.,1.,1.,1.,0.}
    # blue[num] =  {1.,1.,1.,0.,0.,0.}
    # stops[num] = {0.,0.2,0.4,0.6,0.8,1.}*/
    fi = ROOT.TColor::CreateGradientColorTable(num,arr_stops,arr_red,arr_green,arr_blue,bands)
    for i in range(bands):
        colors.append(fi+i)
    arr_colors = array('i', colors)
    ROOT.gStyle.SetNumberContours(bands)
    ROOT.gStyle.SetPalette(bands, arr_colors)


def DrawContours (g2, color, style, leg, name):
    
    out = ROOT.TGraph()
    l = g2.GetContourList(1.)
    if not l:
        return out
    added = False
    max_points = -1
    for i in range(l.GetSize()):
        g = l.At(i)
        if not g:
            continue
        n_points = g.GetN()
        if n_points > max_points:
            out = g
            max_points = n_points
        g.SetLineColor(color)
        g.SetLineStyle(style)
        g.SetLineWidth(5)
        g.Draw("L same")
        if ( not added ) and leg and name != "":
            leg.AddEntry(g, name.c_str(), "l")
            added = true
    return out

SetupColors()

assert filename!=""
vmx = [ ]
vmy = [ ]
vxsec = [ ]
vobs = [ ]
vobsup = [ ]
vobsdown = [ ]
vexp = [ ]
vup = [ ]
vdown = [ ]


start_read = False
for line in open(infile):
    if line.find("---")>=0:
        start_read = True
    elif start_read:
        fields = line[:-1].split()
        pmx,pmy,pxsec,pobs,pobsup,pobsdown,pexp,pup,pdown = [ float(x) for x in line[:-1].split() ]
        vmx.append(pmx)
        vmy.append(pmy)
        vxsec.append(pxsec)
        vobs.append(pobs)
        vobsup.append(pobsup)
        vobsdown.append(pobsdown)
        vexp.append(pexp)
        vup.append(pup)
        vdown.append(pdown)

assert len(vmx) > 2)
assert not ( len(vmx) != len(vmy) \
                 or len(vmx) != len(vxsec) \
                 or len(vmx) != len(vobs) \
                 or len(vmx) != len(vobsup) \
                 or len(vmx) != len(vobsdown) \
                 or len(vmx) != len(vexp) \
                 or len(vmx) != len(vup) \
                 or len(vmx) != len(vdown) ) 
  
  
vlim = [ ]
for i in range(len(vxsec)):
    vlim.append(vxsec[i]*vobs[i])

glim = toGraph2D("glim", "Cross-Section Limit", len(vlim), vmx, vmy, vlim)
gobs = toGraph2D("gobs", "Observed Limit", len(vobs), vmx, vmy, vobs)
gobsup = toGraph2D("gobsup", "Observed +1#sigma Limit", len(vobsup), vmx, vmy, vobsup)
gobsdown = toGraph2D("gobsdown", "Observed -1#sigma Limit", len(vobsdown), vmx, vmy, vobsdown)
gexp = toGraph2D("gexp", "Expected Limit", len(vexp), vmx, vmy, vexp)
gup = toGraph2D("gup", "Expected +1#sigma Limit", len(vup), vmx, vmy, vup)
gdown = toGraph2D("gdown", "Expected -1#sigma Limit", len(vdown), vmx, vmy, vdown)
dots = toGraph(len(vmx), vmx, vmy)

xmin = min(vmx)
xmax = max(vmx)
ymin = min(vmy)
ymax = max(vmy)
bin_size = 12.5
int nxbins = max(1, min(500, int((xmax-xmin+bin_size/100.)/bin_size)))
int nybins = max(1, min(500, int((ymax-ymin+bin_size/100.)/bin_size)))
glim.SetNpx(nxbins)
glim.SetNpy(nybins)

hlim = glim.GetHistogram()
assert hlim
hlim.SetTitle(";m_{gluino} [GeV];m_{LSP} [GeV]")
  
c = ROOT.TCanvas("","",800,800)
c.SetLogz()
hlim.SetMinimum(min(vlim))
hlim.SetMaximum(max(vlim))
glim.Draw("colz")
l = ROOT.TLegend(ROOT.gStyle.GetPadLeftMargin(), 1.-ROOT.gStyle.GetPadTopMargin(), \
                     1.-ROOT.gStyle.GetPadRightMargin(), 1.)
l.SetNColumns(2)
l.SetBorderSize(0)
cup = DrawContours(gup, 2, 2)
cdown = DrawContours(gdown, 2, 2)
cexp = DrawContours(gexp, 2, 1, l, "Expected")
cobsup = DrawContours(gobsup, 1, 2)
cobsdown = DrawContours(gobsdown, 1, 2)
cobs = DrawContours(gobs, 1, 1, l, "Observed")
l.Draw("same")
dots.Draw("p same")
c.Print("limit_scan.pdf")

tfile = ROOT.TFile("limit_scan.root","recreate")
hlim.Write("hXsec_exp_corr")
cobs.Write("graph_smoothed_Obs")
cobsup.Write("graph_smoothed_ObsP")
cobsdown.Write("graph_smoothed_ObsM")
cexp.Write("graph_smoothed_Exp")
cup.Write("graph_smoothed_ExpP")
cdown.Write("graph_smoothed_ExpM")
tfile.Close()

  
