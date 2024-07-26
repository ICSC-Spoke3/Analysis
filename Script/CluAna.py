import numpy as np
from ROOT import *

import sys

TGaxis.SetMaxDigits(3)
gStyle.SetOptFit(1)

print sys.argv

fdata = "../RootData/"+sys.argv[1]

print "Opening data file ", fdata

fr = TFile(fdata)

ct = fr.Get("ClusterTree")

print ct.GetEntries()

h0H = TH1F("h0H", "Horizontal View; Number of Photo-Electrons; Number of Entries", 510, -0.5, 50.5)  
h1H = TH1F("h1H", "Horizontal View - Cluster Size==1; Number of Photo-Electrons; Number of Entries", 510, -0.5, 50.5)  
h2H = TH1F("h2H", "Horizontal View - Cluster Size==2; Number of Photo-Electrons; Number of Entries", 510, -0.5, 50.5)  
h3H = TH1F("h3H", "Horizontal View - Cluster Size==3; Number of Photo-Electrons; Number of Entries", 510, -0.5, 50.5)  

h0V = TH1F("h0V", "Vertical View; Number of Photo-Electrons; Number of Entries", 510, -0.5, 50.5)  
h1V = TH1F("h1V", "Vertical View - Cluster Size==1; Number of Photo-Electrons; Number of Entries", 510, -0.5, 50.5)  
h2V = TH1F("h2V", "Vertical View - Cluster Size==2; Number of Photo-Electrons; Number of Entries", 510, -0.5, 50.5)  
h3V = TH1F("h3V", "Vertical View - Cluster Size==3; Number of Photo-Electrons; Number of Entries", 510, -0.5, 50.5)  

h0H_1 = TH1F("h0H_1", "Horizontal View; Number of Photo-Electrons; Number of Events", 510, -0.5, 50.5)  
h1H_1 = TH1F("h1H_1", "Horizontal View - Cluster Size==1; Number of Photo-Electrons; Number of Events", 510, -0.5, 50.5)  
h2H_1 = TH1F("h2H_1", "Horizontal View - Cluster Size==2; Number of Photo-Electrons; Number of Events", 510, -0.5, 50.5)  
h3H_1 = TH1F("h3H_1", "Horizontal View - Cluster Size==3; Number of Photo-Electrons; Number of Events", 510, -0.5, 50.5)  

h0V_1 = TH1F("h0V_1", "Vertical View; Number of Photo-Electrons; Number of Events", 510, -0.5, 50.5)  
h1V_1 = TH1F("h1V_1", "Vertical View - Cluster Size==1; Number of Photo-Electrons; Number of Events", 510, -0.5, 50.5)  
h2V_1 = TH1F("h2V_1", "Vertical View - Cluster Size==2; Number of Photo-Electrons; Number of Events", 510, -0.5, 50.5)  
h3V_1 = TH1F("h3V_1", "Vertical View - Cluster Size==3; Number of Photo-Electrons; Number of Events", 510, -0.5, 50.5)  

ct.Draw("TnpeHor>>h0H","CsizHor>0&&TcluHor>0","goff")
ct.Draw("TnpeHor>>h1H","CsizHor==1&&TcluHor>0","goff")
ct.Draw("TnpeHor>>h2H","CsizHor==2&&TcluHor>0","goff")
ct.Draw("TnpeHor>>h3H","CsizHor==3&&TcluHor>0","goff")

ct.Draw("TnpeVer>>h0V","CsizVer>0&&TcluVer>0","goff")
ct.Draw("TnpeVer>>h1V","CsizVer==1&&TcluVer>0","goff")
ct.Draw("TnpeVer>>h2V","CsizVer==2&&TcluVer>0","goff")
ct.Draw("TnpeVer>>h3V","CsizVer==3&&TcluVer>0","goff")

c0H = TCanvas("c0H", "c0H")
c0H.SetLogy()

h0H.SetLineColor(8)
h0H.SetFillColor(8)
h0H.SetFillStyle(1001)
h0H.Draw("goff")
gPad.Update()
h0H.SetName("Cluster Size > 0 ")
st0H = h0H.GetListOfFunctions().FindObject("stats")
st0H.SetLineColor(8)
st0H.SetTextColor(8)
st0H.SetX1NDC(0.45)
st0H.SetX2NDC(0.60)
st0H.SetY1NDC(0.70)
st0H.SetY2NDC(0.85)


h1H.SetLineColor(2)
h1H.SetFillColor(2)
h1H.SetFillStyle(3001)
h1H.Draw("goff")
gPad.Update()
h1H.SetName("Cluster Size == 1")
h1H.Draw("goff")
gPad.Update()
st1H = h1H.GetListOfFunctions().FindObject("stats")
st1H.SetLineColor(2)
st1H.SetTextColor(2)
st1H.SetX1NDC(0.65)
st1H.SetX2NDC(0.80)
st1H.SetY1NDC(0.70)
st1H.SetY2NDC(0.85)

h2H.SetLineColor(4)
h2H.SetFillColor(4)
h2H.SetFillStyle(3002)
h2H.Draw("goff")
gPad.Update()
h2H.SetName("Cluster Size == 2")
h2H.Draw("goff")
gPad.Update()
st2H = h2H.GetListOfFunctions().FindObject("stats")
st2H.SetLineColor(4)
st2H.SetTextColor(4)
st2H.SetX1NDC(0.65)
st2H.SetX2NDC(0.80)
st2H.SetY1NDC(0.50)
st2H.SetY2NDC(0.65)

h3H.SetLineColor(6)
h3H.SetFillColor(6)
h3H.SetFillStyle(3002)
h3H.Draw("goff")
gPad.Update()
h3H.SetName("Cluster Size == 3")
h3H.Draw("goff")
gPad.Update()
st3H = h3H.GetListOfFunctions().FindObject("stats")
st3H.SetLineColor(6)
st3H.SetTextColor(6)
st3H.SetX1NDC(0.65)
st3H.SetX2NDC(0.80)
st3H.SetY1NDC(0.30)
st3H.SetY2NDC(0.45)

h0H.Draw()
h1H.Draw("same")
h2H.Draw("same")
h3H.Draw("same")


c0H.Update()

c0V = TCanvas("c0V", "c0V")
c0V.SetLogy()

h0V.SetLineColor(8)
h0V.SetFillColor(8)
h0V.SetFillStyle(1001)
h0V.Draw("goff")
gPad.Update()
h0V.SetName("Cluster Size > 0 ")
st0V = h0V.GetListOfFunctions().FindObject("stats")
st0V.SetLineColor(8)
st0V.SetTextColor(8)
st0V.SetX1NDC(0.45)
st0V.SetX2NDC(0.60)
st0V.SetY1NDC(0.70)
st0V.SetY2NDC(0.85)


h1V.SetLineColor(2)
h1V.SetFillColor(2)
h1V.SetFillStyle(3001)
h1V.Draw("goff")
gPad.Update()
h1V.SetName("Cluster Size == 1")
h1V.Draw("goff")
gPad.Update()
st1V = h1V.GetListOfFunctions().FindObject("stats")
st1V.SetLineColor(2)
st1V.SetTextColor(2)
st1V.SetX1NDC(0.65)
st1V.SetX2NDC(0.80)
st1V.SetY1NDC(0.70)
st1V.SetY2NDC(0.85)

h2V.SetLineColor(4)
h2V.SetFillColor(4)
h2V.SetFillStyle(3002)
h2V.Draw("goff")
gPad.Update()
h2V.SetName("Cluster Size == 2")
h2V.Draw("goff")
gPad.Update()
st2V = h2V.GetListOfFunctions().FindObject("stats")
st2V.SetLineColor(4)
st2V.SetTextColor(4)
st2V.SetX1NDC(0.65)
st2V.SetX2NDC(0.80)
st2V.SetY1NDC(0.50)
st2V.SetY2NDC(0.65)

h3V.SetLineColor(6)
h3V.SetFillColor(6)
h3V.SetFillStyle(3002)
h3V.Draw("goff")
gPad.Update()
h3V.SetName("Cluster Size == 3")
h3V.Draw("goff")
gPad.Update()
st3V = h3V.GetListOfFunctions().FindObject("stats")
st3V.SetLineColor(6)
st3V.SetTextColor(6)
st3V.SetX1NDC(0.65)
st3V.SetX2NDC(0.80)
st3V.SetY1NDC(0.30)
st3V.SetY2NDC(0.45)

h0V.Draw()
h1V.Draw("same")
h2V.Draw("same")
h3V.Draw("same")

c0V.Update()


npeth = 2.5
trasic = ["0", "1", "2", "3", "0 AND 1", "2 AND 3", "(0 OR 1) AND (2 OR 3)", "0 AND 1 AND 2 AND 3"]
h1trig = TH1F("h1trig", "Trigger Configuration (> 2.5 pe); ASIC Number; Fraction of Events", 8, 0., 8.)
h1trig.SetStats(0)

for j in range(8):
    h1trig.Fill(trasic[j], 0.)
    
htit = 'pe > 2.5; Time difference (ns); Number of Events'
h1time = TH1F("h1time", htit, 200, -20., 20.)
htit = 'pe > 2.5; Minimum Time difference (ns); Number of Events'
h1timeMin = TH1F("h1timeMin", htit, 200, -20., 20.)

atri = np.zeros(4)

ngood = 0
for i in range(ct.GetEntries()):

    ct.GetEntry(i)
    Tch = ct.Tch
    if(Tch<=0):
        continue
    ngood += 1
    Ipx = ct.Ipix
    Npe = ct.Npe
    Tns = ct.Tns
    Ipx = np.array(Ipx)
    Npe = np.array(Npe)
    Tns = np.array(Tns)
    Jpx = np.argsort(Ipx)
    Ipx = Ipx[Jpx]
    Npe = Npe[Jpx]
    Tns = Tns[Jpx]
    iV = np.searchsorted(Ipx,64)

    Ichi = ct.Ichip
    Ichi = np.array(Ichi)
    Ichi = Ichi[Jpx]
    #print Ichi

    TsizH = 0.
    TsizV = 0.
    TtimH = 0.
    TtimV = 0.
    TtimHmin = 99999.
    TtimVmin = 99999.
    
    atri *= 0
    for j in range(Tch):
        k = Ichi[j]
        ttim = Tns[j]
        if(Npe[j] > npeth):
            atri[k] = 1
            if(k==0 or k==1):
                TsizH += 1
                TtimH += ttim
                if(ttim<TtimHmin):
                    TtimHmin = ttim
            if(k==2 or k==3):
                TsizV += 1
                TtimV += ttim
                if(ttim<TtimVmin):
                    TtimVmin = ttim

    #print atri
    for j in range(4):
        if(atri[j]>0):
            h1trig.Fill(trasic[j], 1.)

    if(atri[0]*atri[1]>0):
        h1trig.Fill(trasic[4], 1.)

    if(atri[2]*atri[3]>0):
        h1trig.Fill(trasic[5], 1.)

    if( (atri[0]+atri[1])*(atri[2]+atri[3])>0):
        h1trig.Fill(trasic[6], 1.)

    if( atri[0]*atri[1]*atri[2]*atri[3] > 0):
        h1trig.Fill(trasic[7], 1.)

    if (TsizH*TsizV > 0):
        TtimH /= TsizH
        TtimV /= TsizV
        h1time.Fill(TtimH-TtimV)
        h1timeMin.Fill(TtimHmin-TtimVmin)
        
    '''
    print "Event ", i, Tch, np.argmax(Npe), Npe.max(), iV
    print Jpx
    print Ipx
    print Npe
    if(iV >0): 
        print "H --> ", np.argmax(Npe[:iV]), Npe[:iV].max(), Npe[:iV].sum()
    if(iV < Ipx.shape[0]): 
        print "V --> ", np.argmax(Npe[iV:]), Npe[iV:].max(), Npe[iV:].sum() 
    print "---"
    '''
    HTpe = -1
    VTpe = -1
    if(iV >0): 
        HTpe =  Npe[:iV].sum()
        h0H_1.Fill(HTpe)
        if(len(Npe[:iV])==1):
            h1H_1.Fill(HTpe)
          
        if(len(Npe[:iV])==2):
            h2H_1.Fill(HTpe)
        if(len(Npe[:iV])==3):
            h3H_1.Fill(HTpe)
          
    if(iV < Ipx.shape[0]): 
        VTpe =  Npe[iV:].sum()
        h0V_1.Fill(VTpe)
        if(len(Npe[iV:])==1):
            h1V_1.Fill(VTpe)
        if(len(Npe[iV:])==2):
            h2V_1.Fill(VTpe)
        if(len(Npe[iV:])==3):
            h3V_1.Fill(VTpe)
        
h1trig.Scale(1./ngood)

c1H = TCanvas("c1H", "c1H")
c1H.SetLogy()

h0H_1.SetLineColor(8)
h0H_1.SetFillColor(8)
h0H_1.SetFillStyle(1001)
h0H_1.Draw("goff")
gPad.Update()
h0H_1.SetName("Cluster Size > 0 ")
st0H_1 = h0H_1.GetListOfFunctions().FindObject("stats")
st0H_1.SetLineColor(8)
st0H_1.SetTextColor(8)
st0H_1.SetX1NDC(0.45)
st0H_1.SetX2NDC(0.60)
st0H_1.SetY1NDC(0.70)
st0H_1.SetY2NDC(0.85)

h1H_1.SetLineColor(2)
h1H_1.SetFillColor(2)
h1H_1.SetFillStyle(3001)
h1H_1.Draw("goff")
gPad.Update()
h1H_1.SetName("Cluster Size == 1")
h1H_1.Draw("goff")
gPad.Update()
st1H_1 = h1H_1.GetListOfFunctions().FindObject("stats")
st1H_1.SetLineColor(2)
st1H_1.SetTextColor(2)
st1H_1.SetX1NDC(0.65)
st1H_1.SetX2NDC(0.80)
st1H_1.SetY1NDC(0.70)
st1H_1.SetY2NDC(0.85)

h2H_1.SetLineColor(4)
h2H_1.SetFillColor(4)
h2H_1.SetFillStyle(3002)
h2H_1.Draw("goff")
gPad.Update()
h2H_1.SetName("Cluster Size == 2")
h2H_1.Draw("goff")
gPad.Update()
st2H_1 = h2H_1.GetListOfFunctions().FindObject("stats")
st2H_1.SetLineColor(4)
st2H_1.SetTextColor(4)
st2H_1.SetX1NDC(0.65)
st2H_1.SetX2NDC(0.80)
st2H_1.SetY1NDC(0.50)
st2H_1.SetY2NDC(0.65)

h3H_1.SetLineColor(6)
h3H_1.SetFillColor(6)
h3H_1.SetFillStyle(3002)
h3H_1.Draw("goff")
gPad.Update()
h3H_1.SetName("Cluster Size == 3")
h3H_1.Draw("goff")
gPad.Update()
st3H_1 = h3H_1.GetListOfFunctions().FindObject("stats")
st3H_1.SetLineColor(6)
st3H_1.SetTextColor(6)
st3H_1.SetX1NDC(0.65)
st3H_1.SetX2NDC(0.80)
st3H_1.SetY1NDC(0.30)
st3H_1.SetY2NDC(0.45)

h0H_1.Draw()
h1H_1.Draw("same")
h2H_1.Draw("same")
h3H_1.Draw("same")

c1V = TCanvas("c1V", "c1V")
c1V.SetLogy()

h0V_1.SetLineColor(8)
h0V_1.SetFillColor(8)
h0V_1.SetFillStyle(1001)
h0V_1.Draw("goff")
gPad.Update()
h0V_1.SetName("Cluster Size > 0 ")
st0V_1 = h0V_1.GetListOfFunctions().FindObject("stats")
st0V_1.SetLineColor(8)
st0V_1.SetTextColor(8)
st0V_1.SetX1NDC(0.45)
st0V_1.SetX2NDC(0.60)
st0V_1.SetY1NDC(0.70)
st0V_1.SetY2NDC(0.85)

h1V_1.SetLineColor(2)
h1V_1.SetFillColor(2)
h1V_1.SetFillStyle(3001)
h1V_1.Draw("goff")
gPad.Update()
h1V_1.SetName("Cluster Size == 1")
h1V_1.Draw("goff")
gPad.Update()
st1V_1 = h1V_1.GetListOfFunctions().FindObject("stats")
st1V_1.SetLineColor(2)
st1V_1.SetTextColor(2)
st1V_1.SetX1NDC(0.65)
st1V_1.SetX2NDC(0.80)
st1V_1.SetY1NDC(0.70)
st1V_1.SetY2NDC(0.85)

h2V_1.SetLineColor(4)
h2V_1.SetFillColor(4)
h2V_1.SetFillStyle(3002)
h2V_1.Draw("goff")
gPad.Update()
h2V_1.SetName("Cluster Size == 2")
h2V_1.Draw("goff")
gPad.Update()
st2V_1 = h2V_1.GetListOfFunctions().FindObject("stats")
st2V_1.SetLineColor(4)
st2V_1.SetTextColor(4)
st2V_1.SetX1NDC(0.65)
st2V_1.SetX2NDC(0.80)
st2V_1.SetY1NDC(0.50)
st2V_1.SetY2NDC(0.65)

h3V_1.SetLineColor(6)
h3V_1.SetFillColor(6)
h3V_1.SetFillStyle(3002)
h3V_1.Draw("goff")
gPad.Update()
h3V_1.SetName("Cluster Size == 3")
h3V_1.Draw("goff")
gPad.Update()
st3V_1 = h3V_1.GetListOfFunctions().FindObject("stats")
st3V_1.SetLineColor(6)
st3V_1.SetTextColor(6)
st3V_1.SetX1NDC(0.65)
st3V_1.SetX2NDC(0.80)
st3V_1.SetY1NDC(0.30)
st3V_1.SetY2NDC(0.45)

h0V_1.Draw()
h1V_1.Draw("same")
h2V_1.Draw("same")
h3V_1.Draw("same")


ctr = TCanvas("ctr", "ctr")
ctr.SetBottomMargin(0.15)
ctr.SetGridx()
ctr.SetGridy()
h1trig.GetYaxis().SetRangeUser(0,1.025)
h1trig.GetXaxis().SetTitleOffset(2)
h1trig.GetXaxis().CenterTitle()
h1trig.GetYaxis().CenterTitle()

h1trig.Draw()

f1gau = TF1("f1gau", "gaus(0)", -10., 10.)
f1gau.SetParameter(0, 1000.)
f1gau.SetParameter(1, 0.)
f1gau.SetParameter(2, 1.)

f2gau = TF1("f2gau", "gaus(0)+gaus(3)", -10., 10.)
f2gau.SetParameter(0, 1000.)
f2gau.SetParameter(1, 0.)
f2gau.SetParameter(2, 1.)
f2gau.SetParameter(3, 10.)
f2gau.SetParameter(4, 0.)
f2gau.SetParameter(5, 5.)

C4=TCanvas("C4","C4")

yymax = h1time.GetMaximum()*2.
f2gau.SetParLimits(0, 10., yymax)
f2gau.SetParLimits(1, -1.5, 1.5)
f2gau.SetParLimits(2, 0., 3.)

f2gau.SetParLimits(3, 0., yymax/5.)
f2gau.SetParLimits(4, -1.5, 1.5)
f2gau.SetParLimits(5, 1., 10.)

h1time.Fit(f2gau)

C4.Update()

C5=TCanvas("C5","C5")

h1timeMin.Fit(f2gau)

C5.Update()

raw_input("press enter")



