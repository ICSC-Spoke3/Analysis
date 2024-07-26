import sys, os

import numpy as np

import array as ary

from ROOT import *

ran = TRandom3()


# 8 ph/keV, trapping efficiency=5.4%, PDE = 0.4
phkev = 8.*0.92 # including cladding
trap = 5.4e-2
PDE = 0.4
nphthr = 4.

m = 511. # keV
nfiblay = 64*4
nfl = 2
frin = "FibModule_500um_2ribbon00*_Out_Ele.root"
#nfl = 3
#frin = "FibModule_500um_3ribbon00*_Out_Ele.root"

Ss = 250.e-4 # SiPM strip pitch
Rf = 250.e-4 # fiber radius

Sc = 0.5*Rf # middle fiber center wrt the middle strip center 
Hr = Rf - 0.5*Ss + Sc 
Hl = Rf - 0.5*Ss - Sc 

Ar = Rf*Rf*TMath.ACos(1.-Hr/Rf) - (Rf - Hr)*TMath.Sqrt( Rf*Rf - TMath.Power( Rf-Hr, 2.))
Al = Rf*Rf*TMath.ACos(1.-Hl/Rf) - (Rf - Hl)*TMath.Sqrt( Rf*Rf - TMath.Power( Rf-Hl, 2.))

Af = TMath.Pi()*Rf*Rf

Fr = Ar/Af
Fl = Al/Af
Fm = 1. - Fr - Fl

print "Fm = ", Fm
print "Fl = ", Fl
print "Fr = ", Fr



print "Opening ", frin

te = TChain("CFTree")
th = TChain("CFHits")

te.Add(frin)
th.Add(frin)

Nent = te.GetEntries()
print "Tree entries = ", Nent

#Nent = 10000

nemc = 38 #*2
e0mc = 100.
e1mc = 2000.

hPdEdX = TProfile("hdEdX", ";MC Energy (keV); Average Energy loss (MeV/cm)", nemc, e0mc, e1mc)
hPdX = TProfile("hdX", ";MC Energy (keV); Average Track lenght (cm)", nemc, e0mc, e1mc)

for i in range(Nent):
    th.GetEntry(i)
    nh = th.Nhit

    #print "nh = ", nh

    mreg = th.Mreg
    etra = th.Etra
    dtra = th.Dtra
    ctra = th.Ctra
    xtra = th.Xtra
    ytra = th.Ytra
    ztra = th.Ztra
    elos = th.Elos


    oldreg = -1
    ttra = 0.
    telos = 0.
    for j in range(nh):
        if(mreg[j] != oldreg):
            #print "===> oldreg, ttra ", oldreg, mreg[j], xtra[j], ytra[j], ztra[j], ttra, telos
            if(ttra>0.):
                #print "---> dedx = ", telos/ttra
                hPdEdX.Fill(etra[j]*1.e6-m, telos/ttra)
                hPdX.Fill(etra[j]*1.e6-m, ttra)

            ttra = 0.
            telos = 0.
        else:
            ttra += ctra[j]
            telos += elos[j]*1.e3
        oldreg = mreg[j]
    #print "oldreg, ttra ", oldreg, mreg[j], xtra[j], ytra[j], ztra[j], ttra, telos
    if(ttra>0.):
        #print "dedx = ", telos/ttra
        hPdEdX.Fill(etra[j]*1.e6-m, telos/ttra)
        hPdX.Fill(etra[j]*1.e6-m, ttra)
  

gStyle.SetOptStat(0)

c00 = TCanvas("c00","c00")
hPdX.Draw()

c03 = TCanvas("c03","c03")
hPdEdX.Draw()


h1NMC = TH1F("h1NMC", ";MC Energy (keV); Number of Events", nemc, e0mc, e1mc)
h1NTrigX = TH1F("h1NTrigX", ";MC Energy (keV); Number of Trigger X-view", nemc, e0mc, e1mc)
h1NTrigY = TH1F("h1NTrigY", ";MC Energy (keV); Number of Trigger Y-view", nemc, e0mc, e1mc)

h1PTrigX = TH1F("h1PTrigX", ";MC Energy (keV); Trigger Efficiency", nemc, e0mc, e1mc)
h1PTrigY = TH1F("h1PTrigY", ";MC Energy (keV); Trigger Efficiency", nemc, e0mc, e1mc)

h1NTrigX1 = TH1F("h1NTrigX1", ";MC Energy (keV); Number of Trigger X-view", nemc, e0mc, e1mc)
h1NTrigY1 = TH1F("h1NTrigY1", ";MC Energy (keV); Number of Trigger Y-view", nemc, e0mc, e1mc)

h1PTrigX1 = TH1F("h1PTrigX1", ";MC Energy (keV); Trigger Efficiency", nemc, e0mc, e1mc)
h1PTrigY1 = TH1F("h1PTrigY1", ";MC Energy (keV); Trigger Efficiency", nemc, e0mc, e1mc)

hPNphX = TProfile("hPNphX", ";MC Energy (keV); Average p.e.; Number of Events", nemc, e0mc, e1mc, "S")
hPNphy = TProfile("hPNphY", ";MC Energy (keV); Average p.e.; Number of Events", nemc, e0mc, e1mc, "S")

h2NphX = TH2F("h2NphX", ";MC Energy (keV); Average p.e.", nemc, e0mc, e1mc, 101, -0.5, 100.)
h2Nphy = TH2F("h2NphY", ";MC Energy (keV); Average p.e.", nemc, e0mc, e1mc, 101, -0.5, 100.)

hPElosX = TProfile("hPElosX", ";MC Energy (keV); Average Energy loss (keV)", nemc, e0mc, e1mc, "S")
hPElosY = TProfile("hPElosY", ";MC Energy (keV); Average Energy loss (keV)", nemc, e0mc, e1mc, "S")

h2EtEfib_e_x = TH2F("heEtEfib_e_x", ";MC Energy (keV); Fib Energy (keV); Number of Events", nemc, e0mc, e1mc, 200, 0., 1000.) 
h2EtEfib_e_y = TH2F("heEtEfib_e_y", ";MC Energy (keV); Fib Energy (keV); Number of Events", nemc, e0mc, e1mc, 200, 0., 1000.) 


for i in range(Nent):
    te.GetEntry(i)
    ek = te.Ekevt*1.e6
    h1NMC.Fill(ek)
    # print i, te.Xevt*1.e4, te.Yevt*1.e4
    Nfib = te.NFib
    Efib = te.EFib
    Idfib = te.IDFib
    etfibx = 0.
    etfiby = 0.
    nsipmx = 0
    qsipmx = np.zeros((nfiblay+1)*2)
    idsipmx = []
    nsipmy = 0
    qsipmy = np.zeros((nfiblay+1)*2)
    idsipmy = []

    nsipmx1 = 0
    qsipmx1 = np.zeros((nfiblay+1)*2)
    idsipmx1 = []
    nsipmy1 = 0
    qsipmy1 = np.zeros((nfiblay+1)*2)
    idsipmy1 = []
    for j in range(Nfib):
        k = Idfib[j]
        kk = k%nfiblay
        if(k<nfiblay*nfl):
            etfibx += Efib[j]*1.e6
            kkk = int(k/nfiblay)
            isipm = kk*2 + 1 - kkk
            if(qsipmx[isipm]==0.):
                nsipmx += 1
                idsipmx.append(isipm)
                qsipmx[isipm] = Efib[j]*1.e6/2.
            else:
                qsipmx[isipm] += Efib[j]*1.e6/2.
            
            isipml = isipm - 1
            if(isipml >= 0):
                if(qsipmx1[isipml]==0.):
                    nsipmx1 += 1
                    idsipmx1.append(isipml)
                    qsipmx1[isipml] = Efib[j]*1.e6*Fl
                else:
                    qsipmx1[isipml] += Efib[j]*1.e6*Fl

            if(qsipmx1[isipm]==0.):
                nsipmx1 += 1
                idsipmx1.append(isipm)
                qsipmx1[isipm] = Efib[j]*1.e6*Fm
            else:
                qsipmx1[isipm] += Efib[j]*1.e6*Fm

            isipmr = isipm + 1
            if(isipmr < (nfiblay+1)*2):
                if(qsipmx1[isipmr]==0.):
                    nsipmx1 += 1
                    idsipmx1.append(isipmr)
                    qsipmx1[isipmr] = Efib[j]*1.e6*Fr
                else:
                    qsipmx1[isipmr] += Efib[j]*1.e6*Fr
               
            isipm += 1
            if(qsipmx[isipm]==0.):
                nsipmx += 1
                idsipmx.append(isipm)
                qsipmx[isipm] = Efib[j]*1.e6/2.
            else:
                qsipmx[isipm] += Efib[j]*1.e6/2.
        else:
            etfiby += Efib[j]*1.e6
            kkk = int((k-nfl*nfiblay)/nfiblay)
            isipm = kk*2 + 1 - kkk
            if(qsipmy[isipm]==0.):
                nsipmy += 1
                idsipmy.append(isipm)
                qsipmy[isipm] = Efib[j]*1.e6/2.
            else:
                qsipmy[isipm] += Efib[j]*1.e6/2.

            isipml = isipm - 1
            if(isipml >= 0):
                if(qsipmy1[isipml]==0.):
                    nsipmy1 += 1
                    idsipmy1.append(isipml)
                    qsipmy1[isipml] = Efib[j]*1.e6*Fl
                else:
                    qsipmy1[isipml] += Efib[j]*1.e6*Fl

            if(qsipmy1[isipm]==0.):
                nsipmy1 += 1
                idsipmy1.append(isipm)
                qsipmy1[isipm] = Efib[j]*1.e6*Fm
            else:
                qsipmy1[isipm] += Efib[j]*1.e6*Fm

            isipmr = isipm + 1
            if(isipmr < (nfiblay+1)*2):
                if(qsipmy1[isipmr]==0.):
                    nsipmy1 += 1
                    idsipmy1.append(isipmr)
                    qsipmy1[isipmr] = Efib[j]*1.e6*Fr
                else:
                    qsipmy1[isipmr] += Efib[j]*1.e6*Fr
               

            isipm += 1
            if(qsipmy[isipm]==0.):
                nsipmy += 1
                idsipmy.append(isipm)
                qsipmy[isipm] = Efib[j]*1.e6/2.
            else:
                qsipmy[isipm] += Efib[j]*1.e6/2.
            
        #print "### ", j, k, kk, kkk, isipm-1

    idsipmx.sort()
    idsipmy.sort()
    #print "@@ ", nsipmx, idsipmx
    etfibx1 = 0.
    for l in idsipmx:
        #print "---> ", qsipmx[l]
        etfibx1 += qsipmx[l]
        
    #print "  ", etfibx, etfibx1
    #print "$$ ", nsipmy, idsipmy
    etfiby1 = 0.
    for l in idsipmy:
        #print "---> ", qsipmy[l]
        etfiby1 += qsipmy[l]
    #print "  ", etfiby, etfiby1
    
    h2EtEfib_e_x.Fill(ek, etfibx)
    h2EtEfib_e_y.Fill(ek, etfiby)
    hPElosX.Fill(ek, etfibx)  
    hPElosY.Fill(ek, etfiby)
    
    muphx = etfibx*phkev*trap*PDE 
    muphy = etfiby*phkev*trap*PDE 
    hPNphX.Fill(ek, muphx)  
    hPNphY.Fill(ek, muphy)  
    h2NphX.Fill(ek, muphx)  
    h2NphY.Fill(ek, muphy)  
    
    ftrigx = 0
    if(nsipmx>1):
        for l in range(nsipmx):
            ll = idsipmx[l]
            muph0 = qsipmx[ll]*phkev*trap*PDE 
            muph1 = qsipmx[ll+1]*phkev*trap*PDE 
            nph0 = ran.Poisson(muph0)
            nph1 = ran.Poisson(muph1)
            if( (nph0 >= nphthr) and (nph1 >= nphthr) ):
                ftrigx = 1
            # print "% ", l, muph0, muph1, nph0, nph1, ftrigx
    if( ftrigx > 0 ):
        h1NTrigX.Fill(ek)
    
    ftrigy = 0
    if(nsipmy>1):
        for l in range(nsipmy):
            ll = idsipmy[l]
            muph0 = qsipmy[ll]*phkev*trap*PDE 
            muph1 = qsipmy[ll+1]*phkev*trap*PDE 
            nph0 = ran.Poisson(muph0)
            nph1 = ran.Poisson(muph1)
            if( (nph0 >= nphthr) and (nph1 >= nphthr) ):
                ftrigy = 1
            # print "& ", l, muph0, muph1, nph0, nph1, ftrigy

    if( ftrigy > 0 ):
        h1NTrigY.Fill(ek)

    idsipmx1.sort()
    idsipmy1.sort()
    ftrigx = 0
    if(nsipmx1>1):
        for l in range(nsipmx1):
            ll = idsipmx1[l]
            muph0 = qsipmx1[ll]*phkev*trap*PDE 
            muph1 = qsipmx1[ll+1]*phkev*trap*PDE 
            nph0 = ran.Poisson(muph0)
            nph1 = ran.Poisson(muph1)
            if( (nph0 >= nphthr) and (nph1 >= nphthr) ):
                ftrigx = 1
            # print "% ", l, muph0, muph1, nph0, nph1, ftrigx
    if( ftrigx > 0 ):
        h1NTrigX1.Fill(ek)
    
    ftrigy = 0
    if(nsipmy>1):
        for l in range(nsipmy1):
            ll = idsipmy1[l]
            muph0 = qsipmy1[ll]*phkev*trap*PDE 
            muph1 = qsipmy1[ll+1]*phkev*trap*PDE 
            nph0 = ran.Poisson(muph0)
            nph1 = ran.Poisson(muph1)
            if( (nph0 >= nphthr) and (nph1 >= nphthr) ):
                ftrigy = 1
            # print "& ", l, muph0, muph1, nph0, nph1, ftrigy

    if( ftrigy > 0 ):
        h1NTrigY1.Fill(ek)
               
for i in range(h1NTrigX.GetNbinsX()):
    ngen = h1NMC.GetBinContent(i+1)
    if(ngen<=0.):
        continue
    xdx = h1NTrigX.GetBinContent(i+1)
    effx = xdx/ngen
    h1PTrigX.SetBinContent(i+1, effx)
    errx = TMath.Sqrt(effx*(1.-effx)/ngen)
    h1PTrigX.SetBinError(i+1, errx)
    
    xdy = h1NTrigY.GetBinContent(i+1)
    effy = xdy/ngen
    h1PTrigY.SetBinContent(i+1, effy)
    erry = TMath.Sqrt(effy*(1.-effy)/ngen)
    h1PTrigY.SetBinError(i+1, erry)
    # print "Eff : ", i, xdx, xdy, effx, effy

    xdx = h1NTrigX1.GetBinContent(i+1)
    effx = xdx/ngen
    h1PTrigX1.SetBinContent(i+1, effx)
    errx = TMath.Sqrt(effx*(1.-effx)/ngen)
    h1PTrigX1.SetBinError(i+1, errx)
    
    xdy = h1NTrigY1.GetBinContent(i+1)
    effy = xdy/ngen
    h1PTrigY1.SetBinContent(i+1, effy)
    erry = TMath.Sqrt(effy*(1.-effy)/ngen)
    h1PTrigY1.SetBinError(i+1, erry)




c3 = TCanvas("c3","c3")
c3.SetGridx()
c3.SetGridy()
hPElosX.SetMarkerStyle(20)
hPElosY.SetMarkerStyle(20)
hPElosX.SetMarkerColor(1)
hPElosY.SetMarkerColor(2)
hPElosX.SetLineColor(1)
hPElosY.SetLineColor(2)
hPElosY.SetFillColor(2)
hPElosY.GetYaxis().SetRangeUser(0., 700)
hPElosY.Draw("e1")
hPElosX.Draw("e1same")


c0 = TCanvas("c0","c0")
c0.SetGridx()
c0.SetGridy()
hPNphX.SetMarkerStyle(20)
hPNphY.SetMarkerStyle(20)
hPNphX.SetMarkerColor(1)
hPNphY.SetMarkerColor(2)
hPNphX.SetLineColor(1)
hPNphY.SetLineColor(2)
hPNphY.GetYaxis().SetRangeUser(0., 100)
hPNphY.Draw("e1")
hPNphX.Draw("e1same")

c20x = TCanvas("c20x","c20x")
c20x.SetGridx()
c20x.SetGridy()
h2NphX.Draw("colz")

c20y = TCanvas("c20y","c20y")
c20y.SetGridx()
c20y.SetGridy()
h2NphY.Draw("colz")


c1 = TCanvas("c1","c1")
c1.SetGridx()
c1.SetGridy()
h1PTrigX.SetMarkerStyle(20)
h1PTrigY.SetMarkerStyle(20)
h1PTrigX.SetMarkerColor(1)
h1PTrigY.SetMarkerColor(2)
h1PTrigX.SetLineColor(1)
h1PTrigY.SetLineColor(2)
h1PTrigY.GetYaxis().SetRangeUser(0., 1.1)
h1PTrigY.Draw("pe1")
h1PTrigX.Draw("pe1same")

h1PTrigX1.SetMarkerStyle(24)
h1PTrigY1.SetMarkerStyle(24)
h1PTrigX1.SetMarkerColor(1)
h1PTrigY1.SetMarkerColor(2)
h1PTrigX1.SetLineColor(1)
h1PTrigY1.SetLineColor(2)
h1PTrigY1.Draw("pe1same")
h1PTrigX1.Draw("pe1same")



c2 = TCanvas("c2","c2")
c2.SetRightMargin(0.15)
h2EtEfib_e_x.SetMarkerStyle(1)
h2EtEfib_e_y.SetMarkerStyle(1)
h2EtEfib_e_x.SetMarkerColor(1)
h2EtEfib_e_y.SetMarkerColor(2)
h2EtEfib_e_y.Draw("")
h2EtEfib_e_x.Draw("same")


raw_input("Next event ...")
    
