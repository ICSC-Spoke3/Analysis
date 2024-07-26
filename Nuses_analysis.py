
import pylab
import math
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import gridspec
import matplotlib.colors as colors
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
from matplotlib.axes import Axes
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import os
import sys
import time
from ROOT import gStyle, TFile, TTree, TChain, TH1D, TMath, TCanvas, TH1F, TH2F, gDirectory, TColor
import array as ary
import numpy as np
import glob

def NormalizeX(hh):   #Get pointer to TH2F normalized in X slices

    hname = "{0}_normX".format(hh.GetName())
    '''dir = gDirectory.Get(hname)
    print "@@@@@@@ dir" ,dir
    if(dir.Get(hname)):
        print ("NormalizeX::Deleting old {0}".format(str(hname)))
        dir.Get(str(hname)).Delete()'''

    hNorm = TH2F(hh).Clone(hname)
    hNorm.Sumw2()

    for ibinx in range (1, hNorm.GetNbinsX()+1):
        nentries=0
        for ibiny in range (1, hNorm.GetNbinsY()+1):
            nentries+=hNorm.GetBinContent(ibinx,ibiny)
        
        for ibiny in range (1, hNorm.GetNbinsY()+1):
            if( nentries<=0 ): continue
            value = hNorm.GetBinContent(ibinx,ibiny)
            value /= nentries
            hNorm.SetBinContent(ibinx,ibiny,value)
        
    return hNorm

if __name__ == '__main__':
    
    list = sys.argv[1]

    eventNo = []
    particle = []
    gEnergy = []
    gEnergy_um = []
    LayNo = []
    gap_mm = []
    fiber_yield = 8 #ph/keV
    fiber_trapp_eff = 0.054
    SiPM_PDE = 0.4

    bx = 991
    by = 100
    xxMin = 95
    yyMin = -0.5
    xxMax = 10000
    yyMax = 100.5

    h6 = TH2F("h6","Layer 0: ViewY", bx, xxMin, xxMax, by, yyMin, yyMax ) #vista 0 = vista Y
    h6.GetXaxis().SetTitle("Electron energy (keV)")
    h6.GetYaxis().SetTitle("Total Photo Electrons")
    h7 = TH2F("h7","Layer 1: ViewY", bx, xxMin, xxMax, by, yyMin, yyMax ) 
    h7.GetXaxis().SetTitle("Electron energy (keV)")
    h7.GetYaxis().SetTitle("Total Photo Electrons")
    h8 = TH2F("h8","Layer 2: ViewY", bx, xxMin, xxMax, by, yyMin, yyMax )
    h8.GetXaxis().SetTitle("Electron energy (keV)")
    h8.GetYaxis().SetTitle("Total Photo Electrons")
    h9 = TH2F("h9","Layer 0: ViewX", bx, xxMin, xxMax, by, yyMin, yyMax ) 
    h9.GetXaxis().SetTitle("Electron energy (keV)")
    h9.GetYaxis().SetTitle("Total Photo Electrons")
    h10 = TH2F("h10","Layer 1: ViewX", bx, xxMin, xxMax, by, yyMin, yyMax )
    h10.GetXaxis().SetTitle("Electron energy (keV)")
    h10.GetYaxis().SetTitle("Total Photo Electrons")
    h11 = TH2F("h11","Layer 2: ViewX", bx, xxMin, xxMax, by, yyMin, yyMax )
    h11.GetXaxis().SetTitle("Electron energy (keV)")
    h11.GetYaxis().SetTitle("Total Photo Electrons")

    with open(list) as fptr:
        fnamedata = fptr.readlines()
        #print ("fname: ", fnamedata)

        for jf in range(len(fnamedata)):
            rootfile = fnamedata[jf].split('\n')[0]
            print ("[INFO] ROOT file\t",rootfile)
            #print (rootfile.split('_')[5])
            eventNo.append(int((rootfile.split('_')[5]).split('-')[0]))
            particle.append((rootfile.split('_')[1]).split('/')[4])
            gEnergy.append(float((rootfile.split('_')[2]).split('-')[0]))
            gEnergy_um.append(((rootfile.split('_')[2]).split('-')[1]))
            LayNo.append(int((rootfile.split('_')[3]).split('-')[0]))
            gap_mm.append(int((rootfile.split('_')[4]).split('-')[0]))
            #print (gEnergy)
            #print (gEnergy_um)
            #print (LayNo)
            #print (gap_mm)
            #print (particle)
            #print (eventNo)
        
            tfile = TFile(rootfile)
            tree = tfile.Get("TrackerDigi")

            gg = 0
            if gEnergy_um[jf] == "keV":
                gg = 1
            elif gEnergy_um[jf] == "MeV":
                gg = 1000
            else:
                gg = -1            

            ##histo
            Emaxy = gEnergy[jf]*gg

            binx_l = 3
            binx_f = 400
            biny = 10*Emaxy
            minx_l = 0
            minx_f = 0
            miny = 0
            maxx_l = 3
            maxx_f = 400
            maxy = 2*Emaxy+0.5

            '''h1 = TH2F("h1","all Layer", binx_l, minx_l, maxx_l, biny, miny, maxy)
            h1.GetXaxis().SetTitle("Layer number")
            h1.GetYaxis().SetTitle("Edep [keV]")
            #h1.GetYaxis().SetRange(0, 30000)
            h2 = TH2F("h2","all Fiber", binx_f, minx_f, maxx_f, biny, miny, maxy)
            h2.GetXaxis().SetTitle("Fiber number")
            h2.GetYaxis().SetTitle("Edep [keV]")
            h3 = TH2F("h3","LayerNo=0", binx_f, minx_f, maxx_f, biny, miny, maxy)
            h3.GetXaxis().SetTitle("Fiber number")
            h3.GetYaxis().SetTitle("Edep [keV]")
            h4 = TH2F("h4","LayerNo=1", binx_f, minx_f, maxx_f, biny, miny, maxy)
            h4.GetXaxis().SetTitle("Fiber number")
            h4.GetYaxis().SetTitle("Edep [keV]")
            h5 = TH2F("h5","LayerNo=2", binx_f, minx_f, maxx_f, biny, miny, maxy)
            h5.GetXaxis().SetTitle("Fiber number")
            h5.GetYaxis().SetTitle("Edep [keV]")'''

            Nentries = tree.GetEntries()
            print ("------------fatto GEtEntries()")
            
            Edep = []
            LayerNo = []
            ViewNo = []
            FiberNo = []
            photoeleNo = []
            tot_phelectron_0_0 = np.zeros(Nentries) #vista 0 e layer 0
            tot_phelectron_0_1 = np.zeros(Nentries)
            tot_phelectron_0_2 = np.zeros(Nentries)
            tot_phelectron_1_0 = np.zeros(Nentries) #vista 1 e layer 0
            tot_phelectron_1_1 = np.zeros(Nentries)
            tot_phelectron_1_2 = np.zeros(Nentries)

            for i in range(Nentries):
                tree.GetEntry(i)

                Edep.append(tree.Energy_keV)        
                LayerNo.append(tree.LayerNo)
                ViewNo.append(tree.ViewNo)
                #print np.array(ViewNo[i])
                FiberNo.append(tree.FiberNo)                

                photoeleNo.append([])

                for j in range(len(ViewNo[i])):

                    photoeleNo[i].append(Edep[i][j]*fiber_yield*fiber_trapp_eff*SiPM_PDE)

                    #print "V ", ViewNo[i][j], " (elemento ", j, ")"
                    #print "L ", LayerNo[i][j], " (elemento ", j, ")"
                    #print "photoeleNo ", np.array(photoeleNo)

                    if(ViewNo[i][j]==0 and LayerNo[i][j]==0):
                        tot_phelectron_0_0[i] += photoeleNo[i][j]
                        #print "gEnergy ", Emaxy, " keV      Photele00 ", tot_phelectron_0_0[i][j]
                    if(ViewNo[i][j]==0 and LayerNo[i][j]==1):
                        tot_phelectron_0_1[i] += photoeleNo[i][j]                        
                        #print "gEnergy ", Emaxy, " keV      Photele01 ", tot_phelectron_0_1[i][j]
                    if(ViewNo[i][j]==0 and LayerNo[i][j]==2):
                        tot_phelectron_0_2[i] += photoeleNo[i][j]                        
                        #print "gEnergy ", Emaxy, " keV      Photele02 ", tot_phelectron_0_2[i][j]
                    if(ViewNo[i][j]==1 and LayerNo[i][j]==0):
                        tot_phelectron_1_0[i] += photoeleNo[i][j]                        
                        #print "gEnergy ", Emaxy, " keV      Photele10 ", tot_phelectron_1_0[i][j]
                    if(ViewNo[i][j]==1 and LayerNo[i][j]==1):
                        tot_phelectron_1_1[i] += photoeleNo[i][j]                        
                        #print "gEnergy ", Emaxy, " keV      Photele11 ", tot_phelectron_1_1[i][j]
                    if(ViewNo[i][j]==1 and LayerNo[i][j]==2):
                        tot_phelectron_1_2[i] += photoeleNo[i][j]                        
                        #print "gEnergy ", Emaxy, " keV      Photele12 ", tot_phelectron_1_2[i][j]

                h6.Fill(Emaxy, tot_phelectron_0_0[i])
                h7.Fill(Emaxy, tot_phelectron_0_1[i])
                h8.Fill(Emaxy, tot_phelectron_0_2[i])
                h9.Fill(Emaxy, tot_phelectron_1_0[i])
                h10.Fill(Emaxy, tot_phelectron_1_1[i])
                h11.Fill(Emaxy, tot_phelectron_1_2[i])            
            
            print ("Ele_energy = ", gEnergy[jf], " ", gEnergy_um[jf], ":")
            print (" Tot_Ph_00 = ", sum(tot_phelectron_0_0))
            print (" Tot_Ph_01 = ", sum(tot_phelectron_0_1))
            print (" Tot_Ph_02 = ", sum(tot_phelectron_0_2))
            print (" Tot_Ph_10 = ", sum(tot_phelectron_1_0))
            print (" Tot_Ph_11 = ", sum(tot_phelectron_1_1))
            print (" Tot_Ph_12 = ", sum(tot_phelectron_1_2))
            
            '''canv = TCanvas ("c_{0}".format(gEnergy[jf]) ,"canv_{0}".format(gEnergy[jf]))
            canv_fib = TCanvas ("c_fib_{0}".format(gEnergy[jf]) ,"canv_{0}".format(gEnergy[jf]))
    
            canv.Divide(2,1)
            canv.cd(1)
            tree.Draw("Energy_keV:LayerNo >> h1", "", "COLZ")
            canv.cd(2)
            tree.Draw("Energy_keV:FiberNo >> h2", "", "COLZ")
            
            #canv.Print( "./Plots/NUSES/E_vs_layNo_{0}.png".format(gEnergy[jf]) )

            canv.Modified()
            canv.Update()

            canv_fib.Divide(2,2)
            canv_fib.cd(1)
            tree.Draw("Energy_keV:FiberNo >> h3", "LayerNo==0", "COLZ")
            canv_fib.cd(2)
            tree.Draw("Energy_keV:FiberNo >> h4", "LayerNo==1", "COLZ")
            canv_fib.cd(3)
            tree.Draw("Energy_keV:FiberNo >> h5", "LayerNo==2", "COLZ")

            fout.cd()
            #h1.Write()
            #h2.Write()
            canv.Write()
            canv_fib.Write()'''

            
    
    fout = TFile("./Plots/NUSES_Results_{0}evt.root".format(Nentries),"recreate")
    print ("Create output file: ", fout)

    gStyle.SetPalette(55)
    
    canv_photo_0 = TCanvas("c_photo_0" ,"c_photo_0") #layer0
    canv_photo_1 = TCanvas("c_photo_1" ,"c_photo_1")
    canv_photo_2 = TCanvas("c_photo_2" ,"c_photo_2")
    canv_photo_0.Divide(1,2)
    canv_photo_0.cd(1)
    #canv_photo_0.cd(1).SetLogz()
    h6.GetZaxis().SetRangeUser(0, 80)
    h6.SetContour(100)
    h6.Draw("colz")
    '''canv_photo_0.cd(2)
    h6_norm=NormalizeX(h6)
    h6_norm.Draw("COLZ")'''
    canv_photo_0.cd(2)
    h9.GetZaxis().SetRangeUser(0, 80)
    h9.Draw("COLZ")
    '''canv_photo_0.cd(4)
    h9_norm=NormalizeX(h9)
    h9_norm.Draw("COLZ")'''

    canv_photo_1.Divide(1,2)
    canv_photo_1.cd(1)
    h7.GetZaxis().SetRangeUser(0, 80)
    h7.Draw("COLZ")
    '''canv_photo_1.cd(2)
    h7_norm=NormalizeX(h7)
    h7_norm.Draw("COLZ")'''
    canv_photo_1.cd(2)
    h10.GetZaxis().SetRangeUser(0, 80)
    h10.Draw("COLZ")
    '''canv_photo_1.cd(4)
    h10_norm=NormalizeX(h10)
    h10_norm.Draw("COLZ")'''

    canv_photo_2.Divide(1,2)
    canv_photo_2.cd(1)
    h8.GetZaxis().SetRangeUser(0, 80)
    h8.Draw("COLZ")
    '''canv_photo_2.cd(2)
    h8_norm=NormalizeX(h8)
    h8_norm.Draw("COLZ")'''
    canv_photo_2.cd(2)
    h11.GetZaxis().SetRangeUser(0, 80)
    h11.Draw("COLZ")
    '''canv_photo_2.cd(4)
    h11_norm=NormalizeX(h11)
    h11_norm.Draw("COLZ")'''

    fout.cd()
    canv_photo_0.Write("Layer1")
    canv_photo_1.Write("Layer2")
    canv_photo_2.Write("Layer3")

    fout.Close()

            


    '''c1 = TCanvas ("c1" ,"c1", 50,50,800,800)
    c1.Divide(2,1)
    c1.cd(1)
    #hh_lay.Draw("COLTZ")
    tree.Draw("Energy_keV:LayerNo >> hh_lay", "", "COLTZ")
    c1.cd(2)
    #hh_fib.Draw("COLTZ")
    tree.Draw("Energy_keV:FiberNo >> hh_fib","","COLTZ")
    c1.Modified()
    c1.Update()
    c1.Print( "prova.png" )
    #c1.SaveAs("lala.png")

    hh_lay.SetDirectory(0)
    hh_fib.SetDirectory(0)
    tfile.Close ()'''


    

    
    

