import os
import sys
import time
from ROOT import TFile, TTree, TChain, TH1D, TH2D, TMath, TCanvas, TH1F
import array as ary
import numpy as np
import glob

def Clustering(hits, charges):

    nhits = len(hits)
    maxclus = len(hits)
    
    '''istrip0 = np.zeros(maxclus,dtype = int)
    nstrip0 = np.zeros(maxclus,dtype = int)
    charge0 = np.zeros(maxclus)'''
    istrip0 = []
    nstrip0 = []
    charge0 = []
    
    istrip = 0
    nclu = 0
    istrip0.append(hits[istrip])
    charge0.append(charges[istrip])
    nstrip0.append(1)
    nclu += 1
    istrip += 1
    while (istrip < nhits):
        print (istrip)
        print (charge0)
        if ((hits[istrip]==hits[istrip-1]+1)):
            nstrip0[nclu-1] += 1
            charge0[nclu-1] += charges[istrip]
        else:
            
            istrip0.append(hits[istrip])
            nstrip0.append(1)
            charge0.append(charges[istrip])
            nclu += 1                
        istrip += 1
   

    istrip0 = istrip0[:nclu]
    nstrip0 = nstrip0[:nclu]
    charge0 = charge0[:nclu]

    return nclu,istrip0,nstrip0,charge0

if __name__ == '__main__':

    list = sys.argv[1]
    fout = TFile("analysisNUSES.root","recreate")

    totLayers = 3
    nfiber = 130
    fiberyield = 8*0.92#kev/pe
    PDE = 0.4
    trapeff = 0.054
        
    nChan = 1000 
    eventID = ary.array('i',[0])
    Tch = ary.array('i',[0])
    Tchs = ary.array('i',[0])
    ILay = ary.array('i',range(nChan))
    IView = ary.array('i',range(nChan))
    Ipix = ary.array('i',range(nChan))
    Istrip = ary.array('i',range(nChan))
    Npe = ary.array('f',range(nChan))
    NpeStrip = ary.array('f',range(nChan))
    Tns = ary.array('f',range(nChan))
    
    nClusters = 1000 
    Tclu = ary.array('i',[0])
    Tnpe = ary.array('f',range(nClusters))
    Csiz = ary.array('i',range(nClusters))
    Fstr = ary.array('i',range(nClusters))
    Mtim = ary.array('f',range(nClusters))

    TcluHor = ary.array('i',[0])
    CsizHor = ary.array('i',range(nClusters))
    TnpeHor = ary.array('f',range(nClusters))
    FstrHor = ary.array('i',range(nClusters))
    MtimHor = ary.array('f',range(nClusters))

    TcluVer = ary.array('i',[0])
    CsizVer = ary.array('i',range(nClusters))
    TnpeVer = ary.array('f',range(nClusters))
    FstrVer = ary.array('i',range(nClusters))
    MtimVer = ary.array('f',range(nClusters))
    
    OutTree    = TTree("ClusterTree", "ClusterTree")
    
    # Define output tree
    OutTree.Branch("eventID", eventID, "eventID/I")
    OutTree.Branch("Tch", Tch, "Tch/I")
    OutTree.Branch("Tchs", Tchs, "Tchs/I")
    OutTree.Branch("ILay", ILay, "ILay[Tch]/I")
    OutTree.Branch("IView", IView, "IView[Tch]/I")
    OutTree.Branch("Ipix", Ipix, "Ipix[Tch]/I")
    OutTree.Branch("Istrip", Istrip, "Istrip[Tchs]/I")
    OutTree.Branch("Npe", Npe, "Npe[Tch]/F")
    OutTree.Branch("NpeStrip", NpeStrip, "NpeStrip[Tchs]/F")
    

    # Informazioni generali sui cluster
    OutTree.Branch("Tclu", Tclu, "Tclu/I")
    OutTree.Branch("Csiz", Csiz, "Csiz[Tclu]/I")
    OutTree.Branch("Tnpe", Tnpe, "Tnpe[Tclu]/F")
    OutTree.Branch("Fstr", Fstr, "Fstr[Tclu]/I")
    
    # Informazioni sui cluster nelle viste orizzontale(x=1) e verticale(y=0)
    OutTree.Branch("TcluHor", TcluHor, "TcluHor/I")
    OutTree.Branch("CsizHor", CsizHor, "CsizHor[TcluHor]/I")
    OutTree.Branch("TnpeHor", TnpeHor, "TnpeHor[TcluHor]/F")
    OutTree.Branch("FstrHor", FstrHor, "FstrHor[TcluHor]/I")

    OutTree.Branch("TcluVer", TcluVer, "TcluVer/I")
    OutTree.Branch("CsizVer", CsizVer, "CsizVer[TcluVer]/I")
    OutTree.Branch("TnpeVer", TnpeVer, "TnpeVer[TcluVer]/F")
    OutTree.Branch("FstrVer", FstrVer, "FstrVer[TcluVer]/I")
    
    
    canvas = []
    
    hl0y = TH2D("hl0y","layer0 y-view",991,90,10000,200,-0.5,100) 
    hl0y.SetContour(10)
    hl0y.GetZaxis().SetRangeUser(0,80)

    hl0x = TH2D("hl0x","layer0 x-view",991,90,10000,200,-0.5,100) 
    hl0x.SetContour(20)
    hl0x.GetZaxis().SetRangeUser(0,80)


    hl1y = TH2D("hl1y","layer1 y-view",991,90,10000,200,-0.5,100) 
    hl1y.SetContour(30)
    hl1y.GetZaxis().SetRangeUser(0,80)

    hl1x = TH2D("hl1x","layer1 x-view",991,90,10000,200,-0.5,100) 
    hl1x.SetContour(40)
    hl1x.GetZaxis().SetRangeUser(0,80)

    hl2y = TH2D("hl2y","layer2 y-view",991,90,10000,200,-0.5,100) 
    hl2y.SetContour(50)
    hl2y.GetZaxis().SetRangeUser(0,80)
    hl2x = TH2D("hl2x","layer2 x-view",991,90,10000,200,-0.5,100) 
    hl2x.SetContour(60)
    hl2x.GetZaxis().SetRangeUser(0,80)

    for i in range(totLayers):
        canvas.append([])
        canvas[i]=TCanvas("c{0}".format(i),"c{0}".format(i))
        canvas[i].Divide(1,2)
        
    #gEnergy, NLayer, gap_mm, NevtSim, dE_vec,layer_vec,view_vec,fiber_vec,hl0y,hl0x,hl1y,hl1x,hl2y,hl2x = ReadList(list,hl0y,hl0x,hl1y, hl1x,hl2y,hl2x)

    #hl0y,hl0x,hl1y,hl1x,hl2y,hl2x = ReadList(list,hl0y,hl0x,hl1y,hl1x,hl2y,hl2x)

    gEnergy = []
    NLayer = []
    gap_mm = []
    NevtSim = []
    
    with open(list) as fptr:
        fnamedata = fptr.readlines()
        for jf in range(len(fnamedata)):
            
            #print ("File: ",jf)
            
            dE_vec = []
            pe_vec = []
            layer_vec = []
            view_vec = []
            fiber_vec = []
            pe_vec = []
            
            rootfile = fnamedata[jf].split()[0]
            print(fnamedata[jf])
            print (rootfile)

 
            print ("[INFO] ROOT file ",jf," ",rootfile)

            gEnergy.append(float((rootfile.split('_')[2]).split('-')[0]))
            #print(gEnergy)
            energy_mu = (rootfile.split('_')[2]).split('-')[1]
            #print(energy_mu)
            #NLayer.append(int((rootfile.split('_')[3]).split('-')[0]))
            #gap_mm.append(int((rootfile.split('_')[4]).split('-')[0]))
            #NevtSim.append(int((rootfile.split('_')[5]).split('-')[0]))
            
            if energy_mu == "keV":
                sfactor = 1
            elif energy_mu == "MeV":
                sfactor = 1000
            else:
                 sfactor = 0

            Emax  = gEnergy[jf]*sfactor
            print ("[INFO] electron energy",Emax, "keV")
            #nbinE = Emax
            #maxy = Emax +0.2* Emax
            
            tfile = TFile(rootfile)
            tree = tfile.Get("TrackerDigi")
            nhit = tree.GetEntries()
            nhit = 1

            
            totPEl0x = np.zeros(nhit)
            totPEl0y = np.zeros(nhit)
            totPEl1x = np.zeros(nhit)
            totPEl1y = np.zeros(nhit)
            totPEl2x = np.zeros(nhit)
            totPEl2y = np.zeros(nhit)
            
            for i in range(nhit):
                EventID = i
                pe_vec.append([])
                tree.GetEntry(i)     
                dE_vec.append(np.array(tree.Energy_keV))
                layer_vec.append(tree.LayerNo)
                view_vec.append(tree.ViewNo)
                fiber_vec.append(tree.FiberNo)
                
                PEstrip =[[[],[]],[[],[]],[[],[]]]
                for il in range(3):
                    for iv in range(2):
                        PEstrip[il][iv]=np.zeros(nfiber)                
                #print ("hit: ",len(dE_vec[i]))

                Lay = []
                View = []
                Fiber = []
                Strip_vec = []
                PEStrip_vec = []
                npe_vec = []
                Nclu = 0
                Fstrip = []
                Clusiz = []
                Clunpe = []
                Tch[0] = len(dE_vec[i])
                Tchs[0] = 19#Tch[0]*2
                NcluH = 0
                FstripH = []
                ClusizH = []
                ClunpeH = []
                NcluV = 0
                FstripV = []
                ClusizV = []
                ClunpeV = []
                for j in range(len(dE_vec[i])):
                    pe_vec[i].append((dE_vec[i][j]*fiberyield*trapeff*PDE))
               
                    if pe_vec[i][j] != 0:
                        l = layer_vec[i][j]
                        v = view_vec[i][j]
                        f = fiber_vec[i][j]
                        Lay.append(l)
                        View.append(v)
                        Fiber.append(f)
                        npe_vec.append(pe_vec[i][j])

                        if f == 0:
                            index1 = 1
                            index2 = 1
                        elif f == 129:
                            index1 = nfiber-1
                            index2 = nfiber-1
                        else:
                            index1 = f
                            index2 = f-1
                    
                        PEstrip[l][v][index1] += pe_vec[i][j]/2
                        PEstrip[l][v][index2] += pe_vec[i][j]/2
                        Strip_vec.append(index1)
                        Strip_vec.append(index2)
                        PEStrip_vec.append(pe_vec[i][j]/2)
                        PEStrip_vec.append(pe_vec[i][j]/2)

                        '''
                        print ("Layer: ",l )
                        print ("View: ",v )
                        print ("Fiber: ",f,"-> Strips: ",index1,index2)
                        print("PE:", pe_vec[i][j])
                        print("s: ",PEstrip[l][v])'''
                    
                    if (layer_vec[i][j] == 0 and view_vec[i][j] == 0):
                        totPEl0y[i] += pe_vec[i][j]   
                    if (layer_vec[i][j] == 0 and view_vec[i][j] == 1):
                        totPEl0x[i] += pe_vec[i][j]
                    if (layer_vec[i][j] == 1 and view_vec[i][j] == 0):
                        totPEl1y[i] += pe_vec[i][j]
                    if (layer_vec[i][j] ==1 and view_vec[i][j] == 1):
                        totPEl1x[i] += pe_vec[i][j]
                    if (layer_vec[i][j] == 2 and view_vec[i][j] == 0):
                        totPEl2y[i] += pe_vec[i][j]
                    if (layer_vec[i][j] ==2 and view_vec[i][j] == 1):
                        totPEl2x[i] += pe_vec[i][j]

                hl0y.Fill(float(Emax),float(totPEl0y[i]))
                hl0x.Fill(float(Emax),float(totPEl0x[i]))
                hl1y.Fill(float(Emax),float(totPEl1y[i]))
                hl1x.Fill(float(Emax),float(totPEl1x[i]))
                hl2y.Fill(float(Emax),float(totPEl2y[i]))
                hl2x.Fill(float(Emax),float(totPEl2x[i]))        
                
                hits = [[[],[]],[[],[]],[[],[]]]
                PE = [[[],[]],[[],[]],[[],[]]]
                for ilay in range(3):
                    for iview in range(2):
                        hits[ilay][iview]=(np.where(PEstrip[ilay][iview] != 0)[0])
                        PE[ilay][iview]=(PEstrip[ilay][iview][PEstrip[ilay][iview] != 0])
                        #print ("WWWWWWWWWWWWWWWWWWWWWWWWWW",hits[ilay][iview])
                        nclu,istrip0,nstrip0,charge0 = Clustering(hits[ilay][iview],PE[ilay][iview])  
                        Nclu += nclu
                       
                    
                        for kk in range(len(istrip0)):
                            Fstrip.append(istrip0[kk])
                            Clusiz.append(nstrip0[kk])
                            Clunpe.append(charge0[kk])
                            
                        print ("Layer: ",ilay, "View: ", iview)
                        print("fiber hit: ",hits[ilay][iview])
                        print ("nclu",nclu, "firs strip:", istrip0, "cluster siza: ", nstrip0,"npe",charge0 )
                        print ("nclu",Nclu, "firs strip:", Fstrip, "cluster siza: ",Clusiz,"npe", Clunpe )
                       
                    ncluH,istrip0H,nstrip0H,charge0H = Clustering(hits[ilay][0],PE[ilay][0])  
                    NcluH += ncluH 
                    
                    for kk in range(len(istrip0H)):
                        FstripH.append(istrip0H[kk])
                        ClusizH.append(nstrip0H[kk])
                        ClunpeH.append(charge0H[kk])
                    
                    ncluV,istrip0V,nstrip0V,charge0V = Clustering(hits[ilay][1],PE[ilay][1])  
                    NcluV += ncluV 
                         
                    for kk in range(len(istrip0V)):
                        FstripV.append(istrip0V[kk])
                        ClusizV.append(nstrip0V[kk])
                        ClunpeV.append(charge0V[kk])                      
           
            
                Tclu[0] = Nclu
                TcluHor[0] = NcluH
                TcluVer[0] = NcluV
                for ihit in range (Tch[0]):
                    ILay[ihit] = Lay[ihit]
                    IView[ihit] = View[ihit]
                    Ipix[ihit] = Fiber[ihit]
                    Npe[ihit] = npe_vec[ihit]
                for ihit in range (Tchs[0]):
                    print(Tchs[0])
                    Istrip[ihit] = Strip_vec[ihit]
                    NpeStrip[ihit] = PEStrip_vec[ihit]
                for icl in range(Nclu):
                    Csiz[icl] = Clusiz[icl]
                    Tnpe[icl] = Clunpe[icl]
                    Fstr[icl] = Fstrip[icl]
                for icl in range(NcluH):
                    CsizHor[icl] = ClusizH[icl]
                    TnpeHor[icl] = ClunpeH[icl]
                    FstrHor[icl] = FstripH[icl]
                for icl in range(NcluV):
                    CsizVer[icl] = ClusizV[icl]
                    TnpeVer[icl] = ClunpeV[icl]
                    FstrVer[icl] = FstripV[icl]
                '''print ("--------------------------")
                print ("Event:", EventID)
                print ("Tch:", Tch)
                print ("Lay",Lay,"View",View)
                print ("Fiber",Fiber)
                print ("Strip", Strip_vec)
                print ("NPE", Npe)
                print ("PEStrip",PEStrip_vec)'''
                #print ("nclu",Nclu, "firs strip:", Fstr[:Nclu], "cluster siza: ",Csiz[:Nclu],"npe", Tnpe[:Nclu])
                #print ("@@@@@@@@@@@@@@@@@@@@@",np.array(layer_vec[i]))
                OutTree.Fill()
            #fine ciclo sugli eventi
            '''                                     
            print ("PE layer 0 view y ",sum(totPEl0y))
            print ("PE layer 0 view x ",sum(totPEl0x))
            print ("PE layer 1 view y ",sum(totPEl1y))
            print ("PE layer 1 view x ",sum(totPEl1x))
            print ("PE layer 2 view y ",sum(totPEl2y))
            print ("PE layer 2 view x ",sum(totPEl2x))
            '''

            tfile.Close()
            
        #fine ciclo sui file    
    
    
    
    canvas[0].cd(2)
    hl0y.Draw("colz")
    hl0y.GetXaxis().SetTitle("Energy Kev")
    hl0y.GetYaxis().SetTitle("Total Photo Electrons")
    #hl0y.GetYaxis().SetRangeUser(0,maxy)
    canvas[0].cd(1)
    hl0x.Draw("colz")
    hl0x.GetXaxis().SetTitle("Energy Kev")
    hl0x.GetYaxis().SetTitle("Total Photo Electrons")
    #hl0x.GetYaxis().SetRangeUser(0,maxy)
    fout.cd()
    canvas[0].Write()

    canvas[1].cd(2)
    hl1y.Draw("colz")
    hl1y.GetXaxis().SetTitle("Energy Kev")
    hl1y.GetYaxis().SetTitle("Total Photo Electrons")
    #hl0y.GetYaxis().SetRangeUser(0,maxy)
    canvas[1].cd(1)
    hl1x.Draw("colz")
    hl1x.GetXaxis().SetTitle("Energy Kev")
    hl1x.GetYaxis().SetTitle("Total Photo Electrons")
    #hl0x.GetYaxis().SetRangeUser(0,maxy)
    fout.cd()
    canvas[1].Write()

    canvas[2].cd(2)
    hl2y.Draw("colz")
    hl2y.GetXaxis().SetTitle("Energy Kev")
    hl2y.GetYaxis().SetTitle("Total Photo Electrons")
    #hl0y.GetYaxis().SetRangeUser(0,maxy)
    canvas[2].cd(1)
    hl2x.Draw("colz")
    hl2x.GetXaxis().SetTitle("Energy Kev")
    hl2x.GetYaxis().SetTitle("Total Photo Electrons")
    #hl0x.GetYaxis().SetRangeUser(0,maxy)
    fout.cd()
    canvas[2].Write()
        
    OutTree.Write()
    fout.Close()
    #print "OK"
