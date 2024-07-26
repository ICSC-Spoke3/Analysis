import os
import sys
import time
import ROOT
from ROOT import TFile, TTree, TChain, TH1D, TH2D, TMath, TCanvas, TH1F, TRandom3, TObject, TLegend, TColor, gStyle, gPad
import array as ary
import numpy as np
import glob
ran = TRandom3()

def Clustering(hits, charges):
    nhits = len(hits)

    hits = np.array(hits)
    index = hits.argsort()
    
    maxclus = len(hits)
    istrip0 = np.zeros(maxclus,dtype = int)
    nstrip0 = np.zeros(maxclus,dtype = int)
    charge0 = np.zeros(maxclus)

    hits1 = []
    charge1 = []
    hits1.append(hits[index[0]])
    charge1.append(charges[index[0]])
    for i in range (1,nhits):
        if hits[index[i]]==hits1[-1]:
            charge1[-1]+= charges[index[i]]
        else:
            hits1.append(hits[index[i]])
            charge1.append(charges[index[i]])
    istrip = 0
    nclu = 0
    istrip0[nclu] = hits1[istrip]
    charge0[nclu] = charge1[istrip]
    nstrip0[nclu] = 1
    
    nhits1=len(hits1)
    nclu += 1
    istrip += 1
    while (istrip < nhits1):
        if ((hits1[istrip]==hits1[istrip-1]+1)):
            nstrip0[nclu-1] = nstrip0[nclu-1] + 1
            charge0[nclu-1] += charge1[istrip]
        else:
            istrip0[nclu] = hits1[istrip]
            nstrip0[nclu] = 1
            charge0[nclu] = charge1[istrip]
            nclu += 1        
        istrip += 1
   

    istrip0 = istrip0[:nclu]
    nstrip0 = nstrip0[:nclu]
    charge0 = charge0[:nclu]

    return nclu,istrip0,nstrip0,charge0

def encoder(layer,view,fiber):
    pixel = ary.array('i',range(len(layer)))
    for j in range(len(layer)):
        if layer[j] == 0 and view[j] == 0:
            pixel[j]=fiber[j]
        if layer[j] == 0 and view[j] == 1:
            pixel[j]=fiber[j]+nfiber
        if layer[j] == 1 and view[j] == 0:
            pixel[j]=fiber[j]+2*nfiber
        if layer[j] == 1 and view[j] == 1:
            pixel[j]=fiber[j]+3*nfiber
        if layer[j] == 2 and view[j] == 0:
            pixel[j]=fiber[j]+4*nfiber
        if layer[j] == 2 and view[j] == 1:
            pixel[j]=fiber[j]+5*nfiber
    return pixel


if __name__ == '__main__':

    fin = sys.argv[1] 
    fout = TFile("NusesResults.root","recreate")
    
    nview = 2
    totLayers = 3
    nfiber = 130
    fiberyield = 8*0.92#kev/pe
    PDE = 0.4
    trapeff = 0.054
    PEthresh = 3 #pe per strip
        
    nChan = 100000 
    eventID = ary.array('i',[0])
    Tch     = ary.array('i',[0])
    energy  = ary.array('f',[0])
    Ichip   = ary.array('i',range(nChan))
    Ichan   = ary.array('i',range(nChan))
    Ipix    = ary.array('i',range(nChan))
    Npe     = ary.array('f',range(nChan))
    
    nClusters = 1000 
    Tclu = ary.array('i',[0])
    dim = ary.array('i',[0])
    Tnpe = ary.array('f',range(nClusters))
    Csiz = ary.array('i',range(nClusters))
    Fstr = ary.array('i',range(nClusters))
    Mtim = ary.array('f',range(nClusters))

    TstripHor0   = ary.array('i',[0])
    TpeStripHor0 = ary.array('f',[0])
    stripHor0    = ary.array('i',range(nChan))
    PEstripHor0  = ary.array('f',range(nChan))
   
    dimHor0 = ary.array('i',[0])
    TcluHor0 = ary.array('i',[0])
    CsizHor0 = ary.array('i',range(nClusters))
    TnpeHor0 = ary.array('f',range(nClusters))
    FstrHor0 = ary.array('i',range(nClusters))
    
    TstripVer0   = ary.array('i',[0])
    TpeStripVer0 = ary.array('f',[0])
    stripVer0    = ary.array('f',range(nChan))
    PEstripVer0  = ary.array('f',range(nChan))
    
    dimVer0 = ary.array('i',[0])
    TcluVer0 = ary.array('i',[0])
    CsizVer0 = ary.array('i',range(nClusters))
    TnpeVer0 = ary.array('f',range(nClusters))
    FstrVer0 = ary.array('i',range(nClusters))
    
    TstripHor1   = ary.array('i',[0])
    TpeStripHor1 = ary.array('f',[0])
    stripHor1    = ary.array('f',range(nChan))
    PEstripHor1  = ary.array('f',range(nChan))
    
    dimHor1 = ary.array('i',[0])
    TcluHor1 = ary.array('i',[0])
    CsizHor1 = ary.array('i',range(nClusters))
    TnpeHor1 = ary.array('f',range(nClusters))
    FstrHor1 = ary.array('i',range(nClusters))

    TstripVer1   = ary.array('i',[0])
    TpeStripVer1 = ary.array('f',[0])
    stripVer1    = ary.array('i',range(nChan))
    PEstripVer1  = ary.array('f',range(nChan))
    
    dimVer1 = ary.array('i',[0])
    TcluVer1 = ary.array('i',[0])
    CsizVer1 = ary.array('i',range(nClusters))
    TnpeVer1 = ary.array('f',range(nClusters))
    FstrVer1 = ary.array('i',range(nClusters))
    
    TstripHor2   = ary.array('i',[0])
    TpeStripHor2 = ary.array('f',[0])
    stripHor2    = ary.array('i',range(nChan))
    PEstripHor2  = ary.array('f',range(nChan))
    
    dimHor2 = ary.array('i',[0])
    TcluHor2 = ary.array('i',[0])
    CsizHor2 = ary.array('i',range(nClusters))
    TnpeHor2 = ary.array('f',range(nClusters))
    FstrHor2 = ary.array('i',range(nClusters))
    MtimHor2 = ary.array('f',range(nClusters))

    TstripVer2   = ary.array('i',[0])
    TpeStripVer2 = ary.array('f',[0])
    stripVer2    = ary.array('i',range(nChan))
    PEstripVer2  = ary.array('f',range(nChan))
    
    dimVer2 = ary.array('i',[0])
    TcluVer2 = ary.array('i',[0])
    CsizVer2 = ary.array('i',range(nClusters))
    TnpeVer2 = ary.array('f',range(nClusters))
    FstrVer2 = ary.array('i',range(nClusters))
    
    OutTree    = TTree("ClusterTree", "ClusterTree")
    
    # Define output tree
    OutTree.Branch("eventID", eventID, "eventID/I")
    OutTree.Branch("Tch", Tch, "Tch/I")
    OutTree.Branch("energy", energy, "energy_keV/F")
    OutTree.Branch("Ichip", Ichip, "Ichipx[Tch]/I")
    OutTree.Branch("Ichan", Ichan, "Ichan[Tch]/I")
    OutTree.Branch("Ipix", Ipix, "Ipix[Tch]/I")
    OutTree.Branch("Npe", Npe, "Npe[Tch]/F")
   
    # Informazioni generali sui cluster
    OutTree.Branch("Tclu", Tclu, "Tclu/I")
    OutTree.Branch("dim", dim, "dim/I")
    OutTree.Branch("Csiz", Csiz, "Csiz[dim]/I")
    OutTree.Branch("Tnpe", Tnpe, "Tnpe[dim]/F")
    OutTree.Branch("Fstr", Fstr, "Fstr[dim]/I")
    
    # Informazioni sui cluster nelle viste orizzontale(x=1) e verticale(y=0) per i 3 layer
    OutTree.Branch("TstripHor0", TstripHor0, "TstripHor0/I")
    OutTree.Branch("TpeStripHor0", TpeStripHor0, "TpeStripHor0/F")
    OutTree.Branch("stripHor0", stripHor0, "stripHor0[TstripHor0]/I")
    OutTree.Branch("PEstripHor0", PEstripHor0, "PEstripHor0[TstripHor0]/F")
    
    OutTree.Branch("dimHor0", dimHor0, "dimHor0/I")
    OutTree.Branch("TcluHor0", TcluHor0, "TcluHor0/I")
    OutTree.Branch("CsizHor0", CsizHor0, "CsizHor0[dimHor0]/I")
    OutTree.Branch("TnpeHor0", TnpeHor0, "TnpeHor0[dimHor0]/F")
    OutTree.Branch("FstrHor0", FstrHor0, "FstrHor0[dimHor0]/I")

    OutTree.Branch("TstripVer0", TstripVer0, "TstripVer0/I")
    OutTree.Branch("TpeStripVer0", TpeStripVer0, "TpeStripVer0/F")
    OutTree.Branch("stripVer0", stripVer0, "stripVer0[TstripVer0]/I")
    OutTree.Branch("PEstripVer0", PEstripVer0, "PEstripVer0[TstripVer0]/F")
    
    OutTree.Branch("dimVer0", dimVer0, "dimVer0/I")
    OutTree.Branch("TcluVer0", TcluVer0, "TcluVer0/I")
    OutTree.Branch("CsizVer0", CsizVer0, "CsizVer0[dimVer0]/I")
    OutTree.Branch("TnpeVer0", TnpeVer0, "TnpeVer0[dimVer0]/F")
    OutTree.Branch("FstrVer0", FstrVer0, "FstrVer0[dimVer0]/I")
    
    OutTree.Branch("TstripHor1", TstripHor1, "TstripHor1/I")
    OutTree.Branch("TpeStripHor1", TpeStripHor1, "TpeStripHor1/F")
    OutTree.Branch("stripHor1", stripHor1, "stripHor1[TstripHor1]/I")
    OutTree.Branch("PEstripHor1", PEstripHor1, "PEstripHor1[TstripHor1]/F")

    OutTree.Branch("dimHor1", dimHor1, "dimHor1/I")
    OutTree.Branch("TcluHor1", TcluHor1, "TcluHor1/I")
    OutTree.Branch("CsizHor1", CsizHor1, "CsizHor1[dimHor1]/I")
    OutTree.Branch("TnpeHor1", TnpeHor1, "TnpeHor1[dimHor1]/F")
    OutTree.Branch("FstrHor1", FstrHor1, "FstrHor1[dimHor1]/I")

    OutTree.Branch("TstripVer1", TstripVer1, "TstripVer1/I")
    OutTree.Branch("TpeStripVer1", TpeStripVer1, "TpeStripVer1/F")
    OutTree.Branch("stripVer1", stripVer1, "stripVer1[TstripVer1]/I")
    OutTree.Branch("PEstripVer1", PEstripVer1, "PEstripVer1[TstripVer1]/F")
    
    OutTree.Branch("dimVer1", dimVer1, "dimVer1/I")
    OutTree.Branch("TcluVer1", TcluVer1, "TcluVer1/I")
    OutTree.Branch("CsizVer1", CsizVer1, "CsizVer1[dimVer1]/I")
    OutTree.Branch("TnpeVer1", TnpeVer1, "TnpeVer1[dimVer1]/F")
    OutTree.Branch("FstrVer1", FstrVer1, "FstrVer1[dimVer1]/I")
    
    OutTree.Branch("TstripHor2", TstripHor2, "TstripHor2/I")
    OutTree.Branch("TpeStripHor2", TpeStripHor2, "TpeStripHor2/F")
    OutTree.Branch("stripHor2", stripHor2, "stripHor2[TstripHor2]/I")
    OutTree.Branch("PEstripHor2", PEstripHor2, "PEstripHor2[TstripHor2]/F")
    
    OutTree.Branch("Hor2", TcluHor2, "Hor2/I")
    OutTree.Branch("TcluHor2", TcluHor2, "TcluHor2/I")
    OutTree.Branch("CsizHor2", CsizHor2, "CsizHor2[Hor2]/I")
    OutTree.Branch("TnpeHor2", TnpeHor2, "TnpeHor2[Hor2]/F")
    OutTree.Branch("FstrHor2", FstrHor2, "FstrHor2[Hor2]/I")

    OutTree.Branch("TstripVer2", TstripVer2, "TstripVer2/I")
    OutTree.Branch("TpeStripVer2", TpeStripVer2, "TpeStripVer2/F")
    OutTree.Branch("stripVer2", stripVer2, "stripVer2[TstripVer2]/I")
    OutTree.Branch("PEstripVer2", PEstripVer2, "PEstripVer2[TstripVer2]/F")
    
    OutTree.Branch("dimVer2", dimVer2, "dimVer2/I")
    OutTree.Branch("TcluVer2", TcluVer2, "TcluVer2/I")
    OutTree.Branch("CsizVer2", CsizVer2, "CsizVer2[dimVer2]/I")
    OutTree.Branch("TnpeVer2", TnpeVer2, "TnpeVer2[dimVer2]/F")
    OutTree.Branch("FstrVer2", FstrVer2, "FstrVer2[dimVer2]/I")
    
    #definizione istogrammi
    canvas = []
    for i in range(totLayers):
        canvas.append([])
        canvas[i]=TCanvas("c{0}".format(i),"c{0}".format(i))
        canvas[i].Divide(2,2)  
    cef = TCanvas("cef","cef")
    cnum = TCanvas("cnum","cnum")
    cden = TCanvas("cden","cden") 
    
    hl0y = TH2D("hl0y","layer0 y-view",991,90,10000,200,-0.5,100) 
    #hl0y.SetContour(10)
    #hl0y.GetZaxis().SetRangeUser(0,80)

    hl0x = TH2D("hl0x","layer0 x-view",991,90,10000,200,-0.5,100) 
    #hl0x.SetContour(20)
    #hl0x.GetZaxis().SetRangeUser(0,80)
    hl0x.GetXaxis().SetTitle("Energy keV")
    hl0x.GetYaxis().SetTitle("Total Photo Electrons")
    
    hl1y = TH2D("hl1y","layer1 y-view",991,90,10000,200,-0.5,100) 
    #hl1y.SetContour(30)
    #hl1y.GetZaxis().SetRangeUser(0,80)
    hl1y.GetXaxis().SetTitle("Energy keV")
    hl1y.GetYaxis().SetTitle("Total Photo Electrons")
    
    hl1x = TH2D("hl1x","layer1 x-view",991,90,10000,200,-0.5,100) 
    #hl1x.SetContour(40)
    #hl1x.GetZaxis().SetRangeUser(0,80)
    hl1x.GetXaxis().SetTitle("Energy keV")
    hl1x.GetYaxis().SetTitle("Total Photo Electrons")

    hl2y = TH2D("hl2y","layer2 y-view",991,90,10000,200,-0.5,100) 
    #hl2y.SetContour(50)
    #hl2y.GetZaxis().SetRangeUser(0,80)
    hl2y.GetXaxis().SetTitle("Energy keV")
    hl2y.GetYaxis().SetTitle("Total Photo Electrons")
    
    hl2x = TH2D("hl2x","layer2 x-view",991,90,10000,200,-0.5,100) 
    #hl2x.SetContour(60)
    #hl2x.GetZaxis().SetRangeUser(0,80)
    hl2x.GetXaxis().SetTitle("Energy keV")
    hl2x.GetYaxis().SetTitle("Total Photo Electrons")

    hl0yclu = TH2D("hl0yclu","layer0 y-view",991,90,10000,200,-0.5,100) 
    #hl0yclu.SetContour(10)
    #hl0yclu.GetZaxis().SetRangeUser(0,80)

    hl0xclu = TH2D("hl0xclu","layer0 x-view",991,90,10000,200,-0.5,100) 
    #hl0xclu.SetContour(20)
    #hl0xclu.GetZaxis().SetRangeUser(0,80)

    hl1yclu = TH2D("hl1yclu","layer1 y-view",991,90,10000,200,-0.5,100) 
    #hl1yclu.SetContour(30)
    #hl1yclu.GetZaxis().SetRangeUser(0,80)

    hl1xclu = TH2D("hl1xclu","layer1 x-view",991,90,10000,200,-0.5,100) 
    #hl1xclu.SetContour(40)
    #hl1xclu.GetZaxis().SetRangeUser(0,80)

    hl2yclu = TH2D("hl2yclu","layer2 y-view",991,90,10000,200,-0.5,100) 
    #hl2yclu.SetContour(50)
    #hl2yclu.GetZaxis().SetRangeUser(0,80)
    
    hl2xclu = TH2D("hl2xclu","layer2 x-view",991,90,10000,200,-0.5,100) 
    #hl2xclu.SetContour(60)
    #hl2xclu.GetZaxis().SetRangeUser(0,80)

    ekmin = 0
    ekmax = 10001
    nek = 50
    hEkDist = TH1D("hEkDist", "Distribution of kinetic energies (keV)", nek, ekmin, ekmax)
    hNTrigViewVsEk = []
    hEffViewVsEk = []
    for i in range (nview*totLayers):
        title1 = "hNTrigViewVsEk[" + str(i) + "]"
        title2 = "Number of triggers vs kinetic energy - " + str(i)
        hNTrigViewVsEk.append(TH1D(title1, title2, nek, ekmin, ekmax))
        title1 = "hEffViewVsEk[" + str(i) + "]"
        title2 = "Trigger efficiency vs kinetic energy - " + str(i)
        hEffViewVsEk.append(TH1D(title1, title2, nek, ekmin, ekmax))
    
    '''dE_vec = []
    pe_vec = []
    layer_vec = []
    view_vec = []
    fiber_vec = []
    pixel_vec = []
    pe_vec = []'''
    #PrimaryEnergy = []

    rootfile = fin
    
    #print ("[INFO] ROOT file ",rootfile)
    
    tfile = TFile(rootfile)
    tree = tfile.Get("TrackerDigi")
    nhit = tree.GetEntries()
    nhit = 100
    
    
    
    for i in range(nhit):
        pe_vec = []
        Strip_vec = []
        Strip_vec0y = []
        Strip_vec0x = []
        Strip_vec1y = []
        Strip_vec1x = []
        Strip_vec2y = []
        Strip_vec2x = []
        PEstrip_vec = []
        PEstrip_vec0y = []
        PEstrip_vec0x = []
        PEstrip_vec1y = []
        PEstrip_vec1x = []
        PEstrip_vec2y = []
        PEstrip_vec2x = []
                
        totPEl0y = 0
        totPEl0x = 0
        totPEl1y = 0
        totPEl1x = 0
        totPEl2y = 0
        totPEl2x = 0
            
        Lay = []
        View = []
        
        Nclu = 0
        Fstrip = []
        Clusiz = []
        Clunpe = []
        
        NcluH0 = 0
        FstripH0 = []
        ClusizH0 = []
        ClunpeH0 = []
        
        NcluV0 = 0
        FstripV0 = []
        ClusizV0 = []
        ClunpeV0 = []
        
        NcluH1 = 0
        FstripH1 = []
        ClusizH1 = []
        ClunpeH1 = []
        
        NcluV1 = 0
        FstripV1 = []
        ClusizV1 = []
        ClunpeV1 = []
        
        NcluH2 = 0
        FstripH2 = []
        ClusizH2 = []
        ClunpeH2 = []
        
        NcluV2 = 0
        FstripV2 = []
        ClusizV2 = []
        ClunpeV2 = []
        
        eventID[0] = i
        
        #pe_vec.append([])
        
        tree.GetEntry(i)  
        
        #PrimaryEnergy.append(tree.PrimaryParticleEnergy)
        PrimaryEnergy = (tree.PrimaryParticleEnergy)*1e3
        dE_vec        = tree.Energy_keV
        view_vec      = tree.ViewNo
        layer_vec     = tree.LayerNo
        fiber_vec     = tree.FiberNo
        pixel_vec     = encoder(layer_vec,view_vec,fiber_vec)

        '''dE_vec.append(np.array(tree.Energy_keV))
        layer_vec.append(tree.LayerNo)
        view_vec.append(tree.ViewNo)
        fiber_vec.append(tree.FiberNo)
        pixel_vec.append(encoder(layer_vec[i],view_vec[i],fiber_vec[i]))'''

        energy[0] = PrimaryEnergy
        hEkDist.Fill(PrimaryEnergy, 1.0)
        
        
        for j in range(len(dE_vec)):
            #pe_vec[i].append(ran.Poisson(dE_vec[i][j]*fiberyield*trapeff*PDE))
            pe_vec.append(dE_vec[j]*fiberyield*trapeff*PDE)
            Lay.append(layer_vec[j])
            View.append(view_vec[j])
            
            if pixel_vec[j] == 0:
                index1 = 1
                index2 = 1
            elif pixel_vec[j] == 129:
                index1 = nfiber-1
                index2 = nfiber-1
            else:
                index1 = pixel_vec[j]
                index2 = pixel_vec[j]-1
                    
            Strip_vec.append(index1)
            Strip_vec.append(index2)
            PEstrip_vec.append(ran.Poisson(pe_vec[j]/2))
            PEstrip_vec.append(ran.Poisson(pe_vec[j]/2))
            
            if (layer_vec[j] == 0 and view_vec[j] == 0):
                totPEl0y += ran.Poisson(pe_vec[j]) 
                Strip_vec0y.append(index1)
                Strip_vec0y.append(index2)
                PEstrip_vec0y.append(ran.Poisson(pe_vec[j]/2))
                PEstrip_vec0y.append(ran.Poisson(pe_vec[j]/2))
            if (layer_vec[j] == 0 and view_vec[j] == 1):
                totPEl0x += ran.Poisson(pe_vec[j])
                Strip_vec0x.append(index1)
                Strip_vec0x.append(index2)
                PEstrip_vec0x.append(ran.Poisson(pe_vec[j]/2))
                PEstrip_vec0x.append(ran.Poisson(pe_vec[j]/2))
            if (layer_vec[j] == 1 and view_vec[j] == 0):
                totPEl1y += ran.Poisson(pe_vec[j])
                Strip_vec1y.append(index1)
                Strip_vec1y.append(index2)
                PEstrip_vec1y.append(ran.Poisson(pe_vec[j]/2))
                PEstrip_vec1y.append(ran.Poisson(pe_vec[j]/2))
            if (layer_vec[j] ==1 and view_vec[j] == 1):
                totPEl1x += ran.Poisson(pe_vec[j])
                Strip_vec1x.append(index1)
                Strip_vec1x.append(index2)
                PEstrip_vec1x.append(ran.Poisson(pe_vec[j]/2))
                PEstrip_vec1x.append(ran.Poisson(pe_vec[j]/2))
            if (layer_vec[j] == 2 and view_vec[j] == 0):
                totPEl2y += ran.Poisson(pe_vec[j])
                Strip_vec2y.append(index1)
                Strip_vec2y.append(index2)
                PEstrip_vec2y.append(ran.Poisson(pe_vec[j]/2))
                PEstrip_vec2y.append(ran.Poisson(pe_vec[j]/2))
            if (layer_vec[j] ==2 and view_vec[j] == 1):
                totPEl2x += ran.Poisson(pe_vec[j])
                Strip_vec2x.append(index1)
                Strip_vec2x.append(index2)
                PEstrip_vec2x.append(ran.Poisson(pe_vec[j]/2))
                PEstrip_vec2x.append(ran.Poisson(pe_vec[j]/2))
           
            
        #####fine ciclo sulle hit dell'evento i     
        
        TpeStripHor0[0] = totPEl0y
        TpeStripVer0[0] = totPEl0x
        TpeStripHor1[0] = totPEl1y
        TpeStripVer1[0] = totPEl1x
        TpeStripHor2[0] = totPEl2y
        TpeStripVer2[0] = totPEl2x
        totPE = totPEl0y + totPEl0x + totPEl1y + totPEl1x + totPEl2y + totPEl2x

        vtrig= np.zeros(nview*totLayers)

        for istrip in range(len(Strip_vec0y)-1):
            if (PEstrip_vec0y[istrip] >= PEthresh and PEstrip_vec0y[istrip+1] >= PEthresh):
                vtrig[0] = 1
                
        for istrip in range(len(Strip_vec0x)-1):
            if (PEstrip_vec0x[istrip] >= PEthresh and PEstrip_vec0x[istrip+1] >= PEthresh):
                vtrig[1] = 1
                
        for istrip in range(len(Strip_vec1y)-1):
            if (PEstrip_vec1y[istrip] >= PEthresh and PEstrip_vec1y[istrip+1] >= PEthresh):
                vtrig[2] = 1
                
        for istrip in range(len(Strip_vec1x)-1):
            if (PEstrip_vec1x[istrip] >= PEthresh and PEstrip_vec1x[istrip+1] >= PEthresh):
                vtrig[3] = 1
                    
        for istrip in range(len(Strip_vec2y)-1):
            if (PEstrip_vec2y[istrip] >= PEthresh and PEstrip_vec2y[istrip+1] >= PEthresh):
                vtrig[4] = 1

        for istrip in range(len(Strip_vec2x)-1):
            if (PEstrip_vec2x[istrip] >= PEthresh and PEstrip_vec2x[istrip+1] >= PEthresh):
                vtrig[5] = 1

        for ill in range(nview*totLayers):
            if vtrig[ill]>0:
                hNTrigViewVsEk[ill].Fill(PrimaryEnergy, 1.0)
                    
                    
        Tch[0] = (len(PEstrip_vec))
        if len(PEstrip_vec)>0:
            nclu,istrip0,nstrip0,charge0 = Clustering(Strip_vec,PEstrip_vec)  
            Nclu += nclu
                    
            for kk in range(len(istrip0)):
                Fstrip.append(istrip0[kk])
                Clusiz.append(nstrip0[kk])
                Clunpe.append(charge0[kk])
                
        TstripHor0[0] = (len(PEstrip_vec0y))
        if len(PEstrip_vec0y)>0:
            ncluH0,istrip0H0,nstrip0H0,charge0H0 = Clustering(Strip_vec0y,PEstrip_vec0y)
            NcluH0 += ncluH0
            
            for kk in range(len(istrip0H0)):
                FstripH0.append(istrip0H0[kk])
                ClusizH0.append(nstrip0H0[kk])
                ClunpeH0.append(charge0H0[kk])
                
        TstripVer0[0] = (len(PEstrip_vec0x))
        if len(PEstrip_vec0x)>0:    
            ncluV0,istrip0V0,nstrip0V0,charge0V0 = Clustering(Strip_vec0x,PEstrip_vec0x)  
            NcluV0 += ncluV0

            for kk in range(len(istrip0V0)):
                FstripV0.append(istrip0V0[kk])
                ClusizV0.append(nstrip0V0[kk])
                ClunpeV0.append(charge0V0[kk])
                
        TstripHor1[0] = (len(PEstrip_vec1y))
        if len(PEstrip_vec1y)>0:     
            ncluH1,istrip0H1,nstrip0H1,charge0H1 = Clustering(Strip_vec1y,PEstrip_vec1y)  
            NcluH1 += ncluH1
        
            for kk in range(len(istrip0H1)):
                FstripH1.append(istrip0H1[kk])
                ClusizH1.append(nstrip0H1[kk])
                ClunpeH1.append(charge0H1[kk])
                
        TstripVer1[0] = (len(PEstrip_vec1x))
        if len(PEstrip_vec1x)>0:     
            ncluV1,istrip0V1,nstrip0V1,charge0V1 = Clustering(Strip_vec1x,PEstrip_vec1x)  
            NcluV1 += ncluV1
            
            for kk in range(len(istrip0V1)):
                FstripV1.append(istrip0V1[kk])
                ClusizV1.append(nstrip0V1[kk])
                ClunpeV1.append(charge0V1[kk])
                
        TstripHor2[0] = (len(PEstrip_vec2y))
        if len(PEstrip_vec2y)>0:     
            ncluH2,istrip0H2,nstrip0H2,charge0H2 = Clustering(Strip_vec2y,PEstrip_vec2y)  
            NcluH2 += ncluH2

            for kk in range(len(istrip0H2)):
                FstripH2.append(istrip0H2[kk])
                ClusizH2.append(nstrip0H2[kk])
                ClunpeH2.append(charge0H2[kk])
                
        TstripVer2[0] = (len(PEstrip_vec2x))
        if len(PEstrip_vec2x)>0:     
            ncluV2,istrip0V2,nstrip0V2,charge0V2 = Clustering(Strip_vec2x,PEstrip_vec2x)  
            NcluV2 += ncluV2
            
            for kk in range(len(istrip0V2)):
                FstripV2.append(istrip0V2[kk])
                ClusizV2.append(nstrip0V2[kk])
                ClunpeV2.append(charge0V2[kk])
                
        Tclu[0] = Nclu
        TcluHor0[0] = NcluH0
        TcluVer0[0] = NcluV0
        TcluHor1[0] = NcluH1
        TcluVer1[0] = NcluV1
        TcluHor2[0] = NcluH2
        TcluVer2[0] = NcluV2
        dim[0] = Nclu
        dimHor0[0] = NcluH0
        dimVer0[0] = NcluV0
        dimHor1[0] = NcluH1
        dimVer1[0] = NcluV1
        dimHor2[0] = NcluH2
        dimVer2[0] = NcluV2

        
        #print ("ENERGY:", energy[0])
        for ihit in range (Tch[0]):
            #ILay[ihit] = Lay[ihit]
            #IView[ihit] = View[ihit]
            #energy[ihit] = Emax
            Ipix[ihit] = Strip_vec[ihit]
            Npe[ihit] = PEstrip_vec[ihit]
        #print("strip: ",np.array(Ipix[:Tch[0]]))
    
                                        
        for icl in range(Nclu):
            Csiz[icl] = Clusiz[icl]
            Tnpe[icl] = Clunpe[icl]
            Fstr[icl] = Fstrip[icl]
        if Nclu==0: 
            Csiz[0] = 0
            Tnpe[0] = 0
            dim[0] = 1
            '''print ("Nclu =0")
            print ("Csiz: ", Csiz[0], "Tnpe: ", Tnpe[0])'''
                    
                
        for istrip in range (TstripHor0[0]):
            stripHor0[istrip]   = Strip_vec0y[istrip]
            PEstripHor0[istrip] = PEstrip_vec0y[istrip]
        for icl in range(NcluH0):
            CsizHor0[icl] = ClusizH0[icl]
            TnpeHor0[icl] = ClunpeH0[icl]
            FstrHor0[icl] = FstripH0[icl]
        if NcluH0==0: 
            CsizHor0[0] = 0
            TnpeHor0[0] = 0
            dimHor0[0] = 1
            '''print ("NcluH0 =0")
            print ("CsizH0: ", CsizHor0[0], "Tnpe: ", TnpeHor0[0])'''
                
        for istrip in range (TstripVer0[0]):
            stripVer0[istrip]   = Strip_vec0x[istrip]
            PEstripVer0[istrip] = PEstrip_vec0x[istrip]
        for icl in range(NcluV0):
            CsizVer0[icl] = ClusizV0[icl]
            TnpeVer0[icl] = ClunpeV0[icl]
            FstrVer0[icl] = FstripV0[icl]
        if NcluV0==0: 
            CsizVer0[0] = 0
            TnpeVer0[0] = 0
            dimVer0[0] = 1
            '''print ("NcluV0 =0")
            print ("CsizV0: ", CsizVer0[0], "Tnpe: ", TnpeVer0[0])'''
                    
                    
        for istrip in range (TstripHor1[0]):
            stripHor1[istrip]   = Strip_vec1y[istrip]
            PEstripHor1[istrip] = PEstrip_vec1y[istrip]
        for icl in range(NcluH1):
            CsizHor1[icl] = ClusizH1[icl]
            TnpeHor1[icl] = ClunpeH1[icl]
            FstrHor1[icl] = FstripH1[icl]
        if NcluH1==0: 
            CsizHor1[0] = 0
            TnpeHor1[0] = 0
            dimHor1[0] = 1
            '''print ("NcluH1 =0")
            print ("CsizH1: ", CsizHor1[0], "Tnpe: ", TnpeHor1[0])'''
            
        for istrip in range (TstripVer1[0]):
            stripVer1[istrip]   = Strip_vec1x[istrip]
            PEstripVer1[istrip] = PEstrip_vec1x[istrip]
        for icl in range(NcluV1):
            CsizVer1[icl] = ClusizV1[icl]
            TnpeVer1[icl] = ClunpeV1[icl]
            FstrVer1[icl] = FstripV1[icl]
        if NcluV1==0: 
            CsizVer1[0] = 0
            TnpeVer1[0] = 0
            dimVer1[0] = 1
            '''print ("NcluV1 =0")
            print ("CsizV1: ", CsizVer1[0], "Tnpe: ", TnpeVer1[0])'''
                
        for istrip in range (TstripHor2[0]):
            stripHor2[istrip]   = Strip_vec2y[istrip]
            PEstripHor2[istrip] = PEstrip_vec2y[istrip]                    
        for icl in range(NcluH2):
            CsizHor2[icl] = ClusizH2[icl]
            TnpeHor2[icl] = ClunpeH2[icl]
            FstrHor2[icl] = FstripH2[icl]
        if NcluH2==0: 
            CsizHor2[0] = 0
            TnpeHor2[0] = 0
            dimHor2[0] = 1
            '''print ("NcluH2 =0")
            print ("CsizH2: ", CsizHor2[0], "Tnpe: ", TnpeHor2[0])'''   
                    
        for istrip in range (TstripVer2[0]):
            stripVer2[istrip]   = Strip_vec2x[istrip]
            PEstripVer2[istrip] = PEstrip_vec2x[istrip]                    
        for icl in range(NcluV2):
            CsizVer2[icl] = ClusizV2[icl]
            TnpeVer2[icl] = ClunpeV2[icl]
            FstrVer2[icl] = FstripV2[icl]
        if NcluV2==0: 
            CsizVer2[0] = 0
            TnpeVer2[0] = 0
            dimVer2[0] = 1
            '''print ("NcluV2 =0")
            print ("CsizV2: ", CsizVer2[0], "Tnpe: ", TnpeVer2[0])'''

        OutTree.Fill()               
        
    ###########fine ciclo sulle hit del file jl
    
    tfile.Close()
   
    canvas[0].cd(2)
    OutTree.Draw("TpeStripHor0:energy>>hl0y(991,90,10000,100,-0.5,100)","","colz")
    TH2D(htemp) = gPad.GetPrimitive("htemp")
    htemp.GetXaxis().SetTitle("Energy keV")
    htemp.GetYaxis().SetTitle("Total Photo Electrons")
    canvas[0].cd(1)
    OutTree.Draw("TpeStripVer0:energy>>hl0x(991,90,10000,100,-0.5,100)","","colz")
    canvas[0].cd(4)
    OutTree.Draw("CsizHor0:energy>>hl0yclu(991,90,10000,10,-0.5,10)","PEstripHor0 >= 3","colz")
    canvas[0].cd(3)
    OutTree.Draw("CsizVer0:energy>>hl0xclu(991,90,10000,10,-0.5,10)","PEstripVer0 >= 3","colz")

    gStyle.SetPalette(55)
    
    canvas[1].cd(2)
    OutTree.Draw("TpeStripHor1:energy>>hl1y(991,90,10000,100,-0.5,100","","colz")
    canvas[1].cd(1)
    OutTree.Draw("TpeStripVer1:energy>>hl1x(991,90,10000,100,-0.5,100","","colz")
    canvas[1].cd(4)
    OutTree.Draw("CsizHor1:energy>>hl1yclu(991,90,10000,10,-0.5,10","PEstripHor1 >= 3","colz")
    canvas[1].cd(3)
    OutTree.Draw("CsizVer1:energy>>hl1xclu(991,90,10000,10,-0.5,10","PEstripVer1 >= 3","colz")
    
    canvas[2].cd(2)
    OutTree.Draw("TpeStripHor2:energy>>hl2y(991,90,10000,100,-0.5,100","","colz")
    canvas[2].cd(1)
    OutTree.Draw("TpeStripVer2:energy>>hl2x(991,90,10000,100,-0.5,100","","colz")
    canvas[2].cd(4)
    OutTree.Draw("CsizHor2:energy>>hl2yclu(991,90,10000,10,-0.5,10","PEstripHor2 >= 3","colz")
    canvas[2].cd(3)
    OutTree.Draw("CsizVer2:energy>>hl2xclu(991,90,10000,10,-0.5,10","PEstripVer2 >= 3","colz")
        
    cden.cd()
    hEkDist.Draw("")
    
    cnum.cd()
    color = [ROOT.kGreen+1, 4, ROOT.kRed-6, ROOT.kOrange-2, ROOT.kCyan-5, 28]
    label = ["Layer 0 View Y", "Layer 0 View X", "Layer 1 View Y", "Layer 1 View X", "Layer 2 View Y", "Layer 2 View X"]
    for ill in range(nview*totLayers):
        if ill == 0:
            hNTrigViewVsEk[ill].Draw("")
        else: hNTrigViewVsEk[ill].Draw("same")
        cnum.Update()

    cef.cd()
    cef.SetGrid()
    hEffViewVsEk[0].SetTitle("Trigger efficiency vs Ek")
    hEffViewVsEk[0].SetStats(0)
    legend = TLegend(0.66,.32,.89,.59)
    legend.SetBorderSize(1)
   
    for ill in range(nview*totLayers):
        hEffViewVsEk[ill].Divide(hNTrigViewVsEk[ill],hEkDist,nek,nek,"B")
        
        #label = str(ill)
        hEffViewVsEk[ill].GetYaxis().SetTitle("Trigger efficiency")
        hEffViewVsEk[ill].GetYaxis().SetRangeUser(0, 1.05)
        hEffViewVsEk[ill].GetYaxis().SetTitleFont(22)   
        hEffViewVsEk[ill].GetXaxis().SetTitle("Energy [keV]")
        hEffViewVsEk[ill].GetXaxis().SetTitleFont(22)       
        
        hEffViewVsEk[ill].SetLineColor(color[ill])
        hEffViewVsEk[ill].SetLineWidth(1)
        hEffViewVsEk[ill].SetMarkerColor(color[ill])
        hEffViewVsEk[ill].SetMarkerStyle(20)
        if ill == 0:
            hEffViewVsEk[ill].Draw("")
        else: hEffViewVsEk[ill].Draw("same")
        cef.Update()
        
        legend.AddEntry(hEffViewVsEk[ill],label[ill],"lp")
        cef.Update()
    legend.Draw();   
        
    
    fout.cd()
    cnum.Write()
    cden.Write()
    cef.Write()
    canvas[0].Write()
    canvas[1].Write()
    canvas[2].Write()
    OutTree.Write()
    fout.Close()
