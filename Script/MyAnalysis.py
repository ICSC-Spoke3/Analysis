import sys, os
import numpy as np
import array as ary
from ROOT import *

def InitGeom():
    
    # Riempio un vettore con le coordinate dei centri delle fibre e un vettore di flag
    # FiberFlag = 0 ==> Fibra nella vista xz 
    # FiberFlag = 1 ==> Fibra nella vista yz
    # FiberXY = coordinata x (o y) del centro della fibra
    # FiberZ = coordinata z del centro della fibra

    NumberOfLayers = 4
    FibersPerLayer = 256
    FiberRadius = 250.0e-4
    NumberOfFibers = NumberOfLayers*FibersPerLayer
    
    hh = 0.5 * FiberRadius * np.power(3, 0.5) 

    FiberFlag = np.zeros(NumberOfFibers)
    FiberXY = np.zeros(NumberOfFibers)
    FiberZ = np.zeros(NumberOfFibers)

    for ilayer in range (NumberOfLayers):
        for ifiber in range (FibersPerLayer):
            idf = FibersPerLayer*ilayer + ifiber
            if (ilayer == 0):
                FiberFlag[idf] = 0
                FiberZ[idf] = -0.5*FiberRadius - hh
                FiberXY[idf] = 2*FiberRadius*(0.5 + float(ifiber))
            elif (ilayer == 1):
                FiberFlag[idf] = 0
                FiberZ[idf] = -0.5*FiberRadius 
                FiberXY[idf] = 2*FiberRadius*(float(ifiber))
            elif (ilayer == 2):
                FiberFlag[idf] = 1
                FiberZ[idf] = 0.5*FiberRadius 
                FiberXY[idf] = 2*FiberRadius*(0.5 + float(ifiber))
            else:
                FiberFlag[idf] = 1
                FiberZ[idf] = 0.5*FiberRadius + hh
                FiberXY[idf] = 2*FiberRadius*(float(ifiber))

    # Riempio un secondo vettore con le posizioni degli estremi delle strip
    # StripXYcenter = x or y coordinate of the center of the strip
    # StripFlag = flag which says if the strip is horizontal or vertical (0 or 1)

    NumberOfStripViews = 2
    StripsPerView = 520
    NumberOfStrips = NumberOfStripViews*StripsPerView 
    StripSize = 250.0e-4
    StripOffset = -2.5*StripSize 

    StripFlag = np.zeros((NumberOfStrips),dtype=int)
    StripXYcenter = np.zeros(NumberOfStrips)
    
    for iview in range (NumberOfStripViews):
        for istrip in range (StripsPerView):
            idstrip = StripsPerView*iview + istrip
            if (iview == 0):
                StripFlag[idstrip] = 0
            else:
                StripFlag[idstrip] = 1
            StripXYcenter[idstrip] = StripOffset + StripSize*(float(istrip)+0.5)

    # Intersezioni fibre-strip
    FiberNStrips = np.zeros((NumberOfFibers),dtype=int)
    FiberIdStrips = np.zeros((NumberOfFibers,4),dtype=int)
    FiberArea = np.zeros((NumberOfFibers,4))
    FiberFracArea = np.zeros((NumberOfFibers,4))

    afib = np.pi*np.power(FiberRadius, 2)

    for ifib in range (NumberOfFibers):
        flag1 = FiberFlag[ifib]
        x1 = FiberXY[ifib]
        for istrip in range (NumberOfStrips):
            flag2 = StripFlag[istrip]
            x2 = StripXYcenter[istrip]
            if (flag1==flag2 and np.abs(x1-x2)<FiberRadius+0.5*StripSize):
                aa = x1 - FiberRadius
                bb = x1 + FiberRadius
                cc = x2 - 0.5*StripSize
                dd = x2 + 0.5*StripSize
                if (aa < cc):
                    xmin = cc
                    xmax = np.minimum(bb,dd)
                else:
                    xmin = aa
                    xmax = np.minimum(bb,dd)
                cth1 = (xmin-FiberXY[ifib])/FiberRadius
                if (cth1 > 1):
                    cth1 = 1
                if (cth1 < -1):
                    cth1 = -1
                cth2 = (xmax-FiberXY[ifib])/FiberRadius
                if (cth2 > 1):
                    cth2 = 1
                if (cth2 < -1):
                    cth2 = -1               
                th1 = np.arccos(cth1)
                th2 = np.arccos(cth2)
                s1 = np.sin(2*th1)
                s2 = np.sin(2*th2)
                area = np.power(FiberRadius,2)*(th1-th2+0.5*s2-0.5*s1)
                #print aa, bb, cc, dd, xmin, xmax, np.rad2deg(th1), np.rad2deg(th2), s1, s2, area
                # Check to avoid strips with very small area
                if (area > 1.0e-20):
                    FiberIdStrips[ifib][FiberNStrips[ifib]] = istrip
                    FiberArea[ifib][FiberNStrips[ifib]] = area
                    FiberNStrips[ifib] += 1
                    
        
        print ("Fiber ",ifib," Number of strips associated = ",FiberNStrips[ifib])
        for j in range (FiberNStrips[ifib]):
            FiberFracArea[ifib][j] = FiberArea[ifib][j] / afib
            print "Strip n.",FiberIdStrips[ifib][j]," area = ",FiberArea[ifib][j]," fraction = ",FiberFracArea[ifib][j]

    return (NumberOfStrips, StripsPerView, StripFlag, FiberNStrips, FiberIdStrips, FiberArea, FiberFracArea)


# MAIN PROGRAM

if __name__ == '__main__':

    # 8 ph/keV, trapping efficiency=5.4%, PDE = 0.4
    phkev = 8.*0.92 # including cladding
    trap = 5.4e-2
    PDE = 0.4

    nphemin = 4 # Numero di fotoelettroni di soglia nelle due viste

    ran = TRandom3()

    # NumberOfStrips = numero totale di strip
    # StripFlag[istr] = flag che dice se la strip e' in vista x (0) o vista y (1)
    # FiberNStrips[ifib] = numero di strip associate alla fibra ifib 
    # FiberIdStrips[ifib][j] (j in range FiberNStrips[ifib])  = vettore con gli identificativi delle strip associate alla fibra ifib 
    # FiberArea[ifib][j] (j in range FiberNStrips[ifib])  = vettore con le aree della fibra ifib viste da ciascuna delle strip ad essa associate  
    # FiberFracArea[ifib][j] (j in range FiberNStrips[ifib])  = vettore con le frazioni di area della fibra ifib viste da ciascuna delle strip ad essa associate  

    (NumberOfStrips,StripsPerView,StripFlag,FiberNStrips,FiberIdStrips,FiberArea,FiberFracArea) = InitGeom()
    

    # Define data files
    fdat = ["FibModule_500um_2ribbon00*_Out_Ele.root"]
    nfiles = len(fdat)
    print ("Number of files to be analyzed = ",nfiles)
    
    # Analysis trees
    te = TChain("CFTree")
    th = TChain("CFHits")
    for ifile in range (nfiles):
        te.Add(fdat[ifile])
        th.Add(fdat[ifile])

    nevt = te.GetEntries()
    print ("Number of events in the tree = ",nevt)

    # Open root file
    froot = TFile("MyAnalysis.root", "RECREATE")
    froot.cd()

    # Histogram definition
    hEstrip = TProfile("hEstrip", "Average energy (keV) vs strip ID", NumberOfStrips, -0.5, NumberOfStrips-0.5)
    hAvePhElestrip = TProfile("hAvePhElestrip", "Average number of photoelectrons vs strip ID", NumberOfStrips, -0.5, NumberOfStrips-0.5)
    
    hEstripDist = TH1D("hEstripDist", "Distribution of energy releases (keV) in the hit strips", 200, 0., 200.)
    hPhEelestripDist = TH1D("hPhElestripDist", "Distribution of photoelectrons in the hit strips", 51, -0.5, 50.5)

    ekmin = 0.0
    ekmax = 2000.0
    nek = 40
    hEkDist = TH1D("hEkDist", "Distribution of kinetic energies (keV)", nek, ekmin, ekmax)


    hPhEleViewDist = []
    hNTrigViewVsEk = []
    hEffViewVsEk = []
    for i in range (2):
        title1 = "hPhEleViewDist[" + str(i) + "]"
        title2 = "Distribution of photoelectrons in view " + str(i)
        hPhEleViewDist.append(TH1D(title1, title2, 61, -0.5, 60.5))
        title1 = "hNTrigViewVsEk[" + str(i) + "]"
        title2 = "Number of triggers vs kinetic energy - view " + str(i)
        hNTrigViewVsEk.append(TH1D(title1, title2, nek, ekmin, ekmax))
        title1 = "hEffViewVsEk[" + str(i) + "]"
        title2 = "Trigger efficiency vs kinetic energy - view " + str(i)
        hEffViewVsEk.append(TH1D(title1, title2, nek, ekmin, ekmax))


    mytree = TTree("mytree", "Summary")
    uekin = ary.array("f", [0])
    mytree.Branch("ekin", uekin, "ekin/F")
    uidmax = ary.array("i", 2*[0])
    mytree.Branch("idmax[2]", uidmax, "idmax[2]/I")
    utrig = ary.array("i", 2*[0])
    mytree.Branch("trig[2]", utrig, "trig[2]/I")
    ucharge = np.zeros((2,9),dtype=float)
    mytree.Branch("charge[2][9]", ucharge, "charge[2][9]/D")
    
    #nevt = 1000
    # Event loop
    for ievt in range (nevt):
        if (fmod(ievt,1000)==0):
            print "Now at event ",ievt,"/",nevt
        te.GetEntry(ievt)
        ekin = te.Ekevt*1.e6
        hEkDist.Fill(ekin, 1.0)
        nfib = te.NFib
        efib = te.EFib
        idfib = te.IDFib
        estrip = np.zeros(NumberOfStrips)
        for ifib in range (nfib):
            # Energy in keV
            ene = efib[ifib]*1e6
            # Ripartisco l'energia nella fibra tra le strip in base all'area
            jfib = idfib[ifib]
            nst = FiberNStrips[jfib]
            for ist in range (nst):
                jst = FiberIdStrips[jfib][ist]
                estrip[jst] += FiberFracArea[jfib][ist]*ene
        
        # Fill histos
        vnph = np.zeros(NumberOfStrips)
        vnph1 = np.zeros(2)
        vtrig = np.zeros(2)
        for ist in range (NumberOfStrips):
            hEstrip.Fill(float(ist), estrip[ist], 1.0)
            muph = 0
            if (estrip[ist] > 0.):
                hEstripDist.Fill(estrip[ist], 1.0)
                # Convert energy into Photoelectrons
                muph = estrip[ist]*phkev*trap*PDE
                nph = ran.Poisson(muph)
                hPhElestripDist.Fill(nph, 1.0)
                vnph[ist] = nph
                vnph1[StripFlag[ist]] += nph

            # Inside the strip loop, outside the if...
            hAvePhElestrip.Fill(float(ist), muph, 1.0)

        # Outside the strip loop...
        for iview in range (2):
            hPhEleViewDist[iview].Fill(vnph1[iview], 1.0)

            # Check trigger condition
            istrip1 = StripsPerView*iview
            istrip2 = istrip1 + StripsPerView
            for istrip in range (istrip1, istrip2-1):
                if (vnph[istrip] >= nphemin and vnph[istrip+1] >= nphemin):
                    vtrig[iview] = 1

            if (vtrig[iview] > 0):
                hNTrigViewVsEk[iview].Fill(ekin, 1.0)

        # Fill the tree    
        uekin[0] = ekin
        for iview in range (2):
            istrip1 = StripsPerView*iview
            istrip2 = istrip1 + StripsPerView
            qmax = -1.0
            idmax = -999
            for istrip in range (istrip1, istrip2):
                if (vnph[istrip] >= qmax):
                    idmax = istrip
                    qmax = vnph[istrip]
                                        
            uidmax[iview] = idmax
            #print idmax, qmax
            utrig[iview] = int(vtrig[iview])
            kmin = np.maximum(istrip1, idmax-4)
            kmax = np.minimum(istrip2, idmax+5)
            if (idmax < 0):
                for k1 in range (9):
                    ucharge[iview][k1] = 0.0
            else:
                for k in range (kmin, kmax):
                    k1 = k - idmax + 4
                    ucharge[iview][k1] = vnph[k]
                    #print k,k1,vnph[k], ucharge[iview][k1]

        mytree.Fill()

    # Final operations on histos
    for iview in range (2):
        for k in range (nek):
            xx = hNTrigViewVsEk[iview].GetBinContent(k+1)
            yy = hEkDist.GetBinContent(k+1)
            zz = 0.
            dzz = 0.
            if (yy > 0.):
                zz = xx/yy
                dzz = np.power(zz*(1-zz)/yy, 0.5)
            hEffViewVsEk[iview].SetBinContent(k+1, zz)
            hEffViewVsEk[iview].SetBinError(k+1, dzz)


    # Close and write root file
    froot.cd()
    froot.Write()
    froot.Close()
