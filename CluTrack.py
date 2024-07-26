from ROOT import TFile, TTree, TChain, TH1D, TH2D, TMath, TCanvas, TRandom3, TColor, TGraph, TGraphErrors, TF1, TH1F, TObject, gStyle, TLegend
import os
import sys
import time
import array as ary
import numpy as np
import glob
ran = TRandom3()

def evaluateTrackPoint(x0,z0,cx,cz,zP,flag):
    t = (zP-z0)/cz
    x = x0 + cx*t
    if(flag==1):
        print( "t = ",t,"xTr = ",x)
    return x

def evaluateRecoCosineDirector(Mx,My):
    a0 = pow(Mx,2.)
    a1 = pow(My,2.)
    crz = -1/(pow(a0+a1+1,0.5))
    crx = crz*Mx
    cry = crz*My
    return crx,cry,crz

def ScalarProduct(cx,cy,cz,crx,cry,crz):
    ax = cx*crx
    ay = cy*cry
    az = cz*crz
    a = ax+ay+az
    return a

def Clustering(hits, charges):
    nhits = len(hits)
    
    if nhits > 0:
        #print("biii")
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
    
    else: 
        nclu=0
        istrip0 = [-1]
        nstrip0 = [0]
        charge0 = [0]    

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

def read_GeoFile(fgeo):
    f = open(fgeo,'r')
    lines = len(open(fgeo).readlines())

    overall_track_geom = np.zeros(7) # n.of Layer,Views,Fibers, FibLength, FibRadius, TrackerLength, TrackerWidth
    overall_cal_geom   = np.zeros(6) # Pixel Size, Depth, nPix_per_view, TotPixels, CaloLength, CaloWidth
    itr = 0
    ical = 0
    for i in range(lines):
        line = f.readline()
        aaa = line.split()
        Detector = str(aaa[0])
        if(Detector == "TRACKER"):
            if(itr==0):
                layer  = int(aaa[1])
                view   = int(aaa[2])
                fiber  = int(aaa[3])
                lfib   = float(aaa[4])
                rfib   = float(aaa[5])
                TrackerLength = float(aaa[6])
                TrackerWidth  = float(aaa[7])

                overall_track_geom[0] = layer
                overall_track_geom[1] = view
                overall_track_geom[2] = fiber
                overall_track_geom[3] = lfib
                overall_track_geom[4] = rfib
                overall_track_geom[5] = TrackerLength
                overall_track_geom[6] = TrackerWidth

                gxcfib = np.zeros((layer,view,fiber))
                gycfib = np.zeros((layer,view,fiber))
                gzcfib = np.zeros((layer,view,fiber))
            else:
                layer = int(aaa[1])
                view  = int(aaa[2])
                fiber = int(aaa[3])
                xc = float(aaa[4])
                yc = float(aaa[5])
                zc = float(aaa[6])
                if(view==0):
                    gxcfib[layer][view][fiber] = xc
                    gycfib[layer][view][fiber] = yc
                    gzcfib[layer][view][fiber] = zc
                else:
                    gxcfib[layer][view][fiber] = yc
                    gycfib[layer][view][fiber] = xc
                    gzcfib[layer][view][fiber] = zc
            itr +=1

        if(Detector=="CALORIMETER"):
            if(ical==0):
                pixSize = float(aaa[1])
                pixDepth = float(aaa[2])
                nPix = int(aaa[3])
                totPix = nPix*nPix
                CaloLength = float(aaa[4])
                CaloWidth  = float(aaa[5])

                overall_cal_geom[0] = pixSize
                overall_cal_geom[1] = pixDepth
                overall_cal_geom[2] = nPix
                overall_cal_geom[3] = totPix
                overall_cal_geom[4] = CaloLength
                overall_cal_geom[5] = CaloWidth

                gxcpix = np.zeros(totPix)
                gycpix = np.zeros(totPix)
                gzcpix = np.zeros(totPix)

            else:
                iXPix = int(aaa[1])
                iYPix = int(aaa[2])
                iPix  = int(aaa[3])
                xcPix = float(aaa[4])
                ycPix = float(aaa[5])
                zcPix = float(aaa[6])

                gxcpix[iPix] = xcPix
                gycpix[iPix] = ycPix
                gzcpix[iPix] = zcPix

            ical +=1

    f.close()
    return overall_track_geom, overall_cal_geom,gxcfib,gycfib,gzcfib,gxcpix,gycpix,gzcpix

if __name__ == '__main__':
    print('something goes wrong 1')
    #InputFile
    try:
        fin = sys.argv[1]
        #fin     = "/home/nadia_work/Scrivania/NUSESe-Pow_0.1-5000_1000000-evt.root"
        #fgeo = "/home/nadia_work/Scrivania/NUSESe-Pow_0.1-5000_1000000-evt.txt"
        fgeo = sys.argv[2]
        
    except:
        print("ERROR MESSAGE: \n =====> Insert first argument = ROOTfile.root and second argument = Geofile.txt <=====")
        sys.exit(1)
        
    tfile   = TFile(fin)

    #OutTreeFile
    OutputName = "Tree_"+(fin.split("/")[-1]).split(".r")[0] + ".root"
    OutputFile = OutputName
    root_file = TFile(OutputFile, "RECREATE")

    OutTree           = TTree("ClusterTree", "ClusterTree")
    BigClusterTree    = TTree("BigClusterTree", "BigClusterTree")

    Track_info,Calo_info,gxcfib,gycfib,gzcfib,gxcpix,gycpix,gzcpix = read_GeoFile(fgeo)
    Layers               = Track_info[0]
    Views                = Track_info[1]
    Fibers               = Track_info[2]
    FibLength            = Track_info[3]
    FibRadius            = Track_info[4]
    TrackerLength        = Track_info[5]
    TrackerWidth         = Track_info[6]
    print ("")
    print ("\t --- Tracker (DIM",TrackerLength,"x",TrackerLength,"x",TrackerWidth,"mm^3) ---")
    print ("[TRACKER_GEO_INFO]: Number of Layers in the Tracker = ", Layers)
    print ("[TRACKER_GEO_INFO]: Number of Views per Layer = ", Views)
    print ("[TRACKER_GEO_INFO]: Number of Fibers per View = ", Fibers)
    print ("[TRACKER_GEO_INFO]: Fiber Length = ",FibLength,"mm","Fiber Radius = ",FibRadius,"mm")
    #print ("x",gxcfib[0][0])
    #print ("y",gycfib[0][1])
    #print ("z",gzcfib[0])

    nview = int(Views)
    totLayers = int(Layers)
    nfiber  = int(Fibers)
    rfiber  = FibRadius #mm
    pitch   = 0.250 #mm
    stripNo = int((nfiber*rfiber+rfiber)/pitch) #128
    xoff    = 0 #rfiber*0.5
    yoff    = 0 #rfiber*0.5
    sigma_x = pitch/2 #pitch/np.sqrt(12)
    ################################
    #ho 128 strip per l'elettronica di lettura: da 0 a 127
    #quindi l'ultima fibra e' letta dalla 127 fibra
    #le fibre sono 130 da 0 a 129 quindi l'ultima fibra non viene letta
    #Complessivamente NON leggiamo meta' della fibra zero, meta' della fibra 128 e la fibra 129
    ################################
    fiberyield = 8*0.92#kev/pe (con cladding)
    PDE = 0.4
    trapeff = 0.054
    PEthresh = 3 #pe per strip

    # Define output tree
    nChan = 400000 
    eventID = ary.array('i',[0])
    Tch = ary.array('i',[0])
    energy  = ary.array('f',[0.])
    Ichip = ary.array('i',range(nChan))
    Ichan = ary.array('i',range(nChan))
    Ipix = ary.array('i',range(nChan))
    Npe = ary.array('f',range(nChan))

    OutTree.Branch("eventID", eventID, "eventID/I")
    OutTree.Branch("Tch", Tch, "Tch/I")
    OutTree.Branch("energy", energy, "energy_keV/F")
    OutTree.Branch("Ichip", Ichip, "Ichip[Tch]/I")
    OutTree.Branch("Ichan", Ichan, "Ichan[Tch]/I")
    OutTree.Branch("Ipix", Ipix, "Ipix[Tch]/I")
    OutTree.Branch("Npe", Npe, "Npe[Tch]/F")

    # Informazioni generali sui cluster
    nClusters = 1000 
    Tclu_tot = ary.array('i',[0])
    dim_tot = ary.array('i',[0])
    Csiz_tot = ary.array('i',range(nClusters))
    Tnpe_tot = ary.array('f',range(nClusters))
    Fstr_tot = ary.array('i',range(nClusters))

    OutTree.Branch("Tclu_tot", Tclu_tot, "Tclu_tot/I")
    OutTree.Branch("dim_tot", dim_tot, "dim_tot/I")
    OutTree.Branch("Csiz_tot", Csiz_tot, "Csiz_tot[dim_tot]/I")
    OutTree.Branch("Tnpe_tot", Tnpe_tot, "Tnpe_tot[dim_tot]/F")
    OutTree.Branch("Fstr_tot", Fstr_tot, "Fstr_tot[dim_tot]/I")

    # Informazioni sui cluster nelle viste orizzontale (x=0) e verticale(y=1) per i 3 layer
    for il in range(totLayers):    
        for iv in range(nview):
            if iv == 0: 
                v = "Hor"
            else:
                v="Ver"    
            sTstrip    = "Tstrip{0}{1}".format(v,il)
            sTstrip2   = "Tstrip{0}{1}/I".format(v,il)
            sTpeStrip  = "TpeStrip{0}{1}".format(v,il)
            sTpeStrip2 = "TpeStrip{0}{1}/F".format(v,il)
            sstrip     = "strip{0}{1}".format(v,il)
            sstrip2    = "strip{0}{1}[{2}]/I".format(v,il,sTstrip)
            sPEstrip   = "PEstrip{0}{1}".format(v,il)
            sPEstrip2  = "PEstrip{0}{1}[{2}]/F".format(v,il,sTstrip)

            sdim       = "dim{0}{1}".format(v,il)
            sdim2      = "dim{0}{1}/I".format(v,il)
            sTclu      = "Tclu{0}{1}".format(v,il)
            sTclu2     = "Tclu{0}{1}/I".format(v,il)
            sCsiz      = "Csiz{0}{1}".format(v,il)
            sCsiz2     = "Csiz{0}{1}[{2}]/I".format(v,il,sdim)
            sTnpe      = "Tnpe{0}{1}".format(v,il)
            sTnpe2     = "Tnpe{0}{1}[{2}]/F".format(v,il,sdim)
            sFstr      = "Fstr{0}{1}".format(v,il)
            sFstr2     = "Fstr{0}{1}[{2}]/I".format(v,il,sdim)

            exec("%s = ary.array('i',[0])" % sTstrip)
            exec("%s = ary.array('f',[0])" % sTpeStrip)
            exec("%s = ary.array('i',range(nChan))" % sstrip)
            exec("%s = ary.array('f',range(nChan))" % sPEstrip)

            exec("%s = ary.array('i',[0])" % sdim)
            exec("%s = ary.array('i',[0])" % sTclu)
            exec("%s = ary.array('i',range(nClusters))" % sCsiz)
            exec("%s = ary.array('f',range(nClusters))" % sTnpe)
            exec("%s = ary.array('i',range(nClusters))" % sFstr)     

            OutTree.Branch(sTstrip, eval(sTstrip),sTstrip2)
            OutTree.Branch(sTpeStrip, eval(sTpeStrip), sTpeStrip2)
            OutTree.Branch(sstrip, eval(sstrip), sstrip2)
            OutTree.Branch(sPEstrip, eval(sPEstrip), sPEstrip2)

            OutTree.Branch(sdim, eval(sdim), sdim2)
            OutTree.Branch(sTclu, eval(sTclu), sTclu2)
            OutTree.Branch(sCsiz, eval(sCsiz), sCsiz2)
            OutTree.Branch(sTnpe, eval(sTnpe), sTnpe2)
            OutTree.Branch(sFstr, eval(sFstr), sFstr2)

    # Informazioni sul tracciamento
    nlayers = ary.array('i',[0])
    trk_flag = ary.array('i',[0])
    x_vec = ary.array('f',range(totLayers))
    dx_vec = ary.array('f',range(totLayers))
    y_vec = ary.array('f',range(totLayers))
    dy_vec = ary.array('f',range(totLayers))
    zx_vec = ary.array('f',range(totLayers))
    dzx_vec = ary.array('f',range(totLayers))
    zy_vec = ary.array('f',range(totLayers))
    dzy_vec = ary.array('f',range(totLayers))
    cosx = ary.array('f',range(1))
    cosy = ary.array('f',range(1))
    cosz = ary.array('f',range(1))
    
    qx = ary.array('f',range(1))
    mx = ary.array('f',range(1))
    chi2x = ary.array('f',range(1))
    x_residual = ary.array('f',range(totLayers))
    qy = ary.array('f',range(1))
    my = ary.array('f',range(1))
    chi2y = ary.array('f',range(1))
    y_residual = ary.array('f',range(totLayers))

    xMC_vec = ary.array('f',range(totLayers))
    yMC_vec = ary.array('f',range(totLayers))
    zxMC_vec = ary.array('f',range(totLayers))
    zyMC_vec = ary.array('f',range(totLayers))
    cosxMC =  ary.array('f',range(1))
    cosyMC =  ary.array('f',range(1))
    coszMC =  ary.array('f',range(1))

    angle_deg = ary.array('f',range(1))

    OutTree.Branch("nlayers", nlayers, "nlayers/I")
    OutTree.Branch("trk_flag", trk_flag, "trk_flag/I")
    OutTree.Branch("x_vec", x_vec, "x_vec[nlayers]/F")
    OutTree.Branch("dx_vec", dx_vec, "dx_vec[nlayers]/F")
    OutTree.Branch("y_vec", y_vec, "y_vec[nlayers]/F")
    OutTree.Branch("dy_vec", dy_vec, "dy_vec[nlayers]/F")
    OutTree.Branch("zx_vec", zx_vec, "zx_vec[nlayers]/F")
    OutTree.Branch("dzx_vec", dzx_vec, "dzx_vec[nlayers]/F")
    OutTree.Branch("zy_vec", zy_vec, "zy_vec[nlayers]/F")
    OutTree.Branch("dzy_vec", dzy_vec, "dzy_vec[nlayers]/F")
    OutTree.Branch("cosx", cosx, "cosx[trk_flag]/F")
    OutTree.Branch("cosy", cosy, "cosy[trk_flag]/F")
    OutTree.Branch("cosz", cosz, "cosz[trk_flag]/F")

    OutTree.Branch("qx", qx, "qx[trk_flag]/F")
    OutTree.Branch("mx", mx, "mx[trk_flag]/F")
    OutTree.Branch("chi2x", chi2x, "chi2x[trk_flag]/F")
    OutTree.Branch("x_residual", x_residual, "x_residual[nlayers]/F")
    OutTree.Branch("qy", qy, "qy[trk_flag]/F")
    OutTree.Branch("my", my, "my[trk_flag]/F")
    OutTree.Branch("chi2y", chi2y, "chi2y[trk_flag]/F")
    OutTree.Branch("y_residual", y_residual, "y_residual[nlayers]/F")

    OutTree.Branch("xMC_vec", xMC_vec, "xMC_vec[nlayers]/F")
    OutTree.Branch("yMC_vec", yMC_vec, "yMC_vec[nlayers]/F")
    OutTree.Branch("zxMC_vec", zxMC_vec, "zxMC_vec[nlayers]/F")
    OutTree.Branch("zyMC_vec", zyMC_vec, "zyMC_vec[nlayers]/F")
    OutTree.Branch("cosxMC", cosxMC, "cosxMC[trk_flag]/F")
    OutTree.Branch("cosyMC", cosyMC, "cosyMC[trk_flag]/F")
    OutTree.Branch("coszMC", coszMC, "coszMC[trk_flag]/F")

    OutTree.Branch("angle_deg", angle_deg, "angle_deg[trk_flag]/F")

    tree = tfile.Get("TrackerDigi")
    ptree = tfile.Get("Primary")
    nhit = tree.GetEntries()
    #nhit = 1

    root_file.cd()

    for i in range(nhit):
        eventID[0] = i
        if i%1000==0:
            print("Event: ", eventID[0])
        #print("Event: ", eventID[0])  
        #event[0]   = i
        #--- Primary Track ---#
        ptree.GetEntry(i)
        x0MC = ptree.PrimaryParticlePositionX
        y0MC = ptree.PrimaryParticlePositionY
        z0MC = ptree.PrimaryParticlePositionZ
        cx = ptree.PrimaryParticleDirectionX
        cy = ptree.PrimaryParticleDirectionY
        cz = ptree.PrimaryParticleDirectionZ
        #lettura DigiTrackerTree 
        tree.GetEntry(i) 
        xTMC = tree.TrackParticlePositionX
        yTMC = tree.TrackParticlePositionY
        zTMC = tree.TrackParticlePositionZ
        PrimaryEnergy = (tree.PrimaryParticleEnergy)*1e3
        dE_vec        = tree.Energy_keV
        view_vec      = tree.ViewNo
        layer_vec     = tree.LayerNo
        fiber_vec     = tree.FiberNo

        energy[0] = PrimaryEnergy

        Strip_vec   = []
        PEstrip_vec = []
        Strips      = []
        Strip_r     = []
        eStrips     = []
        PEStrips    = []
        z           = []
        totPE       = []
        String_zs   = []

        Lay  = []
        View = []

        Nclu_tot = 0
        Fstrip_tot = []
        Clusiz_tot = []
        Clunpe_tot = []

        Nclu   = []
        Fstrip = []
        Clusiz = []
        Clunpe = []

        #ccc.append([])    

        for il in range(totLayers):
            Strips.append([])
            Strip_r.append([])
            eStrips.append([])
            PEStrips.append([])
            z.append([])
            String_zs.append([])

            Nclu.append([])
            Fstrip.append([])
            Clusiz.append([])
            Clunpe.append([])

            totPE.append([])

            for iv in range(nview):
                Strips[il].append([])
                Strip_r[il].append([])
                eStrips[il].append([])
                PEStrips[il].append([])
                z[il].append([])
                String_zs[il].append([])

                Nclu[il].append(0)
                Fstrip[il].append([])
                Clusiz[il].append([])
                Clunpe[il].append([])

                totPE[il].append(0)

        for j in range(len(fiber_vec)):
            mpe = (dE_vec[j]*fiberyield*trapeff*PDE)
            npe1 = ran.Poisson(mpe/2)
            npe2 = ran.Poisson(mpe/2)

            #Lay.append(layer_vec[j])
            #View.append(view_vec[j])

            #Fibra all'estrema sx e estrema dx dei layer lette a meta' sull'upper ribble
            #l'ultima non e' letta
            #le centrali su ogni strato sono lette da 2 strip ciascuna:
            #dx
            if (fiber_vec[j] == 0):
                index1 = 0
                Strip_vec.append(index1) 
                PEstrip_vec.append(npe1)
                Lay.append(layer_vec[j])
                View.append(view_vec[j])
                for il in range(totLayers):
                    for iv in range(nview):
                        if (layer_vec[j] == il and view_vec[j] == iv):
                            if index1 not in Strips[il][iv]:
                                Strips[il][iv].append(index1)
                                if iv == 0:
                                    Strip_r[il][iv].append(gxcfib[il][iv][fiber_vec[j]]-pitch/2)
                                elif iv == 1:
                                    Strip_r[il][iv].append(gycfib[il][iv][fiber_vec[j]]-pitch/2)
                                z[il][iv].append(gzcfib[il][iv][fiber_vec[j]])
                                eStrips[il][iv].append(dE_vec[j]/2)
                                PEStrips[il][iv].append(npe1)
                            else:
                                a=np.where(np.array(Strips[il][iv])==index1)[0][0]
                                eStrips[il][iv][a]+=(dE_vec[j]/2)
                                PEStrips[il][iv][a]+=npe1
                            totPE[il][iv] += npe1
            #sx
            elif (fiber_vec[j] == stripNo): 
                index1 = (stripNo-1)
                Strip_vec.append(index1)
                PEstrip_vec.append(npe1)
                Lay.append(layer_vec[j])
                View.append(view_vec[j])
                for il in range(totLayers):
                    for iv in range(nview):
                        if (layer_vec[j] == il and view_vec[j] == iv):
                            if index1 not in Strips[il][iv]:
                                Strips[il][iv].append(index1)
                                if iv == 0:
                                    Strip_r[il][iv].append(gxcfib[il][iv][fiber_vec[j]]+pitch/2)
                                elif iv == 1:
                                    Strip_r[il][iv].append(gycfib[il][iv][fiber_vec[j]]+pitch/2)
                                z[il][iv].append(gzcfib[il][iv][fiber_vec[j]]) 
                                eStrips[il][iv].append(dE_vec[j]/2)
                                PEStrips[il][iv].append(npe1)
                            else:
                                a=np.where(np.array(Strips[il][iv])==index1)[0][0]
                                eStrips[il][iv][a]+=(dE_vec[j]/2)
                                PEStrips[il][iv][a]+=npe1    
                            totPE[il][iv] += npe1
            #fibre centrali
            elif(fiber_vec[j] > 0 and fiber_vec[j] < stripNo):
                index1 = fiber_vec[j]
                index2 = index1-1
                Strip_vec.append(index1)
                Strip_vec.append(index2)
                PEstrip_vec.append(npe1)
                PEstrip_vec.append(npe2)
                Lay.append(layer_vec[j])
                Lay.append(layer_vec[j])
                View.append(view_vec[j])
                View.append(view_vec[j])
                for il in range(totLayers):  
                    for iv in range(nview):
                        if (layer_vec[j] == il and view_vec[j] == iv):
                            if index1 not in Strips[il][iv]:
                                Strips[il][iv].append(index1)
                                if iv == 0:
                                    Strip_r[il][iv].append(gxcfib[il][iv][fiber_vec[j]]+pitch/2)
                                elif iv == 1:
                                    Strip_r[il][iv].append(gycfib[il][iv][fiber_vec[j]]+pitch/2)
                                z[il][iv].append(gzcfib[il][iv][fiber_vec[j]])
                                eStrips[il][iv].append(dE_vec[j]/2)
                                PEStrips[il][iv].append(npe1)
                            elif index1 in Strips[il][iv]:
                                a=np.where(np.array(Strips[il][iv])==index1)[0][0]
                                eStrips[il][iv][a]+=(dE_vec[j]/2)
                                PEStrips[il][iv][a]+=npe1    
                            if index2 not in Strips[il][iv]:
                                Strips[il][iv].append(index2)
                                if iv == 0:
                                    Strip_r[il][iv].append(gxcfib[il][iv][fiber_vec[j]]-pitch/2)
                                elif iv == 1:
                                    Strip_r[il][iv].append(gycfib[il][iv][fiber_vec[j]]-pitch/2)
                                z[il][iv].append(gzcfib[il][iv][fiber_vec[j]])
                                eStrips[il][iv].append(dE_vec[j]/2)
                                PEStrips[il][iv].append(npe2)
                            elif index2 in Strips[il][iv]:
                                a=np.where(np.array(Strips[il][iv])==index2)[0][0]
                                eStrips[il][iv][a]+=(dE_vec[j]/2)
                                PEStrips[il][iv][a]+=npe2   
                            totPE[il][iv] += npe1
                            totPE[il][iv] += npe2
            #ultima fibra saltata
            elif (fiber_vec[j] == nfiber -1): 
                continue
        
        ordered_index_list = []#np.argsort(Strips)
        ordered_strip=[]
        ordered_centre=[]
        ordered_energy=[]
        ordered_PE=[]

        for il in range(totLayers):
            ordered_index_list.append([])
            ordered_strip.append([])
            ordered_centre.append([])
            ordered_energy.append([])
            ordered_PE.append([])
            for iv in range(nview):
                ordered_index_list[il].append([])
                ordered_strip[il].append([])
                ordered_centre[il].append([])
                ordered_energy[il].append([])
                ordered_PE[il].append([])
                ordered_index_list[il][iv]=np.argsort(Strips[il][iv])
                for ind in range (len(ordered_index_list[il][iv])):
                    ordered_strip[il][iv].append(Strips[il][iv][ordered_index_list[il][iv][ind]])
                    ordered_centre[il][iv].append(Strip_r[il][iv][ordered_index_list[il][iv][ind]])
                    ordered_energy[il][iv].append(eStrips[il][iv][ordered_index_list[il][iv][ind]])
                    ordered_PE[il][iv].append(PEStrips[il][iv][ordered_index_list[il][iv][ind]])
        
        Strips   = ordered_strip
        Strip_r  = ordered_centre
        eStrips  = ordered_energy
        PEStrips = ordered_PE
        ######fine ciclo sulle hit dell'evento
        #print("strip_r0x",Strip_r[0][0])
        #print("strip_r0y",Strip_r[0][1])
        #print("strip_r1x",Strip_r[1][0])
        #print("strip_r1y",Strip_r[1][1])

        vtrig= np.zeros((totLayers,nview))
        Tch[0] = (len(PEstrip_vec))
        #Tch2[0] = (len(PEstrip_vec))
        nclu,istrip0,nstrip0,charge0 = Clustering(Strip_vec,PEstrip_vec)  
        Nclu_tot += nclu

        for kk in range(len(istrip0)):
            Fstrip_tot.append(istrip0[kk])
            Clusiz_tot.append(nstrip0[kk])
            Clunpe_tot.append(charge0[kk])
        for il in range(totLayers):
            for iv in range(nview):
                if iv == 0: 
                    v = "Hor"
                else:
                    v="Ver"    
                sTstrip    = "Tstrip{0}{1}".format(v,il)
                sTpeStrip  = "TpeStrip{0}{1}".format(v,il)
                sstrip     = "strip{0}{1}".format(v,il)
                sPEstrip   = "PEstrip{0}{1}".format(v,il)

                sdim       = "dim{0}{1}".format(v,il)
                sTclu      = "Tclu{0}{1}".format(v,il)
                sCsiz      = "Csiz{0}{1}".format(v,il)
                sTnpe      = "Tnpe{0}{1}".format(v,il)
                sFstr      = "Fstr{0}{1}".format(v,il)

                nclu,istrip0,nstrip0,charge0 = Clustering(Strips[il][iv],PEStrips[il][iv])
                Nclu[il][iv] = nclu
                for kk in range(len(istrip0)):
                    Fstrip[il][iv].append(istrip0[kk])
                    Clusiz[il][iv].append(nstrip0[kk])
                    Clunpe[il][iv].append(charge0[kk])

                exec("%s[0] =len(PEStrips[il][iv])" % (sTstrip))
                exec("%s[0] = totPE[il][iv]" % (sTpeStrip))
                exec("%s[0] = Nclu[il][iv]" % (sTclu))
                exec("%s[0] = Nclu[il][iv]" % (sdim))

                for istrip in range (eval(sTstrip)[0]):
                    exec("%s[istrip] = Strips[il][iv][istrip]" % (sstrip))
                    exec("%s[istrip] = PEStrips[il][iv][istrip]" % (sPEstrip))
                for icl in range(Nclu[il][iv]):
                    exec("%s[icl] = Clusiz[il][iv][icl]" % (sCsiz))
                    exec("%s[icl] = Clunpe[il][iv][icl]" % (sTnpe))
                    exec("%s[icl] = Fstrip[il][iv][icl]" % (sFstr))
                if Nclu[il][iv]==0: 
                    exec("%s[0] = 0" % (sCsiz))
                    exec("%s[0] = 0" % (sTnpe))
                    exec("%s[0] = 1" % (sdim))

        Tclu_tot[0] = Nclu_tot
        dim_tot[0] = Nclu_tot
        for ihit in range (Tch[0]):
            #ILay[ihit] = Lay[ihit]
            #IView[ihit] = View[ihit]
            #energy[ihit] = Emax
            Ipix[ihit] = Strip_vec[ihit]
            Npe[ihit] = PEstrip_vec[ihit]

        for icl in range(Nclu_tot):
            Csiz_tot[icl] = Clusiz_tot[icl]
            Tnpe_tot[icl] = Clunpe_tot[icl]
            Fstr_tot[icl] = Fstrip_tot[icl]
        if Nclu_tot==0: 
            Csiz_tot[0] = 0
            Tnpe_tot[0] = 0
            dim_tot[0] = 1                             
        
        #cerco il cluster di carica massima su ciascun piano:
        flag=0
        istrip_vec = [] 
        bc = []
        debc = []
        x = []
        y = []
        dx = []
        dy = []
        zx = []
        zy = []
        dzx = []
        dzy = []
        primary_track = []
        xtrack = []
        ytrack = []

        #ciclo sul baricentro di carica
        for il in range(totLayers):
            istrip_vec.append([])
            bc.append([])
            debc.append([])
            primary_track.append([])
            for iv in range(nview):
                bc[il].append([])
                debc[il].append([])
                primary_track[il].append([])
                istrip_vec[il].append([])
                if iv == 0: 
                    v = "Hor"
                else:
                    v="Ver"
                sTclu      = "Tclu{0}{1}".format(v,il)
                sTnpe      = "Tnpe{0}{1}".format(v,il)
                sstrip     = "strip{0}{1}".format(v,il)
                sFstr      = "Fstr{0}{1}".format(v,il)
                sCsiz      = "Csiz{0}{1}".format(v,il)

                if eval(sTclu)[0] > 0:
                    String_zs[il][iv] = np.mean(z[il][iv])
                    maxClu    = max(eval(sTnpe)[:eval(sTclu)[0]])
                    indexmax  = eval(sTnpe)[:eval(sTclu)[0]].index(maxClu)
                    fstripmax = eval(sFstr)[indexmax]
                    
                    indexfstripmax    = eval(sstrip).index(fstripmax)

                    for istr in range(eval(sCsiz)[indexmax]):
                        istrip_vec[il][iv].append(indexfstripmax+istr)

                    x0 = istrip_vec[il][iv][0]
                    x1 = istrip_vec[il][iv][-1]+1
                    num = sum(np.multiply(Strip_r[il][iv][x0:x1],PEStrips[il][iv][x0:x1]))
                    den = sum(PEStrips[il][iv][x0:x1])

                    if den >= 4:  #almeno 4 fotoelettroni nel cluster

                        de1 = np.array(PEStrips[il][iv][x0:x1])*sigma_x/den
                        sumde2 = np.array(Strip_r[il][iv][x0:x1])*den-np.array(Strip_r[il][iv][x0:x1])*np.array(PEStrips[il][iv][x0:x1])
                        de2 = np.array(sumde2) *pow(np.array(PEStrips[il][iv][x0:x1]),1/2)/pow(den,2)

                        bc[il][iv]= num/den
                        debc[il][iv] = np.sqrt(sum(pow(de1,2))+sum(pow(de2,2)))
                        primary_track[il][iv] = evaluateTrackPoint(y0MC,z0MC,cy,cz,String_zs[il][iv],flag)
                        if iv==0:
                            x.append(bc[il][iv])
                            zx.append(String_zs[il][iv])
                            dx.append(debc[il][iv])
                            dzx.append(2*rfiber)
                            xtrack.append(primary_track[il][iv])
                            
                            x_vec[il] = bc[il][iv] 
                            dx_vec[il] = debc[il][iv]
                            zx_vec[il] = String_zs[il][iv]
                            dzx_vec[il] = 2*rfiber 
                            xMC_vec[il] = primary_track[il][iv]
                            zxMC_vec[il] = String_zs[il][iv]

                        elif iv==1:
                            y.append(bc[il][iv])
                            zy.append(String_zs[il][iv])
                            dy.append(debc[il][iv])
                            dzy.append(2*rfiber)
                            ytrack.append(primary_track[il][iv])

                            y_vec[il] = bc[il][iv]
                            dy_vec[il] = debc[il][iv]
                            zy_vec[il] = String_zs[il][iv]
                            dzy_vec[il] = 2*rfiber
                            yMC_vec[il] = primary_track[il][iv]
                            zyMC_vec[il] = String_zs[il][iv]

        #reco vs MC track
        x = np.array(x).astype(float)
        xtrack = np.array(xtrack).astype(float)                                      
        dx = np.array(dx).astype(float)
        zx = np.array(zx).astype(float)
        dzx= np.array(dzx).astype(float)
        y = np.array(y).astype(float)
        dy = np.array(dy).astype(float)
        ytrack = np.array(ytrack).astype(float)
        zy = np.array(zy).astype(float)
        dzy= np.array(dzy).astype(float)
        x_res = []
        y_res = []

        if (np.shape(x)[0] == totLayers and np.shape(y)[0] == totLayers): 
            
            nlayers[0] = totLayers
            trk_flag[0] = 1
            
            x_reco = TGraphErrors(totLayers,zx,x,dzx,dx)
            x_reco.SetTitle("Hor")
            x_true = TGraph(3,zx,xtrack)
            #x_true = TGraph(len(xTMC),zTMC,xTMC) #posizione MC layer per layer
            print('something goes wrong')
            y_reco = TGraphErrors(totLayers,zy,y,dzy,dy)
            y_reco.SetTitle("Ver")
            y_true = TGraph(3,zy,ytrack)  
            ##y_true = TGraph(len(xTMC),zTMC,yTMC)  
            for iv in range(nview):
                #ccc[i].append(TCanvas("ccc{0}_{1}".format(iv,i),"ccc{0}_{1}".format(iv,i)))
                #ccc[i][iv].cd()
                if iv == 0:
                    xmin = min(zx)-0.05
                    xmax = max(zx)+0.05
                    rettax = TF1("rettax","[0]*x+[1]",xmin,xmax)
                    rettaMC = TF1("rettaMC","[0]*x+[1]",min(xtrack),max(xtrack))
                    x_reco.SetMarkerStyle(8)
                    #x_reco.GetXaxis().SetLimits(-1,1)
                    x_reco.GetYaxis().SetRangeUser(-1.5,1.5)
                    x_reco.Fit("rettax","Q0")
                    x_true.Fit("rettaMC","Q0")
                    #x_reco.Draw("AP")
                    #x_true.Draw("*")
                if iv == 1:
                    xmin = min(zy)-0.05
                    xmax = max(zy)+0.05
                    rettay = TF1("rettay","[0]*x+[1]",xmin,xmax)
                    rettaMC = TF1("rettaMC","[0]*x+[1]",min(ytrack),max(ytrack))
                    y_reco.SetMarkerStyle(11)
                    #y_reco.GetXaxis().SetLimits(-1,1)
                    y_reco.GetYaxis().SetRangeUser(-1.5,1.5)
                    y_reco.Fit("rettay","Q0")
                    y_true.Fit("rettaMC","Q0")
                    #y_reco.Draw("AP")
                    #y_true.Draw("*")
                #ccc[i][iv].Write()

            #PSF (Point Spread Function)
            xchi = rettax.GetChisquare()
            xprob = rettax.GetProb()
            mxrec = rettax.GetParameter(0)
            qxrec = rettax.GetParameter(1)

            qx[0] = qxrec
            mx[0] = mxrec
            chi2x[0] = xchi
            
            ychi = rettay.GetChisquare()
            yprob = rettay.GetProb()
            myrec = rettay.GetParameter(0)
            qyrec = rettay.GetParameter(1)

            qy[0] = qyrec
            my[0] = myrec
            chi2y[0] = ychi

            
            for il in range(totLayers):
                x_res.append(x[il]-(mxrec*zx[il]+qxrec))
                x_residual[il] = x_res[il]

                y_res.append(y[il]-(myrec*zy[il]+qyrec))
                y_residual[il] = y_res[il]

            crFx,crFy,crFz = evaluateRecoCosineDirector(mxrec, myrec)
            #print("Centri", crFx,crFy,crFz)
            cosx[0] = crFx
            cosy[0] = crFy
            cosz[0] = crFz
            cosxMC[0] = cx
            cosyMC[0] = cy
            coszMC[0] = cz
            scalarF = ScalarProduct(cx,cy,cz,crFx,crFy,crFz)   
            thetaF = TMath.ACos(scalarF)*180/TMath.Pi()
            angle_deg[0] = thetaF
            #print("ScalarProduct", scalarF, "PSF",thetaF)
        
        else:
            nlayers[0] = 0
            trk_flag[0] = 0

        OutTree.Fill()

        ####################################################
        #definizione bigstrip:informazioni generali   
        ####################################################

    ###########fine ciclo sugli eventi 
    tfile.Close()
    
    OutTree.Write()
    root_file.Write()
    root_file.Close()














