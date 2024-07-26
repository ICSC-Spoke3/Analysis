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

def move_figure(f, x, y):
    """Move figure's upper left corner to pixel (x, y)"""
    backend = mpl.get_backend()
    if backend == 'TkAgg':
        f.canvas.manager.window.wm_geometry("+%d+%d" % (x, y))
    elif backend == 'WXAgg':
        f.canvas.manager.window.SetPosition((x, y))
    else:
        # This works for QT and GTK
        # You can also use window.setGeometry
        f.canvas.manager.window.move(x, y)

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
    
    #InputFile
    try:
        fin = sys.argv[1]
        #fin     = "/lustrehome/llorusso/Sim_Geant/Builds/JobsOutput/NUSESe-Pow_0.1-5000/rootOutput/NUSESe-Pow_0.1-5000_400000-evt.root"
        #fgeo = "/lustrehome/llorusso/Sim_Geant/Builds/JobsOutput/NUSESe-Pow_0.1-5000/rootOutput/NUSESe-Pow_0.1-5000_400000-evt.txt"
        fgeo = sys.argv[2]
        
    except:
        print("ERROR MESSAGE: \n =====> Insert first argument = ROOTfile.root and second argument = Geofile.txt <=====")
        sys.exit(1)
        
    tfile   = TFile(fin)

    #OutTreeFile
    OutputName = (fin.split("/")[-1]).split(".r")[0] + ".root"
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
    xoff    = 0#rfiber*0.5
    yoff    = 0#rfiber*0.5
    sigma_x = pitch/np.sqrt(12)
    ################################
    #ho 128 strip per l'elettronica di lettura: da 0 a 127
    #quindi l'ultima fibra e' letta dalla 127 fibra
    #le fibre sono 130 da 0 a 129 quindi l'ultima fibra non viene letta
    #Complessivamente non leggiamo meta' della fibra zero, meta' della fibra 128 e la fibra 129
    ################################
    fiberyield = 8*0.92#kev/pe (con cladding)
    PDE = 0.4
    trapeff = 0.054
    PEthresh = 3 #pe per strip

    #definizione canvas
    #ccc = []
    canvas = []
    canvaspsf = []
    for il in range(totLayers):
        canvas.append([])
        canvaspsf.append([])
        canvas[il]=TCanvas("c{0}".format(il),"c{0}".format(il))
        canvas[il].Divide(2,2)
        for iv in range(nview):
            if iv==0:
                v= "x"
            if iv==1:
                v= "y"
            canvaspsf[il].append(TCanvas("c{0}psf{1}".format(il,v),"c{0}psf{1}".format(il,v)))

    cefficiency = TCanvas("cefficiency","cefficiency") 
    cnum        = TCanvas("cnum","cnum") 
    cenergy     = TCanvas("cenergy","cenergy")
    cscalarF    = TCanvas("cscalarF","cscalarF")
    cthetaF     = TCanvas("cthetaF","cthetaF")
    cPSFTH2D    = TCanvas("cPSFTH2D","cPSFTH2D")
    PSFTGraph   = TCanvas("PSFTGraph", "PSFTGraph")
    ccrFz       = TCanvas("ccrFz","ccrFz")
    cQuant      = TCanvas("cQuant", "cQuant")
    #cQuant.SetRightMargin(0.15)
    cQuant.SetLogx()
    cQuant.SetLogy()

    #definizione istogrammi
    '''for i in range(totLayers):
        h.append([])
        for j in range(nview):
            if j == 0:
                v="y"
            if j == 1:
                v="x"
            h[i].append(TH2D("hl{0}v{1}","layer{0} {1}-view;Energy (KeV); Total Photo Electrons",991,90,10000,100,-0.5,500))
    '''
    hdelta = []
    for il in range(totLayers):
        hdelta.append([])
        for iv in range(nview):
            if iv == 0:
                v = "y"
            if iv == 1:
                v = "x"
            hdelta[il].append(TH2D("hdeltal{0}v{1}".format(il,v),"layer{0} {1}-view;Energy (KeV); Delta (mm)".format(il,v),991,90,10000,1000,-50.5,49.5))

    hl0y = TH2D("hl0y","layer0 y-view;Energy (KeV); Total Photo Electrons",9999,90,5000000,100,-0.5,500) #(991,90,10000,100,-0.5,500)
    hl0x = TH2D("hl0x","layer0 x-view;Energy (KeV); Total Photo Electrons",9999,90,5000000,100,-0.5,500) 
    hl1y = TH2D("hl1y","layer1 y-view;Energy (KeV); Total Photo Electrons",9999,90,5000000,100,-0.5,500) 
    hl1x = TH2D("hl1x","layer1 x-view;Energy (KeV); Total Photo Electrons",9999,90,5000000,100,-0.5,500) 
    hl2y = TH2D("hl2y","layer2 y-view;Energy (KeV); Total Photo Electrons",9999,90,5000000,100,-0.5,500) 
    hl2x = TH2D("hl2x","layer2 x-view;Energy (KeV); Total Photo Electrons",9999,90,5000000,100,-0.5,500)

    hl0yclu = TH2D("hl0yclu","layer0 x-view;Energy (KeV); Cluster Size (>3PE)",9999,90,5000000,10,-0.5,10.5) 
    hl0xclu = TH2D("hl0xclu","layer0 y-view;Energy (KeV); Cluster Size (>3PE)",9999,90,5000000,10,-0.5,10.5) 
    hl1yclu = TH2D("hl1yclu","layer1 x-view;Energy (KeV); Cluster Size (>3PE)",9999,90,5000000,10,-0.5,10.5) 
    hl1xclu = TH2D("hl1xclu","layer1 y-view;Energy (KeV); Cluster Size (>3PE)",9999,90,5000000,10,-0.5,10.5) 
    hl2yclu = TH2D("hl2yclu","layer2 x-view;Energy (KeV); Cluster Size (>3PE)",9999,90,5000000,10,-0.5,10.5) 
    hl2xclu = TH2D("hl2xclu","layer2 y-view;Energy (KeV); Cluster Size (>3PE)",9999,90,5000000,10,-0.5,10.5)

    '''hdx0 = TH2D("hdx0","layer0 x-view;Energy (KeV); Delta (mm)",991,90,10000,1000,-50.5,49.5) 
    hdy0 = TH2D("hdy0","layer0 y-view;Energy (KeV); Delta (mm)",991,90,10000,1000,-50.5,49.5) 
    hdx1 = TH2D("hdx1","layer1 x-view;Energy (KeV); Delta (mm)",991,90,10000,1000,-50.5,49.5) 
    hdy1 = TH2D("hdy1","layer1 y-view;Energy (KeV); Delta (mm)",991,90,10000,1000,-50.5,49.5) 
    hdx2 = TH2D("hdx2","layer2 x-view;Energy (KeV); Delta (mm)",991,90,10000,1000,-50.5,49.5) 
    hdy2 = TH2D("hdy2","layer2 y-view;Energy (KeV); Delta (mm)",991,90,10000,1000,-50.5,49.5)
    '''
    ekmin = 0
    ekmax = 5000000.0
    nek = 10000
    hEkDist = TH1D("hEkDist", "Distribution of kinetic energies (keV)", nek, ekmin, ekmax)
    hNTrigViewVsEk = []
    hEffViewVsEk = []
    for il in range (totLayers):
        hNTrigViewVsEk.append([])
        hEffViewVsEk.append([])
        for iv in range(nview):
            title1 = "hNTrigViewVsEk[" + str(il) + "][" + str(iv) + "]"
            title2 = "Number of triggers vs kinetic energy - " + str(il)+str(iv)
            hNTrigViewVsEk[il].append(TH1D(title1, title2, nek, ekmin, ekmax))
            title1 = "hEffViewVsEk[" + str(il) + "][" + str(iv) + "]"
            title2 = "Trigger efficiency vs kinetic energy - "  + str(il)+str(iv)
            hEffViewVsEk[il].append(TH1D(title1, title2, nek, ekmin, ekmax))

    hcrFz    = TH1D("hcrFz", "cz Fibers", 500, -1, 0) #-1<cz<0
    hscalarF = TH1D("hscalarF","C_{reco} #cdot C_{True} with Fibers",500,0,1)
    hthetaF  = TH1D("hthetaF","PSF Fiber",180,0,180)
    hPSFvdEk = TH2D("hPSFvdEk","PSF Fiber vs electron Energy",10000,90,5000000, 180,0,180)
    gra50 = TGraph()
    gra68 = TGraph()
    gra95 = TGraph()

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
    nlayers = ary.array('i',[totLayers])
    chi_flag = ary.array('i',[0])
    x_vec = ary.array('f',range(totLayers))
    y_vec = ary.array('f',range(totLayers))
    zx_vec = ary.array('f',range(totLayers))
    zy_vec = ary.array('f',range(totLayers))
    cosx = ary.array('f',[0.])
    cosy = ary.array('f',[0.])
    cosz = ary.array('f',[0.])
    
    qx = ary.array('f',[0.])
    mx = ary.array('f',[0.])
    chi2x = ary.array('f',[0.])
    qy = ary.array('f',[0.])
    my = ary.array('f',[0.])
    chi2y = ary.array('f',[0.])

    xMC_vec = ary.array('f',range(totLayers))
    yMC_vec = ary.array('f',range(totLayers))
    zxMC_vec = ary.array('f',range(totLayers))
    zyMC_vec = ary.array('f',range(totLayers))
    cosxMC =  ary.array('f',[0.])
    cosyMC =  ary.array('f',[0.])
    coszMC =  ary.array('f',[0.])

    OutTree.Branch("nlayers", nlayers, "nlayers/I")
    OutTree.Branch("chi_flag", chi_flag, "chi_flag/I")
    OutTree.Branch("x_vec", x_vec, "x_vec[nlayers]/F")
    OutTree.Branch("y_vec", y_vec, "y_vec[nlayers]/F")
    OutTree.Branch("zx_vec", zx_vec, "zx_vec[nlayers]/F")
    OutTree.Branch("zy_vec", zy_vec, "zy_vec[nlayers]/F")
    OutTree.Branch("cosx", cosx, "cosx/F")
    OutTree.Branch("cosy", cosy, "cosy/F")
    OutTree.Branch("cosz", cosz, "cosz/F")

    OutTree.Branch("qx", qx, "qx/F")
    OutTree.Branch("mx", mx, "mx/F")
    OutTree.Branch("chi2x", chi2x, "chi2x/F")
    OutTree.Branch("qy", qy, "qy/F")
    OutTree.Branch("my", my, "my/F")
    OutTree.Branch("chi2y", chi2y, "chi2y/F")

    OutTree.Branch("xMC_vec", xMC_vec, "xMC_vec[nlayers]/F")
    OutTree.Branch("yMC_vec", yMC_vec, "yMC_vec[nlayers]/F")
    OutTree.Branch("zxMC_vec", zxMC_vec, "zxMC_vec[nlayers]/F")
    OutTree.Branch("zyMC_vec", zyMC_vec, "zyMC_vec[nlayers]/F")
    OutTree.Branch("cosxMC", cosxMC, "cosxMC/F")
    OutTree.Branch("cosyMC", cosyMC, "cosyMC/F")
    OutTree.Branch("coszMC", coszMC, "coszMC/F")

    tree = tfile.Get("TrackerDigi")
    ptree = tfile.Get("Primary")
    nhit = tree.GetEntries()
    #nhit = 100000

    root_file.cd()
    
    theta=[]
    eK=[]
    
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
        #pixel_vec     = encoder(layer_vec,view_vec,fiber_vec)    

        energy[0] = PrimaryEnergy
        hEkDist.Fill(PrimaryEnergy, 1.0)

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

        #for j in range(len(pixel_vec)):
        for j in range(len(fiber_vec)):
            #print("j",j," pixel", pixel_vec[j])
            mpe = (dE_vec[j]*fiberyield*trapeff*PDE)
            npe1 = ran.Poisson(mpe/2)
            npe2 = ran.Poisson(mpe/2)

            #Lay.append(layer_vec[j])
            #View.append(view_vec[j])

            #Fibra all'estrema sx e estrema dx dei layer lette a metà sull'upper ribble
            #l'ultima non è letta
            #le centrali su ogni strato sono lette da 2 strip ciascuna:
            for rr in range(totLayers*nview):
                #sx
                #if (pixel_vec[j] == rr*nfiber):
                if (fiber_vec[j] == 0):
                    #index1 = rr*stripNo
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
                                        Strip_r[il][iv].append(gxcfib[il][iv][fiber_vec[j]]+xoff)
                                    elif iv == 1:
                                        Strip_r[il][iv].append(gycfib[il][iv][fiber_vec[j]]+yoff)
                                    z[il][iv].append(gzcfib[il][iv][fiber_vec[j]])
                                    eStrips[il][iv].append(dE_vec[j]/2)
                                    PEStrips[il][iv].append(npe1)
                                else:
                                    a=np.where(np.array(Strips[il][iv])==index1)[0][0]
                                    eStrips[il][iv][a]+=(dE_vec[j]/2)
                                    PEStrips[il][iv][a]+=npe1
                                totPE[il][iv] += npe1
                #dx
                #elif (pixel_vec[j] == stripNo + rr*nfiber): 
                elif (fiber_vec[j] == stripNo): 
                    #index1 = (stripNo-1) + rr*stripNo
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
                                        Strip_r[il][iv].append(gxcfib[il][iv][fiber_vec[j]]-xoff)
                                    elif iv == 1:
                                        Strip_r[il][iv].append(gycfib[il][iv][fiber_vec[j]]-yoff)
                                    z[il][iv].append(gzcfib[il][iv][fiber_vec[j]]) 
                                    eStrips[il][iv].append(dE_vec[j]/2)
                                    PEStrips[il][iv].append(npe1)
                                else:
                                    a=np.where(np.array(Strips[il][iv])==index1)[0][0]
                                    eStrips[il][iv][a]+=(dE_vec[j]/2)
                                    PEStrips[il][iv][a]+=npe1    
                                totPE[il][iv] += npe1
                #fibre centrali
                #elif(pixel_vec[j] > rr*nfiber and pixel_vec[j] < stripNo + rr*nfiber ):
                elif(fiber_vec[j] > 0 and fiber_vec[j] < stripNo):
                    #index1 = (pixel_vec[j]-rr*nfiber)+rr*stripNo
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
                                        Strip_r[il][iv].append(gxcfib[il][iv][fiber_vec[j]]+xoff)
                                    elif iv == 1:
                                        Strip_r[il][iv].append(gycfib[il][iv][fiber_vec[j]]+yoff)
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
                                        Strip_r[il][iv].append(gxcfib[il][iv][fiber_vec[j]]-xoff)
                                    elif iv == 1:
                                        Strip_r[il][iv].append(gycfib[il][iv][fiber_vec[j]]-yoff)
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
                #elif (rr > 0 and pixel_vec[j] == rr*nfiber -1): 
                elif (fiber_vec[j] == nfiber -1): 
                    continue
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

                '''for kk in range(1,len(z[il][iv])):
                    #print(z[il][iv])
                    if z[il][iv][kk] != z[il][iv][kk-1]:
                        String_zs[il][iv]=(z[il][iv][kk]+z[il][iv][kk-1])/2
                    else:
                        String_zs[il][iv]=z[il][iv][kk]
                    #print("String_zs[0]", String_zs[il][0])
                    #print("String_zs[1]", String_zs[il][1])'''
                for istrip in range(len(Strips[il][iv])-1):
                    if (PEStrips[il][iv][istrip] >= PEthresh and PEStrips[il][iv][istrip+1] >= PEthresh):
                        vtrig[il][iv] = 1

                if vtrig[il][iv]>0:
                    hNTrigViewVsEk[il][iv].Fill(PrimaryEnergy, 1.0)

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
                    x1 = istrip_vec[il][iv][-1]
                    num = sum(np.multiply(Strip_r[il][iv][x0:x1],PEStrips[il][iv][x0:x1]))
                    den = sum(PEStrips[il][iv][x0:x1])

                    if den >= 4:  #almeno 4 fotoelettroni a strip

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
                            dzx.append(1/np.sqrt(12))
                            xtrack.append(primary_track[il][iv])
                            
                            x_vec[il] = bc[il][iv]                            
                            zx_vec[il] =  String_zs[il][iv]   
                            xMC_vec[il] = primary_track[il][iv]
                            zxMC_vec[il] = String_zs[il][iv]

                        elif iv==1:
                            y.append(bc[il][iv])
                            zy.append(String_zs[il][iv])
                            dy.append(debc[il][iv])
                            dzy.append(1/np.sqrt(12))
                            ytrack.append(primary_track[il][iv])

                            y_vec[il] = bc[il][iv]
                            zy_vec[il] = String_zs[il][iv]
                            yMC_vec[il] = primary_track[il][iv]
                            zyMC_vec[il] = String_zs[il][iv]

                        delta = bc[il][iv] - y0MC
                        hdelta[il][iv].Fill(PrimaryEnergy,delta)


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

        if (np.shape(x)[0] == totLayers and np.shape(y)[0] == totLayers):        
            x_reco = TGraphErrors(totLayers,zx,x,dzx,dx)
            x_reco.SetTitle("Hor")
            x_true = TGraph(3,zx,xtrack)
            #x_true = TGraph(len(xTMC),zTMC,xTMC) #posizione MC layer per layer
            y_reco = TGraphErrors(totLayers,zy,y,dzy,dy)
            y_reco.SetTitle("Ver")
            y_true = TGraph(3,zy,ytrack)  
            #y_true = TGraph(len(xTMC),zTMC,yTMC)  

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
                    x_reco.Draw("AP")
                    x_true.Draw("*")
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
                    y_reco.Draw("AP")
                    y_true.Draw("*")
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
            #print("ScalarProduct", scalarF, "PSF",thetaF)
                
            if xchi < 2 and ychi < 2:

                chi_flag[0] = 1

                hcrFz.Fill(crFz)
                hscalarF.Fill(scalarF)
                hthetaF.Fill(thetaF)
                hPSFvdEk.Fill(PrimaryEnergy, thetaF)

                theta.append(thetaF)
                eK.append(PrimaryEnergy)

            #reco vs MC track (con due spezzate)
            else:
                #ccc[i].append(TCanvas("ccc{0}_{1}".format(iv,i),"ccc{0}_{1}".format(iv,i)))
                #ccc[i][iv].cd()
                x_reco2 = TGraphErrors(totLayers,zx,x,dzx,dx)
                x_reco2.SetTitle("Hor")
                #x_true = TGraph(3,zx,xtrack)
                y_reco2 = TGraphErrors(totLayers,zy,y,dzy,dy)
                y_reco2.SetTitle("Ver")
                #y_true = TGraph(3,zy,ytrack)  
               
                r1 = TF1("r1","[0]*x+[1]",zx[0],zx[1])
                r2 = TF1("r2","[0]*x+[1]",zx[1],zx[2])
               
                x_reco2.SetMarkerStyle(8)
                x_reco2.GetYaxis().SetRangeUser(-1.5,1.5)
                x_reco2.Fit("r1","Q0")
                x_reco2.Fit("r2","Q0")
                x_reco2.Draw("AP")
               
                
        

        OutTree.Fill()

        ####################################################
        #definizione bigstrip:informazioni generali   
        ####################################################

    ###########fine ciclo sugli eventi 
    tfile.Close()

    '''nPSF=len(theta)
    PSF = np.array(theta).astype(float)
    eK_prim = np.array(eK).astype(float)
    PSFvsE=TGraph(nPSF, eK_prim,PSF)
    PSFTGraph.cd()
    PSFvsE.Draw("AP*")'''

    #aggiungo quantili
    sorte = sorted(eK)
    net1 = len(eK)
    Et = ary.array('f',range(net1+1))

    Et[0] = sorte[0]
    j = 1
    for i in range(1,net1):
        if(sorte[i] > sorte[i-1]):
            Et[j] = sorte[i]
            j += 1

    net = j-1

    h1Ngen = TH1F("h1Ngen","; MC Energy (MeV); Number of Generated Events", net, Et)
    for i in range(net1):
        h1Ngen.Fill(eK[i], nhit)

    nq = 3
    xq = ary.array("d",range(nq))
    yq = ary.array("d",range(nq))
    xq[0] = 0.50
    xq[1] = 0.68
    xq[2] = 0.95

    g = 0
    for k in range(len(eK)):
        xg = h1Ngen.GetBinContent(k+1)
        #x  = hPSFvdEk.GetXaxis().GetBinCenter(k+1)
        x  = hPSFvdEk.GetXaxis().GetBinLowEdge(k+1)
        if(xg>0.):
            htmp_py = hPSFvdEk.ProjectionY("htmp_py", k, k+1)
            if(htmp_py.Integral()>0.):
                htmp_py.GetQuantiles(nq, yq, xq)
                gra50.SetPoint(g, x, yq[0])
                gra68.SetPoint(g, x, yq[1])
                gra95.SetPoint(g, x, yq[2])

                #print ("quantiles: ", k, "LowEdge", x, "quant_vec", yq)
                g += 1

    cQuant.cd()
    cQuant.SetGrid()
    leg = TLegend(0.18,.66,.38,.86)
    leg.SetBorderSize(1)
    hPSFvdEk.GetXaxis().SetTitle("Energy [keV]")
    hPSFvdEk.GetYaxis().SetTitle("#theta (#circ)")
    hPSFvdEk.Draw("colz")
    gra50.SetLineWidth(3)
    gra68.SetLineWidth(3)
    gra95.SetLineWidth(3)
    gra50.SetLineStyle(8)
    gra68.SetLineStyle(8)
    gra95.SetLineStyle(8)
    gra50.SetLineColor(616+2)
    gra68.SetLineColor(800+7)
    gra95.SetLineColor(820-9)
    gra50.Draw("lsame")
    gra68.Draw("lsame")
    gra95.Draw("lsame")
    #comportamento che tiene conto dello scattering multiplo coulombiano
    ll = 0.8 #mm
    kk = 400 #mm
    ww = 0.1 #mm
    pp = 25 #mm
    l_fisica = TF1("l_fisica", "pow(pow([0],2)+pow([1]*1/x,2),0.5)", 100,5000000) #radice(dtheta^2_msc+dtheta^2_risoluzione)
    l_fisica.SetParameter(0, ww/pp)
    l_fisica.SetParameter(1, 13.6e3*pow(ll/ww,0.5))
    l_fisica.SetLineColor(2)
    l_fisica.Draw("same")
    leg.AddEntry(gra50,"50%","l")
    leg.AddEntry(gra68,"68%","l")
    leg.AddEntry(gra95,"95%","l")
    leg.AddEntry(l_fisica,"teorica","l")
    leg.Draw()
    cQuant.Update() 

    '''
    canvas[0].cd(1)
    OutTree.Draw("TpeStripHor0:energy>>hl0y","","colz")
    canvas[0].cd(2)
    OutTree.Draw("TpeStripVer0:energy>>hl0x","","colz")
    canvas[0].cd(3)
    OutTree.Draw("CsizHor0:energy>>hl0yclu","PEstripHor0 >= 3 || CsizHor0 == 0","colz")
    canvas[0].cd(4)
    OutTree.Draw("CsizVer0:energy>>hl0xclu","PEstripVer0 >= 3 || CsizVer0 == 0","colz")

    canvas[1].cd(1)
    OutTree.Draw("TpeStripHor1:energy>>h[1][0]","","colz")
    canvas[1].cd(2)
    OutTree.Draw("TpeStripVer1:energy>>h[1][1]","","colz")
    canvas[1].cd(3)
    OutTree.Draw("CsizHor1:energy>>hl1yclu","PEstripHor1 >= 3 || CsizHor1 == 0","colz")
    canvas[1].cd(4)
    OutTree.Draw("CsizVer1:energy>>hl1xclu","PEstripVer1 >= 3 || CsizVer1 == 0 ","colz")

    canvas[2].cd(1)
    OutTree.Draw("TpeStripHor2:energy>>h[2][0]","","colz")
    canvas[2].cd(2)
    OutTree.Draw("TpeStripVer2:energy>>h[2][1]","","colz")
    canvas[2].cd(3)
    OutTree.Draw("CsizHor2:energy>>hl2yclu","PEstripHor2 >= 3 || CsizHor2 == 0 ","colz")
    canvas[2].cd(4)
    OutTree.Draw("CsizVer2:energy>>hl2xclu","PEstripVer2 >= 3 || CsizVer2 == 0 ","colz")'''

    cenergy.cd()
    hEkDist.Draw()
    color=[[416+1, 4], [632-6, 800-2], [432-5, 28]]
    label=[]

    hEffViewVsEk[0][0].SetTitle("Trigger efficiency vs Ek")
    hEffViewVsEk[0][0].SetStats(0)
    legend = TLegend(0.66,.32,.89,.59)
    legend.SetBorderSize(1)

    for il in range(totLayers):
        label.append([])

        for iv in range(nview):
            if iv == 0:
                v = "Hor"
                coo = "x"
            else:
                v= "Ver"
                coo="y"
                
            sCsiz       = "Csiz{0}{1}".format(v,il)
            sTpeStrip = "TpeStrip{0}{1}".format(v,il)
            sPEstrip    = "PEstrip{0}{1}".format(v,il)
            
            canvas[il].cd(iv+1)
            title = sTpeStrip + ":energy>>hl"+str(il)+coo
            cut = ""
            OutTree.Draw(title, cut, "colz")
            
            canvas[il].cd(iv+3)
            title = sCsiz + ":energy>>hl"+str(il)+coo +"clu"
            cut = sPEstrip + ">=3 || " + sCsiz + " == 0"
            OutTree.Draw(title, cut, "colz")

            label[il].append("Layer {0} View {1}".format(il,iv))
            cnum.cd()
            if il == 0 and iv == 0:
                hNTrigViewVsEk[il][iv].Draw("")
            else: hNTrigViewVsEk[il][iv].Draw("same")

            cefficiency.cd()
            cefficiency.SetGrid()

            hEffViewVsEk[il][iv].Divide(hNTrigViewVsEk[il][iv],hEkDist,nek,nek,"B")    

            hEffViewVsEk[il][iv].GetYaxis().SetTitle("Trigger efficiency")
            hEffViewVsEk[il][iv].GetYaxis().SetRangeUser(0, 1.05)
            hEffViewVsEk[il][iv].GetYaxis().SetTitleFont(22)   
            hEffViewVsEk[il][iv].GetXaxis().SetTitle("Energy [keV]")
            hEffViewVsEk[il][iv].GetXaxis().SetTitleFont(22)       

            hEffViewVsEk[il][iv].SetLineColor(color[il][iv])
            hEffViewVsEk[il][iv].SetLineWidth(1)
            hEffViewVsEk[il][iv].SetMarkerColor(color[il][iv])
            hEffViewVsEk[il][iv].SetMarkerStyle(20)

            if il == 0 and iv == 0:
                hEffViewVsEk[il][iv].Draw("")
            else: 
                hEffViewVsEk[il][iv].Draw("same")
            cefficiency.Update()

            legend.AddEntry(hEffViewVsEk[il][iv],label[il][iv],"lp")
            cefficiency.Update()
    legend.Draw();   

    ccrFz.cd()
    hcrFz.Draw()
    cscalarF.cd()
    hscalarF.Draw()
    cthetaF.cd()
    hthetaF.GetXaxis().SetTitle("#theta (#circ)")
    hthetaF.Draw()
    cPSFTH2D.cd()
    hPSFvdEk.Draw("colz")

    for il in range(totLayers):
        for iv in range(nview):
            canvaspsf[il][iv].cd()
            hdelta[il][iv].Draw("colz") 
            canvaspsf[il][iv].Write()

    cenergy.Write()
    #cnum.Write()
    cefficiency.Write()
    canvas[0].Write()
    canvas[1].Write()
    canvas[2].Write()
    cscalarF.Write()
    cthetaF.Write()
    cPSFTH2D.Write()
    #PSFTGraph.Write()
    ccrFz.Write()
    cQuant.Write()
    
    OutTree.Write()
    root_file.Write()
    root_file.Close()

