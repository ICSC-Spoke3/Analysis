from turtle import shape
from ROOT import TFile, TTree, TChain, TH1D, TH2D, TMath, TCanvas, TH1F,TRandom3,TGraphErrors, TGraph,TF1
import os
import sys
import time
import array as ary
import numpy as np
import glob
import math
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
        #print("nobiiii")
        nclu=0
        istrip0 = [-1]
        nstrip0 = [0]
        charge0 = [0]    
        #print (istrip0)

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

def cstripCalc(Rfib,c0,offset,pitch,iFib,cfib,StripNo):
    strip_centres = []
    strip_index = []
    
    xminfib = cfib - Rfib
    xmaxfib = cfib + Rfib

    #print("FIBRA:", iFib)
    #print ("centro fibra:",cfib)
    #print ("xminfib",xminfib)
    #print ("xmaxfib",xmaxfib)
    
    NStripPerFiber = 2*Rfib/pitch
    #print("Numero di strip per fibra:",NStripPerFiber)
    
    stripFirstfib = (2*Rfib - offset)/pitch
    #print("Strip nella prima fibra:",stripFirstfib)
    
    if iFib == 0:
        StripIndex = 0
    else:
        StripIndex = math.ceil(NStripPerFiber*(iFib/2-1)+ stripFirstfib)
    #print("StripIndex",StripIndex)

    cstrip = c0 + pitch*(0.5+StripIndex)    
    #print("CENTRO STRIP:",cstrip)

    xminstrip = cstrip - pitch/2
    xmaxstrip = cstrip + pitch/2 
    #print ("xminstrip",xminstrip)
    #print ("xmaxstrip",xmaxstrip)
    
    strip_centres.append(cstrip)
    strip_index.append(StripIndex)
    
    xminstrip1 = xminstrip
    xmaxstrip1 = xmaxstrip
    StripIndex1 = StripIndex
    cstrip1 = cstrip
    while round(xminstrip1,2) > round(xminfib,2) and  xminstrip1 <= xmaxfib and xmaxstrip1 >= xminfib and StripIndex > 0 :
        cstrip1 = cstrip1 - pitch
        xminstrip1 = cstrip1 - pitch/2
        xmaxstrip1 = cstrip1 + pitch/2
        StripIndex1 -= 1
        
        strip_centres.append(cstrip1)
        strip_index.append(StripIndex1)

        #print("cstrip1",cstrip1)
        #print("xminstrip1",xminstrip1)
        # #print("strip index:",StripIndex1)
        
    xminstrip2 = xminstrip
    xmaxstrip2 = xmaxstrip
    StripIndex2 = StripIndex 
    cstrip1 = cstrip
    while round(xmaxstrip2,2) < round(xmaxfib,2) and xminstrip2 <= xmaxfib and xmaxstrip2 >= xminfib and StripIndex < StripNo:
        cstrip1 = cstrip1 + pitch
        xminstrip2 = cstrip1 - pitch/2
        xmaxstrip2 = cstrip1 + pitch/2 
        StripIndex2 += 1
        
        strip_centres.append(cstrip1)
        strip_index.append(StripIndex2)

        #print("cstrip1",cstrip1)
        #print("xmaxstrip2",xmaxstrip2)
        #print("strip index:",StripIndex2)

    return strip_centres, strip_index

def FracArea(Rfib,pitch,cfib,cstrip):
    #print("(#######################################################)")
    AFib = math.pi*pow(Rfib,2)
        
    '''print("cfib",cfib)
    print("cstrip",cstrip)
    print("pitch",pitch)
    print("Rfib",Rfib)'''
    
    cfib = cfib
    xminfib = cfib - Rfib
    xmaxfib = cfib + Rfib
    xminstrip = cstrip - pitch/2
    xmaxstrip = cstrip + pitch/2 
    
    '''print("xminfib-xmaxfib",xminfib,xmaxfib)
    print("xminstrip-xmaxstrip",xminstrip,xmaxstrip)
    '''
    if xminfib < xminstrip:
        #print("xminfib < xminstrip",xminfib,"<",xminstrip)
        xmin = xminstrip
    elif xminfib >= xminstrip:
        #print("xminfib >= xminstrip",xminfib,">=",xminstrip)
        xmin = xminfib

    xmax = min(xmaxfib,xmaxstrip)
    #print ("xmin-xmax",xmin,xmax)

    costh1 = round((xmin - cfib)/Rfib,2)#cosalpha
    costh2 = round((xmax - cfib)/Rfib,2)#cosbetha
    '''print("xmin - cfib",xmin - cfib)
    print("xmax - cfib",xmax - cfib)
    print("costh1",costh1)
    print("costh2",costh2)'''


    th1 = math.acos(costh1)#alpha
    th2 = math.acos(costh2)#betha

    sin1 = math.sin(2*th1)
    sin2 = math.sin(2*th2)

    AStrip = pow(Rfib,2)*((th1-th2) + sin2/2 - sin1/2)
    fracA=AStrip/AFib

    '''print("fracA=", fracA)
    print("(#######################################################)")
   '''
    return fracA

def CalcResiduals (n,x,z,m,q):
    r = []
    for i in range(n):
        x_fit = m*z[i]+q
        r.append(x[i]-x_fit)
    return r

if __name__ == "__main__":

    #InputFile
    try:
        fin = sys.argv[1]
        #fin     = "/home/nadia_work/Scrivania/NUSESe-Pow_0.1-5000_1000000-evt.root"
        #fgeo = "/home/nadia_work/Scrivania/NUSESe-Pow_0.1-5000_1000000-evt.txt"
        fgeo = sys.argv[2]
        
    except:
        print("ERROR MESSAGE: \n =====> Insert first argument = ROOTfile.root - second argument = Geofile.txt - third argument = pitch (mm)<=====")
        sys.exit(1)
        
    tfile   = TFile(fin)

    #OutTreeFile
    pitch   = sys.argv[3] #0.25 #mm
    OutputName = "Results/Tree/Tree_"+(fin.split("/")[-1]).split(".r")[0] + "_pitch"+str(pitch)+".root"
    OutputFile = OutputName
    print("OUTPUT:", OutputFile)
    root_file = TFile(OutputFile, "RECREATE")
    
    OutTree           = TTree("ClusterTree", "ClusterTree")
    BigClusterTree    = TTree("BigClusterTree", "BigClusterTree")

    Track_info,Calo_info,gxcfib,gycfib,gzcfib,gxcpix,gycpix,gzcpix = read_GeoFile(fgeo)
    Layers        = Track_info[0]
    Views         = Track_info[1]
    Fibers        = Track_info[2]
    FibLength     = Track_info[3]
    FibRadius     = Track_info[4]
    TrackerLength = Track_info[5]
    TrackerWidth  = Track_info[6]
    
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
    print("FIBER RADIUS:", rfiber)
    xoff    = rfiber#rfiber*0.5
    yoff    = rfiber#rfiber*0.5
    cx0strip = gxcfib[0][0][0] - rfiber + xoff
    cy0strip = gycfib[0][1][0] - rfiber + yoff
    pitch   = float(pitch) #0.25 #mm
    print("[INFO]Strip pitch:",pitch)
    stripNo = int((nfiber*rfiber+rfiber)/pitch)#128
    sigma_x = pitch/2#sqrt(12)
    ################################
    #ho 128 strip per l'elettronica di lettura: da 0 a 127
    #quindi l'ultima fibra e' letta dalla 127 fibra
    #le fibre sono 130 da 0 a 129 quindi l'ultima fibra non viene letta
    #Complessivamente non leggiamo meta' della fibra zero, meta' della fibra 128 e la fibra 129
    ################################
    fiberyield = 8*0.92#kev/pe
    PDE = 0.4
    trapeff = 0.054
    PEthresh = 3 #pe per strip

    # Define output tree
    nChan = 400000 
    eventID = ary.array('i',[0])
    Tch = ary.array('i',[0])
    energy  = ary.array('f',[0])
    Ipix = ary.array('i',range(nChan))
    Npe = ary.array('f',range(nChan))

    OutTree.Branch("eventID", eventID, "eventID/I")
    OutTree.Branch("Tch", Tch, "Tch/I")
    OutTree.Branch("energy", energy, "energy_keV/F")
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
                pos = "x"
            else:
                v="Ver"
                pos = "y"
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

    ####informazioni sul tracciamento
        
    nlayers  = ary.array('i',[0])
    nlayersx  = ary.array('i',[0])
    nlayersy  = ary.array('i',[0])
    trk_flag = ary.array('i',[0])
    Ilayerx = ary.array('i',range(nChan))
    Ilayery = ary.array('i',range(nChan))
    x_vec = ary.array('f',range(totLayers))
    y_vec = ary.array('f',range(totLayers))
    zx_vec = ary.array('f',range(totLayers))
    zy_vec = ary.array('f',range(totLayers))
    dx_vec = ary.array('f',range(totLayers))
    dy_vec = ary.array('f',range(totLayers))
    dzx_vec = ary.array('f',range(totLayers))
    dzy_vec = ary.array('f',range(totLayers))
    
    xMC_vec = ary.array('f',range(totLayers))
    yMC_vec = ary.array('f',range(totLayers))
    zxMC_vec = ary.array('f',range(totLayers))
    zyMC_vec = ary.array('f',range(totLayers))

    mx = ary.array('f',[0.])
    qx = ary.array('f',[0.])
    chi2x = ary.array('f',[0.])
    my = ary.array('f',[0.])
    qy = ary.array('f',[0.])
    chi2y = ary.array('f',[0.])

    cosx= ary.array('f',[0.])
    cosy= ary.array('f',[0.])
    cosz= ary.array('f',[0.])

    cosxMC= ary.array('f',[0.])
    cosyMC= ary.array('f',[0.])
    coszMC= ary.array('f',[0.])
    
    angle_deg= ary.array('f',[0.])
   
    OutTree.Branch("nLayers", nlayers, "nlayers/I")
    OutTree.Branch("nLayersx", nlayersx, "nlayersx/I")
    OutTree.Branch("nLayersy", nlayersy, "nlayersy/I")
    OutTree.Branch("trk_flag", trk_flag, "trk_flag/I")
    OutTree.Branch("Ilayerx", Ilayerx, "Ilayerx[nlayersx]/I")
    OutTree.Branch("Ilayery", Ilayery, "Ilayery[nlayersy]/I")
    OutTree.Branch("x_vec", x_vec, "x_vec[nlayersx]/F")
    OutTree.Branch("zx_vec", zx_vec, "zx_vec[nlayersx]/F")
    OutTree.Branch("y_vec", y_vec, "y_vec[nlayersy]/F")
    OutTree.Branch("zy_vec", zy_vec, "zy_vec[nlayersy]/F")
    OutTree.Branch("dx_vec", dx_vec, "dx_vec[nlayersx]/F")
    OutTree.Branch("dzx_vec", dzx_vec, "dzx_vec[nlayersx]/F")
    OutTree.Branch("dy_vec", dy_vec, "dy_vec[nlayersy]/F")
    OutTree.Branch("dzy_vec", dzy_vec, "dzy_vec[nlayersy]/F")
    
    OutTree.Branch("xMC_vec", xMC_vec, "xMC_vec[nlayersx]/F")
    OutTree.Branch("zxMC_vec", zxMC_vec, "zxMC_vec[nlayersx]/F")
    OutTree.Branch("yMC_vec", yMC_vec, "yMC_vec[nlayersy]/F")
    OutTree.Branch("zyMC_vec", zyMC_vec, "zyMC_vec[nlayersy]/F")

    OutTree.Branch("intercept_x", qx, "qx[trk_flag]/F")
    OutTree.Branch("slope_x", mx, "mx[trk_flag]/F")
    OutTree.Branch("Chi2_x", chi2x, "chi2x[trk_flag]/F")

    OutTree.Branch("intercept_y", qy, "qy[trk_flag]/F")
    OutTree.Branch("slope_y", my, "my[trk_flag]/F")
    OutTree.Branch("Chi2_y", chi2y, "chi2y[trk_flag]/F")

    OutTree.Branch("Cosx", cosx, "cosx[trk_flag]/F")
    OutTree.Branch("Cosy", cosy, "cosy[trk_flag]/F")
    OutTree.Branch("Cosz", cosz, "cosz[trk_flag]/F")

    OutTree.Branch("Cosx_MC", cosxMC, "cosxMC[trk_flag]/F")
    OutTree.Branch("Cosy_MC", cosyMC, "cosyMC[trk_flag]/F")
    OutTree.Branch("Cosz_MC", coszMC, "coszMC[trk_flag]/F")

    OutTree.Branch("Theta_deg", angle_deg, "angle_deg[trk_flag]/F")
      
    for iv in range(nview):
        if iv == 0: 
            v = "Hor"
            pos = "x"
        else:
            v="Ver"
            pos = "y"

        sResiduals    = "Residuals{0}".format(pos)
        sResiduals2    = "Residuals{0}[nlayers]/F".format(pos)

        exec("%s = ary.array('f',range(totLayers))" % sResiduals)

        OutTree.Branch(sResiduals, eval(sResiduals),sResiduals2)

    tree = tfile.Get("TrackerDigi")
    ptree = tfile.Get("Primary")
    nhit = tree.GetEntries()
    #nhit = 10000
    
    root_file.cd()
    e =[]
    theta =[]
    tetha =[]

    for i in range(nhit):
        eventID[0] = i
        if i%10000==0:
            print("Event: ", eventID[0])
        #print("Event: ", eventID[0])
        #event[0]   = i
        #--- Primary Track ---#
        ptree.GetEntry(i)
        xmc = ptree.PrimaryParticlePositionX
        ymc = ptree.PrimaryParticlePositionY
        zmc = ptree.PrimaryParticlePositionZ
        cx = ptree.PrimaryParticleDirectionX
        cy = ptree.PrimaryParticleDirectionY
        cz = ptree.PrimaryParticleDirectionZ
        #lettura DigiTrackerTree 
        tree.GetEntry(i) 
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
        Layx  = []
        Layy  = []
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
        
        '''print("Layers:",list(layer_vec))
        print("View:",list(view_vec))
        print("Fibers:",list(fiber_vec))'''
        
        for j in range(len(fiber_vec)):
            
            mpe_fiber = (dE_vec[j]*fiberyield*trapeff*PDE)#numero medio di PE nella fibra

            if view_vec[j] == 0:
                cfib = gxcfib[layer_vec[j]][view_vec[j]][fiber_vec[j]]
                strip_centres, strip_index = cstripCalc(rfiber,cx0strip,xoff,pitch,fiber_vec[j],gxcfib[layer_vec[j]][view_vec[j]][fiber_vec[j]],stripNo)
            if view_vec[j] == 1:
                cfib = gycfib[layer_vec[j]][view_vec[j]][fiber_vec[j]]
                strip_centres, strip_index = cstripCalc(rfiber,cy0strip,yoff,pitch,fiber_vec[j],gycfib[layer_vec[j]][view_vec[j]][fiber_vec[j]],stripNo)
            Strip_vec.extend(strip_index)
            

            for istrip in range(len(strip_index)):
                Afract = round(FracArea(rfiber,pitch,cfib,strip_centres[istrip]),1)
                mpe = dE_vec[j]*Afract*fiberyield*trapeff*PDE
                npe = ran.Poisson(mpe)#numero medio di PE nella strip
                PEstrip_vec.append(npe)###########################
            
            #print("#########fiber",j,":",fiber_vec[j], "Strip_vec",Strip_vec)

            '''for kk in range(len(strip_index)):
                Lay.append(layer_vec[j])
                View.append(view_vec[j])'''
            il = layer_vec[j]
            iv = view_vec[j]

            if (iv == 0): 
                strip_centres, strip_index = cstripCalc(rfiber,cx0strip,xoff,pitch,fiber_vec[j],gxcfib[il][iv][fiber_vec[j]],stripNo)
                cfib=gxcfib[il][iv][fiber_vec[j]]
            if (iv == 1): 
                strip_centres, strip_index = cstripCalc(rfiber,cy0strip,yoff,pitch,fiber_vec[j],gycfib[il][iv][fiber_vec[j]],stripNo)
                cfib=gycfib[il][iv][fiber_vec[j]]
            for istrip in range(len(strip_index)):
                Afract = round(FracArea(rfiber,pitch,cfib,strip_centres[istrip]),1)
                mpe = dE_vec[j]*Afract*fiberyield*trapeff*PDE
                npe = ran.Poisson(mpe)
                if strip_index[istrip] not in Strips[il][iv]:
                    Strips[il][iv].append(strip_index[istrip])
                    Strip_r[il][iv].append(strip_centres[istrip])
                    z[il][iv].append(gzcfib[il][iv][fiber_vec[j]])
                    eStrips[il][iv].append(dE_vec[j]*Afract)
                    PEStrips[il][iv].append(npe) 
                else:
                    a=np.where(np.array(Strips[il][iv])==strip_index[istrip])[0][0]
                    eStrips[il][iv][a]+=(dE_vec[j]*Afract)
                    PEStrips[il][iv][a]+=npe
                totPE[il][iv] += npe
        #print(">>>>>>>>>>>>>>>>>>>>eStrips", eStrips)
        #print("dE_vec", list(dE_vec))
        ordered_index_list =[]# np.argsort(Strips)
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

                '''print("STRIP:",Strips[il][iv]) 
                print(ordered_strip[il][iv])
                print("CENTRE:",Strip_r[il][iv])
                print(ordered_centre[il][iv])
                print("ENERGY:",eStrips[il][iv])
                print(ordered_energy[il][iv])
                print("PE",PEStrips[il][iv])
                print(ordered_PE[il][iv])'''
        
        Strips   = ordered_strip
        Strip_r  = ordered_centre
        eStrips  = ordered_energy
        PEStrips = ordered_PE 
        #####fine ciclo sulle hit dell'evento i    
        '''print("strip_r0x",Strip_r[0][0])
        print("strip_r0y",Strip_r[0][1])
        print("strip_r1x",Strip_r[1][0])
        print("strip_r1y",Strip_r[1][1])'''
        #print("fine ciclo sulle hit dell'evento",i)
       
        vtrig= np.zeros((totLayers,nview))
        Tch[0] = (len(PEstrip_vec))
        
        ####CLUSTERING
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
                #print("istrip0",istrip0)
                Nclu[il][iv] = nclu
                for kk in range(len(istrip0)):
                    Fstrip[il][iv].append(istrip0[kk])
                    Clusiz[il][iv].append(nstrip0[kk])
                    Clunpe[il][iv].append(charge0[kk])
                #print("FSTRIP:", Fstrip)
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
                    #print("éééé",Fstrip[il][iv][icl])
                #print("Fstr",eval(sFstr)[0])
                if Nclu[il][iv]==0: 
                    exec("%s[0] = 0" % (sCsiz))
                    exec("%s[0] = 0" % (sTnpe))
                    exec("%s[0] = -1" % (sFstr))
                    exec("%s[0] = 1" % (sdim))

        Tclu_tot[0] = Nclu_tot
        dim_tot[0] = Nclu_tot
        for ihit in range (Tch[0]):
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
        view_flag = []
        primary_track = []
        xtrack = []
        ytrack = []
        for il in range(totLayers):
            istrip_vec.append([])
            bc.append([])
            debc.append([])
            view_flag.append([])
            primary_track.append([])
            for iv in range(nview):
                #bc[il].append([])
                #debc[il].append([])
                view_flag[il].append([])
                primary_track[il].append([])
                #print("view",iv)
                istrip_vec[il].append([])
                if iv == 0: 
                    v = "Hor"
                    coo ="x"
                else:
                    v="Ver"
                    coo ="y"
                sTclu      = "Tclu{0}{1}".format(v,il)
                sTnpe      = "Tnpe{0}{1}".format(v,il)
                sTpe       = "Tpe{0}{1}".format(v,il)
                sstrip     = "strip{0}{1}".format(v,il)
                sFstr      = "Fstr{0}{1}".format(v,il)
                sCsiz      = "Csiz{0}{1}".format(v,il)
                sPEstrip   = "PEstrip{0}{1}".format(v,il)

                if eval(sTclu)[0] > 0 :

                    #print(i,"Tclu{0}{1} > 0 :".format(v,il),eval(sTclu)[0])

                    if iv == 0: 
                        Layx.append(il)
                    elif iv == 1:
                        Layy.append(il)
                                        
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

                    if den >= 3:
                        #print(i,"den>3:",den)
                        de1 = np.array(PEStrips[il][iv][x0:x1])*sigma_x/den
                        sumde2 = np.array(Strip_r[il][iv][x0:x1])*den-np.array(Strip_r[il][iv][x0:x1])*np.array(PEStrips[il][iv][x0:x1])
                        de2 = np.array(sumde2) *pow(np.array(PEStrips[il][iv][x0:x1]),1/2)/pow(den,2)
                        bc[il].append(num/den)
                        debc[il].append(pow(sum(pow(de1,2))+sum(pow(de2,2)),0.5))
                        
                        if iv==0:
                            primary_track[il][iv] = (evaluateTrackPoint(xmc,zmc,cx,cz,String_zs[il][iv],flag))
                            
                        elif iv==1:
                            primary_track[il][iv] = (evaluateTrackPoint(ymc,zmc,cy,cz,String_zs[il][iv],flag))
        
                    
                    #print("iv", iv, "primary_track", len(primary_track))

        #print("String", String_zs)             
        #print("EVENT",i)
        #print("bc",bc)
        indexl = 0
        indexv = []
        hitview = []
        for il in range(totLayers):
            #print("il",il, "indexl",indexl, "LENBC",len(bc[il]))
            #print("LEN",len(bc[il]))
            if len(bc[il])>0:
                indexl+=1
                indexv.append(il)
                hitview.append(len(bc[il]))
        #print(indexl,indexv)
        for il in range(indexl):
            #print("IL", il)
            for iv in range(hitview[il]):
                #print("IV", iv)
                if hitview[il]==1 and primary_track[il][1]!=[]:
                    #print("1")
                    #print("vistay")
                    y_vec[il] = bc[indexv[il]][iv]
                    dy_vec[il] = debc[indexv[il]][iv]
                    zy_vec[il] = String_zs[indexv[il]][1]
                    dzy_vec[il] = 2*rfiber
                    yMC_vec[il] = primary_track[indexv[il]][1]
                    zyMC_vec[il] = String_zs[indexv[il]][1]

                    y.append(bc[indexv[il]][iv])
                    dy.append(debc[indexv[il]][iv])
                    zy.append(String_zs[indexv[il]][1])
                    dzy.append(2*rfiber)
                    ytrack.append(primary_track[indexv[il]][1])

                elif hitview[il]==2:
                    if iv == 0:
                        #print("vistax")
                        zx_vec[il] = String_zs[indexv[il]][iv]
                        dzx_vec[il] = 2*rfiber
                        x_vec[il] = bc[indexv[il]][iv]
                        dx_vec[il] = debc[indexv[il]][iv]
                        xMC_vec[il] = primary_track[indexv[il]][iv]
                        zxMC_vec[il] = String_zs[indexv[il]][iv]

                        x.append(bc[indexv[il]][iv])
                        dx.append(debc[indexv[il]][iv])
                        zx.append(String_zs[indexv[il]][iv])
                        dzx.append(2*rfiber)                                                      
                        xtrack.append(primary_track[indexv[il]][iv])

                    elif iv == 1:
                        #print("vistay")
                        y_vec[il]=bc[indexv[il]][iv]
                        dy_vec[il] = debc[indexv[il]][iv]
                        zy_vec[il] = String_zs[indexv[il]][iv]
                        dzy_vec[il] = 2*rfiber
                        yMC_vec[il] = primary_track[indexv[il]][iv]
                        zyMC_vec[il] = String_zs[indexv[il]][iv]

                        y.append(bc[indexv[il]][iv])
                        dy.append(debc[indexv[il]][iv])
                        zy.append(String_zs[indexv[il]][iv])
                        dzy.append(2*rfiber)
                        ytrack.append(primary_track[indexv[il]][iv])
        
        #print("x_vec",x_vec)
        #print("y_vec",y_vec)
        #print("zx_vec",zx_vec)
        for ihit in range(len(Layx)):
            Ilayerx[ihit] = Layx[ihit]
        #print("layx",Layx)
        for ihit in range(len(Layy)):
            Ilayery[ihit] = Layy[ihit]                                      
        #print("layy",Layy)

        #reco vs MC track
        x = np.array(x).astype(float)
        y = np.array(y).astype(float)
        
        #print ("x",x)
        #print ("y",y)
        
        '''print ("zx",zx, dzx)
        print ("zy",zy, dzy)'''
        
        '''if np.shape(x)[0] < totLayers and np.shape(y)[0] < totLayers: 
            nlayers[0]  = 0
            trk_flag[0] = 0'''
        #elif np.shape(x)[0] == totLayers and np.shape(y)[0] == totLayers: 
        if np.shape(x)[0] == totLayers and np.shape(y)[0] == totLayers: 

            xtrack = np.array(xtrack).astype(float)                                      
            dx = np.array(dx).astype(float)
            zx = np.array(zx).astype(float)
            dzx= np.array(dzx).astype(float)
            dy = np.array(dy).astype(float)
            ytrack = np.array(ytrack).astype(float)
            zy = np.array(zy).astype(float)
            dzy= np.array(dzy).astype(float)

            nlayers[0]  = totLayers
            nlayersx[0]  = totLayers
            nlayersy[0]  = totLayers
            trk_flag[0] = 1
            x_reco = TGraphErrors(totLayers,zx,x,dzx,dx)
            x_reco.SetTitle("Hor")
            x_true = TGraph(3,zx,xtrack)  
            y_reco = TGraphErrors(totLayers,zy,y,dzy,dy)
            y_reco.SetTitle("Ver")
            y_true = TGraph(3,zy,ytrack)  
            #root_file.cd()    

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
                    x_reco.GetYaxis().SetRangeUser(min(x),max(x))
                    x_reco.Fit("rettax","Q")
                    x_true.Fit("rettaMC","Q")
                    #x_reco.Draw("AP")
                    #x_true.Draw("*")
                if iv == 1:
                    xmin = min(zy)-0.05
                    xmax = max(zy)+0.05
                    rettay = TF1("rettay","[0]*x+[1]",xmin,xmax)
                    rettaMC = TF1("rettaMC","[0]*x+[1]",min(ytrack),max(ytrack))
                    y_reco.SetMarkerStyle(11)
                    #y_reco.GetXaxis().SetLimits(-1,1)
                    y_reco.GetYaxis().SetRangeUser(min(y),max(y))
                    y_reco.Fit("rettay","Q")
                    y_true.Fit("rettaMC","Q")
                    #y_reco.Draw("AP")
                    #y_true.Draw("*")
                #ccc[i][iv].Write() 
            xchi = rettax.GetChisquare() 
            xprob = rettax.GetProb()
            mxrec = (rettax.GetParameter(0))
            qxrec = rettax.GetParameter(1)
            ychi = rettay.GetChisquare() 
            yprob = rettay.GetProb() 
            myrec = (rettay.GetParameter(0))
            qyrec = rettay.GetParameter(1)

            mx[0] = mxrec
            qx[0] = qxrec
            chi2x[0] = xchi
            my[0] = myrec
            qy[0] = qyrec
            chi2y[0] = ychi
            
            crFx,crFy,crFz = evaluateRecoCosineDirector(mxrec, myrec)
            scalarF = ScalarProduct(cx,cy,cz,crFx,crFy,crFz)   
            thetaF = TMath.ACos(scalarF)*180/TMath.Pi()
            
            cosx[0] = crFx
            cosy[0] = crFy
            cosz[0] = crFz

            cosxMC[0] = cx
            cosyMC[0] = cy
            coszMC[0] = cz
            
            angle_deg[0] = thetaF
            
            rx = CalcResiduals(totLayers,x,zx,mxrec,qxrec)
            ry = CalcResiduals(totLayers,y,zy,myrec,qyrec)
            
            for il in range(totLayers):
                for iv in range(nview):
                    if iv == 0: 
                        pos = "x"
                    else:
                        pos = "y"
                    sResiduals = "Residuals{0}".format(pos)
                    exec("%s[il] = r%s[il]" % (sResiduals,pos))
                           
        else:
            nlayers[0] = 0
            nlayersx[0]  = np.shape(x)[0]
            nlayersy[0]  = np.shape(y)[0]

            trk_flag[0] = 0
        
        OutTree.Fill()
        
       
    ###########fine ciclo sugli eventi 
    tfile.Close()
    OutTree.Write()
    root_file.Write()
    root_file.Close()
