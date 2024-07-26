from ROOT import TFile, TTree, TChain, TH1D, TH2D, TMath, TCanvas, TColor, TGraph, TGraphErrors, TF1, TH2F, TObject, gStyle, TLegend
import os
import sys
import time
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

    hNorm = TH2D(hh).Clone(hname)
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

def GenerateLogBinning(nbins, min, max):
    X = np.zeros(nbins+1, float)
    #print(X)
    log_interval = (TMath.Log10(max) - TMath.Log10(min))/nbins
    #print(log_interval)
    for i in range(nbins+1):
        X[i]=pow(10, TMath.Log10(min)+i*log_interval)
    #print (X)
    return X

def calcMCcoord(x0,cx,cz,z):
    x = x0+z*(cx/cz)
    return x

def plot_PSFvsEk (c, hen, hhpsf, emin, emax, eK, root_file, pitch, gap, rFiber, labelY, fflag = 1):
    #fflag ==1 -> Lin fflag==0 ->Log#
    root_file.cd()

    gra50 = TGraph()
    gra68 = TGraph()
    gra95 = TGraph()
    
    nq = 3
    xq = ary.array("d",range(nq))
    yq = ary.array("d",range(nq))
    xq[0] = 0.50
    xq[1] = 0.68
    xq[2] = 0.95
    
    g = 0
    for k in range(len(eK)):
        xg = hen.GetBinContent(k+1) #histo prmary energy
        x  = hhpsf.GetXaxis().GetBinLowEdge(k+1) #histo psf vs energy
        if(xg>0.):
            htmp_py = hhpsf.ProjectionY("htmp_py", k, k+1)
            if(htmp_py.Integral()>0.):
                htmp_py.GetQuantiles(nq, yq, xq)
                gra50.SetPoint(g, x, yq[0])
                gra68.SetPoint(g, x, yq[1])
                gra95.SetPoint(g, x, yq[2])

                #print ("quantiles: ", k, "LowEdge", x, "quant_vec", yq)
                g += 1
                
    c.cd()
    c.SetGrid()
    c.SetLogx()
    c.SetLogy()
    leg = TLegend(0.12,.66,.30,.86)
    leg.SetBorderSize(1)
    hhpsf.SetStats(0)
    hhpsf.GetXaxis().SetLabelFont(42)
    hhpsf.GetXaxis().SetLabelSize(0.04)
    hhpsf.GetXaxis().SetTitle("Primary energy [MeV]")
    hhpsf.GetXaxis().SetTitleSize(0.05)
    hhpsf.GetXaxis().SetTitleFont(22)
    hhpsf.GetXaxis().SetTitleOffset(0.80)
    hhpsf.GetYaxis().SetLabelFont(42)
    hhpsf.GetYaxis().SetLabelSize(0.04)
    hhpsf.GetYaxis().SetTitle(labelY)
    hhpsf.GetYaxis().SetTitleSize(0.05)
    hhpsf.GetYaxis().SetTitleFont(22)
    hhpsf.GetYaxis().SetTitleOffset(0.70)
    hhpsf.GetZaxis().SetLabelFont(42)
    hhpsf.GetZaxis().SetLabelSize(0.04)
    hhpsf.GetZaxis().SetTitleSize(0.05)
    hhpsf.GetZaxis().SetTitleFont(22)
    hhpsf.GetZaxis().SetTitleOffset(0.60)
    hhpsf.Draw("colz")
    gra50.SetLineWidth(4)
    gra68.SetLineWidth(4)
    gra95.SetLineWidth(4)
    gra50.SetLineStyle(7)
    gra68.SetLineStyle(7)
    gra95.SetLineStyle(7)
    gra50.SetLineColor(432-9)
    gra68.SetLineColor(800+1)
    gra95.SetLineColor(820-9)
    gra50.Draw("lsame")
    gra68.Draw("lsame")
    gra95.Draw("lsame")

    #comportamento che tiene conto dello scattering multiplo coulombiano
    h = pow(3.,0.5)*rFiber    
    ll = 2*rFiber + h #mm (module width) 
    X0 = 400 #mm (radiation length)
    ww = pitch/pow(12,0.5) #mm
    pp = gap #mm 
    l_fisica = TF1("l_fisica", "pow(pow([0],2)+pow([1]*1/x,2),0.5)*[2]", emin,emax) #radice(dtheta^2_msc+dtheta^2_risoluzione)
    l_fisica.SetParameter(0, ww/pp)
    l_fisica.SetParameter(1, 13.6*pow(ll/X0,0.5))
    l_fisica.SetParameter(2, 180/TMath.Pi())
    l_fisica.SetLineColor(2)
    l_fisica.Draw("same")
    leg.AddEntry(l_fisica,"theoretical","l")
    
    leg.AddEntry(gra50,"50%","l")
    leg.AddEntry(gra68,"68%","l")
    leg.AddEntry(gra95,"95%","l")
    leg.Draw()
    c.Update()
    c.Write()
    
    return c

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

    overall_track_geom = np.zeros(8) # n.of Layer,Views,Fibers, FibLength, FibRadius, TrackerLength, TrackerWidth
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
                gapLay = float(aaa[8])

                overall_track_geom[0] = layer
                overall_track_geom[1] = view
                overall_track_geom[2] = fiber
                overall_track_geom[3] = lfib
                overall_track_geom[4] = rfib
                overall_track_geom[5] = TrackerLength
                overall_track_geom[6] = TrackerWidth
                overall_track_geom[7] = gapLay

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
        #fin     = "/lustrehome/llorusso/Sim_Geant/Analysis/Tree_NUSESe-Pow_0.1-5000_1000000-evt.root"
        #fgeo = "/lustrehome/llorusso/Sim_Geant/Builds/JobsOutput/NUSESe-Pow_0.1-5000/rootOutput/NUSESe-Pow_0.1-5000_1000000-evt.txt"
        fgeo = sys.argv[2]
        
    except:
        print("ERROR MESSAGE: \n =====> Insert first argument = Tree_ROOTfile.root   &  second argument: GeoFile.txt <=====")
        sys.exit(1)
        
    tfile   = TFile(fin)

    #OutFile
    ###Giuls
    #pitch = (((fin.split("/")[-1]).split(".r")[0]).split("pitch")[-1]).split("u")[0] #um
    #pitch = int(pitch)/1000 #mm
    ###Nadia
    pitch = ((fin.split("/")[-1]).split(".root")[0]).split("pitch")[-1]
    #pitch = sys.argv[3]
    print(pitch)
    OutputName = "Results/Res"+(fgeo.split("/")[-1]).split(".txt")[0] + "_pitch" +str(pitch)+ ".root"
    #OutputName = "CONTROLLORes"+(fgeo.split("/")[-1]).split(".txt")[0] + ".root"
    OutputFile = OutputName
    root_file = TFile(OutputFile, "RECREATE")

    Track_info,Calo_info,gxcfib,gycfib,gzcfib,gxcpix,gycpix,gzcpix = read_GeoFile(fgeo)
    Layers               = Track_info[0]
    Views                = Track_info[1]
    Fibers               = Track_info[2]
    FibLength            = Track_info[3]
    FibRadius            = Track_info[4]
    TrackerLength        = Track_info[5]
    TrackerWidth         = Track_info[6]
    GapLayers            = Track_info[7]
    print ("")
    print ("\t --- Tracker (DIM",TrackerLength,"x",TrackerLength,"x",TrackerWidth,"mm^3) ---")
    print ("[TRACKER_GEO_INFO]: Number of Layers in the Tracker = ", Layers)
    print ("[TRACKER_GEO_INFO]: Number of Views per Layer = ", Views)
    print ("[TRACKER_GEO_INFO]: Number of Fibers per View = ", Fibers)
    print ("[TRACKER_GEO_INFO]: Fiber Length = ",FibLength,"mm","Fiber Radius = ",FibRadius,"mm")
    #print ("x",gxcfib[0][0])
    #print ("y",gycfib[0][1])
    #print ("z",gzcfib[0])

    particle = str((((fgeo.split("/")[-1]).split(".txt")[0]).split("Pow")[0]).split("NUSES")[-1])
    type_track = str((((fgeo.split("/")[-1]).split("_")[2])))
    nview = int(Views)
    totLayers = int(Layers)
    nfiber  = int(Fibers)
    rfiber  = FibRadius #mm
    pitch   = float(pitch)
    ekmax   = (((fgeo.split("/")[-1]).split(".txt")[0]).split("Pow")[-1]).split("-")[1]
    stripNo = int((nfiber*rfiber+rfiber)/pitch) #128
    sigma_x = pitch/2 #pitch/np.sqrt(12)
    print(">>>>",particle, "radius", rfiber, "pitch",pitch, type_track,"<<<<<<")
    ################################
    #ho 128 strip per l'elettronica di lettura: da 0 a 127
    #quindi l'ultima fibra e' letta dalla 127 fibra
    #le fibre sono 130 da 0 a 129 quindi l'ultima fibra non viene letta
    #Complessivamente non leggiamo meta' della fibra zero, meta' della fibra 128 e la fibra 129
    ################################
    fiberyield = 8*0.92 #kev/pe (con cladding)
    PDE = 0.4
    trapeff = 0.054
    PEthresh = 3 #pe per strip

    #definizione canvas
    #ccc = []
    canvas = []
    canvasdelta = []
    canvasposreco = []
    cresidual = []
    for il in range(totLayers):
        canvas.append([])
        canvasdelta.append([])
        canvasposreco.append([])
        cresidual.append([])
        canvas[il]=TCanvas("c{0}".format(il),"c{0}".format(il))
        canvas[il].Divide(2,2)
        canvas[il].SetLogx()
        for iv in range(nview):
            if iv==0:
                v= "x"
            if iv==1:
                v= "y"
            canvasposreco[il].append(TCanvas("c{0}posreco{1}".format(il,v),"c{0}posreco{1}".format(il,v)))
            canvasdelta[il].append(TCanvas("c{0}delta{1}".format(il,v),"c{0}delta{1}".format(il,v)))
            cresidual[il].append(TCanvas("c{0}resid{1}".format(il,v),"c{0}resid{1}".format(il,v)))
            cresidual[il][iv].Divide(2,2)

    cefficiency  = TCanvas("cefficiency","cefficiency")
    ceffTwoPlanes = TCanvas("ceffTwoPlanes","ceffTwoPlanes")
    cnum         = TCanvas("cnum","cnum") 
    cenergy      = TCanvas("cenergy","cenergy")
    cQuant       = TCanvas("cQuant", "cQuant")
    cQuant10     = TCanvas("cQuant10", "cQuant10")
    cQuant2      = TCanvas("cQuant2", "cQuant2")
    cQuant_log   = TCanvas("cQuant_log", "cQuant_log")
    cQuant_log10 = TCanvas("cQuant_log10", "cQuant_log10")
    cQuant_log2  = TCanvas("cQuant_log2", "cQuant_log2")
    cchi2Ek      = TCanvas("cchi2Ek", "cchi2Ek")
    cchi2Ek.Divide(1,2)
    cthetachi2   = TCanvas("cthetachi2", "cthetachi2")
    cthetachi2.Divide(1,2)
    cXchi2Ek_quant = TCanvas("cXchi2Ek_quant","cXchi2Ek_quant")
    cYchi2Ek_quant = TCanvas("cYchi2Ek_quant","cYchi2Ek_quant")

    #definizione istogrammi
    ekmin = 0.1 #MeV
    ekmax = float(ekmax) #MeV
    if ekmax == 20:
        nek = 100
        Nlogbins = 230 #per avere deltaX logaritmico = 100
    elif ekmax==5000:
        nek = 550
        Nlogbins = 470 #per avere deltaX logaritmico = 100

    X = GenerateLogBinning(Nlogbins, ekmin, ekmax)
    if(particle=="proton"):
        bbin = 1001
        mmin = -0.5
        mmax = 1000.5
    else:
        bbin = 5001
        mmin = -0.5
        mmax = 500.5

    hdelta = []
    hposreco = []
    hresidual = []
    hresidual2 = []
    hresidual10 = []
    h_norm = []
    hclu_norm = []
    for il in range(totLayers):
        hdelta.append([])
        hposreco.append([])
        hresidual.append([])
        hresidual2.append([])
        hresidual10.append([])
        h_norm.append([])
        hclu_norm.append([])
        for iv in range(nview):
            if iv == 0:
                v = "x"
            if iv == 1:
                v = "y"
            hposreco[il].append(TH2D("hposreco{0}v{1}".format(il,v),"layer{0} {1}-view;Energy [MeV]; |reco-MC| [mm]".format(il,v),Nlogbins, X, 1000,-50.5,49.5))
            hdelta[il].append(TH2D("hdeltal{0}v{1}".format(il,v),"layer{0} {1}-view;Energy [MeV]; |reco-MC| [mm]".format(il,v),Nlogbins, X, 1000,-50.5,49.5))
            hresidual[il].append(TH1D("hresidual{0}v{1}".format(il,v),"layer{0} {1}-view; {1}i - {1}i_fit [mm]; Number of events".format(il,v),500,-1,1))
            hresidual10[il].append(TH1D("hresidual10{0}v{1}".format(il,v),"layer{0} {1}-view - chi2 < 10; {1}i - {1}i_fit [mm]; Number of events".format(il,v),500,-1,1))
            hresidual2[il].append(TH1D("hresidual2{0}v{1}".format(il,v),"layer{0} {1}-view - chi2 < 2; {1}i - {1}i_fit [mm]; Number of events".format(il,v),500,-1,1))
            h_norm[il].append(TH2D("h_norml{0}v{1}".format(il,v),"layer{0} {1}-view;Energy [MeV]; % Photo Electrons".format(il,v), Nlogbins, X, bbin, mmin, mmax))
            hclu_norm[il].append(TH2D("hclu_norml{0}v{1}".format(il,v),"layer{0} {1}-view;Energy [MeV]; Cluster Size (>3pe)".format(il,v), Nlogbins, X, 21,-0.5,20.5))


    hl0y = TH2D("hl0y","layer0 y-view;Energy [MeV]; Total Photo Electrons",Nlogbins, X,bbin, mmin, mmax) #nek,ekmin,ekmax(991,90,10000,100,-0.5,500)
    hl0x = TH2D("hl0x","layer0 x-view;Energy [MeV]; Total Photo Electrons",Nlogbins, X,bbin, mmin, mmax) 
    hl1y = TH2D("hl1y","layer1 y-view;Energy [MeV]; Total Photo Electrons",Nlogbins, X,bbin, mmin, mmax) 
    hl1x = TH2D("hl1x","layer1 x-view;Energy [MeV]; Total Photo Electrons",Nlogbins, X,bbin, mmin, mmax) 
    hl2y = TH2D("hl2y","layer2 y-view;Energy [MeV]; Total Photo Electrons",Nlogbins, X,bbin, mmin, mmax) 
    hl2x = TH2D("hl2x","layer2 x-view;Energy [MeV]; Total Photo Electrons",Nlogbins, X,bbin, mmin, mmax)

    hl0yclu = TH2D("hl0yclu","layer0 y-view;Energy [MeV]; Cluster Size (>3pe)",Nlogbins, X,21,-0.5,20.5) 
    hl0xclu = TH2D("hl0xclu","layer0 x-view;Energy [MeV]; Cluster Size (>3pe)",Nlogbins, X,21,-0.5,20.5) 
    hl1yclu = TH2D("hl1yclu","layer1 y-view;Energy [MeV]; Cluster Size (>3pe)",Nlogbins, X,21,-0.5,20.5) 
    hl1xclu = TH2D("hl1xclu","layer1 x-view;Energy [MeV]; Cluster Size (>3pe)",Nlogbins, X,21,-0.5,20.5) 
    hl2yclu = TH2D("hl2yclu","layer2 y-view;Energy [MeV]; Cluster Size (>3pe)",Nlogbins, X,21,-0.5,20.5) 
    hl2xclu = TH2D("hl2xclu","layer2 x-view;Energy [MeV]; Cluster Size (>3pe)",Nlogbins, X,21,-0.5,20.5)
    
    hEkDist = TH1D("hEkDist", "Distribution of kinetic energies [MeV]", nek, ekmin, ekmax)
    hEkDist_log = TH1D("hEkDist_log", "Distribution of kinetic energies [MeV] - log binning", Nlogbins, X)
    hNTrigViewVsEk = []
    hEffViewVsEk = []
    for il in range (totLayers):
        hNTrigViewVsEk.append([])
        hEffViewVsEk.append([])
        for iv in range(nview):
            title1 = "hNTrigViewVsEk[" + str(il) + "][" + str(iv) + "]"
            title2 = "Number of triggers vs kinetic energy - " + str(il)+str(iv)
            hNTrigViewVsEk[il].append(TH1D(title1, title2, Nlogbins, X))
            title1 = "hEffViewVsEk[" + str(il) + "][" + str(iv) + "]"
            title2 = "Overall efficiency vs kinetic energy - "  + str(il)+str(iv)
            hEffViewVsEk[il].append(TH1D(title1, title2, Nlogbins, X))

    #hcrFz    = TH1D("hcrFz", "cz Fibers", 500, -1, 0) #-1<cz<0
    hthetaF  = TH1D("hthetaF","PSF Fiber",900,0,90) #180,0,180
    hthetaF.GetXaxis().SetTitle("#theta [#circ]")
    hPSFvsEk = TH2D("hPSFvsEk","PSF Fiber vs primary Energy",nek,ekmin,ekmax,900,0,90)
    hPSFvsEk10 = TH2D("hPSFvsEk10","PSF Fiber vs primary Energy - chi2 < 10",nek,ekmin,ekmax, 900,0,90)
    hPSFvsEk2 = TH2D("hPSFvsEk2","PSF Fiber vs Energy - chi2 < 2",nek,ekmin,ekmax, 900,0,90)
    hPSFvsEk_log = TH2D("hPSFvsEk_log","PSF Fiber vs primary Energy - log binning", Nlogbins, X, 900,0,90)
    hPSFvsEk_log10 = TH2D("hPSFvsEk_log10","PSF Fiber vs primary Energy - log binning - chi2 < 10", Nlogbins, X, 900,0,90)
    hPSFvsEk_log2 = TH2D("hPSFvsEk_log2","PSF Fiber vs primary Energy - log binning - chi2 < 2", Nlogbins, X, 900,0,90)
    hchi2xVsEk = TH2D("hchi2xVsEk", "x-view #chi^{2} vs primary Energy", Nlogbins, X, 3000, 0, 300)
    hchi2yVsEk = TH2D("hchi2yVsEk", "y-view #chi^{2} vs primary Energy", Nlogbins, X, 3000, 0, 300)
    hthetaVschi2x = TH2D("hthetaVschi2x", "x-view #theta vs #chi^{2}", 3000, 0, 300, 900,0,90)
    hthetaVschi2y = TH2D("hthetaVschi2y", "y-view #theta vs #chi^{2}", 3000, 0, 300, 900,0,90)
    hchi2x = TH1D("hchi2x","x view", 3000, 0, 300)
    hchi2y = TH1D("hchi2y","y view", 3000, 0, 300)
    htrigpitch05 = TH1D("htrigpitch05","number of trigger", Nlogbins, X)
    heffpitch05 = TH1D("heffpitch05","Overall efficiency vs kinetic energy", Nlogbins, X)

    tree = tfile.Get("ClusterTree")
    nhit = tree.GetEntries()
    #nhit = 1000

    root_file.cd()
    
    theta=[]
    eK=[]
    
    for i in range(nhit):
        if i%1000==0:
            print("Event: ", i)

        #print("Event:", i)
        tree.GetEntry(i)
        
        PrimaryEnergy = (tree.energy)*1e-3
        hEkDist.Fill(PrimaryEnergy, 1.0)
        hEkDist_log.Fill(PrimaryEnergy, 1.0)
        
        trk_flag = tree.trk_flag
        if trk_flag == 1:
            thetaF = tree.angle_deg
            thetaF = list(thetaF)[0]
            hthetaF.Fill(thetaF)
            hPSFvsEk.Fill(PrimaryEnergy, thetaF)
            hPSFvsEk_log.Fill(PrimaryEnergy, thetaF)
            theta.append(thetaF)
            eK.append(PrimaryEnergy)

            chi2X = tree.chi2x
            chi2X = list(chi2X)[0]
            chi2Y = tree.chi2y
            chi2Y = list(chi2Y)[0]
            hchi2x.Fill(chi2X)
            hchi2y.Fill(chi2Y)
            hchi2xVsEk.Fill(PrimaryEnergy, chi2X)
            hchi2yVsEk.Fill(PrimaryEnergy, chi2Y)
            hthetaVschi2x.Fill(chi2X,thetaF)
            hthetaVschi2y.Fill(chi2Y,thetaF)
            if (chi2X < 10 and chi2Y < 10):
                hPSFvsEk10.Fill(PrimaryEnergy, thetaF)
                hPSFvsEk_log10.Fill(PrimaryEnergy, thetaF)
            if (chi2X < 2 and chi2Y < 2):
                hPSFvsEk2.Fill(PrimaryEnergy, thetaF)
                hPSFvsEk_log2.Fill(PrimaryEnergy, thetaF)

        Strips = []
        PEStrips = []
        vtrig= np.zeros((totLayers,nview))
        efftrig = 0
        residual = []

        xpos = tree.x_vec
        xposMC = tree.xMC_vec
        ypos = tree.y_vec
        yposMC = tree.yMC_vec

        nlayX = tree.nlayersx
        nlayY = tree.nlayersy

        #print("nlayX", nlayX )
        #print("xpos", list(xpos))
        #print("xposMC", list(xposMC))

        LayX = tree.Ilayerx
        #print("LayX", list(LayX))
        LayY = tree.Ilayery
        #print("LayY", list(LayY))

        for il in range(nlayX):
                    
            xpos_reco = xpos[il] - xposMC[il]
            hposreco[LayX[il]][0].Fill(PrimaryEnergy, xpos_reco)

        for il in range(nlayY):
                    
            ypos_reco = ypos[il] - yposMC[il]
            hposreco[LayY[il]][1].Fill(PrimaryEnergy, ypos_reco)

        for il in range(totLayers):
            Strips.append([])
            PEStrips.append([])
            residual.append([])
            for iv in range(nview):
                Strips[il].append([])
                PEStrips[il].append([])
                residual[il].append([])
                if iv == 0: 
                    v = "Hor"
                    p = "x"
                else:
                    v = "Ver" 
                    p = "y"
                sstrip     = "strip{0}{1}".format(v,il)
                sPEstrip   = "PEstrip{0}{1}".format(v,il)
                
                exec("Strips[il][iv]   = tree.%s" % sstrip)
                exec("PEStrips[il][iv]   = tree.%s" % sPEstrip)
                
                for istrip in range(len(Strips[il][iv])-1):
                    if (PEStrips[il][iv][istrip] >= PEthresh and PEStrips[il][iv][istrip+1] >= PEthresh):
                        vtrig[il][iv] = 1
                        
                if vtrig[il][iv]>0:
                    hNTrigViewVsEk[il][iv].Fill(PrimaryEnergy, 1.0)

                if trk_flag == 1:
                    sx = "{0}_vec".format(p)
                    sxMC = "{0}MC_vec".format(p)
                    exec("x  = tree.%s" % sx)
                    exec("xMC  = tree.%s" % sxMC)

                    delta = x[il] - xMC[il]
                    hdelta[il][iv].Fill(PrimaryEnergy,delta)

                    #sres = "{0}_residual".format(p)
                    sres = "Residuals{0}".format(p)
                    exec("res  = tree.%s" % sres)
                    
                    hresidual[il][iv].Fill(res[il])
                    if (chi2X < 10 and chi2Y < 10):
                        hresidual10[il][iv].Fill(res[il])
                    if (chi2X < 2 and chi2Y < 2):
                        hresidual2[il][iv].Fill(res[il])

        #overall efficiency for pitch = 0.5 mm
        #if(pitch == 0.50):
        if(nlayX >= 2 and nlayY >= 2):
                efftrig = 1

        if efftrig > 0:
            htrigpitch05.Fill(PrimaryEnergy, 1.0)
            efftrig = 0
    ##########fine ciclo eventi#########

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

    '''h1Ngen = TH1F("h1Ngen","; MC Energy (MeV); Number of Generated Events", net, Et)
    for i in range(net1):
        h1Ngen.Fill(eK[i], nhit)'''

    cQuant = plot_PSFvsEk(cQuant, hEkDist, hPSFvsEk, ekmin, ekmax, eK, root_file,pitch,GapLayers,FibRadius,"#theta [#circ]", fflag = 1)
    cQuant10 = plot_PSFvsEk(cQuant10,hEkDist, hPSFvsEk10, ekmin, ekmax, eK, root_file,pitch,GapLayers,FibRadius,"#theta [#circ]", fflag = 1)
    cQuant2 = plot_PSFvsEk(cQuant2,hEkDist, hPSFvsEk2, ekmin, ekmax, eK, root_file,pitch,GapLayers,FibRadius,"#theta [#circ]", fflag = 1)
    cQuant_log = plot_PSFvsEk(cQuant_log, hEkDist_log, hPSFvsEk_log, ekmin, ekmax, eK, root_file,pitch,GapLayers,FibRadius,"#theta [#circ]", fflag = 0)
    cQuant_log10 = plot_PSFvsEk(cQuant_log10, hEkDist_log, hPSFvsEk_log10, ekmin, ekmax, eK, root_file,pitch,GapLayers,FibRadius,"#theta [#circ]", fflag = 0)
    cQuant_log2 = plot_PSFvsEk(cQuant_log2, hEkDist_log, hPSFvsEk_log2, ekmin, ekmax, eK, root_file,pitch,GapLayers,FibRadius,"#theta [#circ]", fflag = 0)

    cenergy.cd()
    hEkDist_log.GetXaxis().SetLabelFont(42)
    hEkDist_log.GetXaxis().SetLabelSize(0.04)
    hEkDist_log.GetXaxis().SetTitle("Primary energy [MeV]")
    hEkDist_log.GetXaxis().SetTitleSize(0.06)
    hEkDist_log.GetXaxis().SetTitleFont(22)
    hEkDist_log.GetXaxis().SetTitleOffset(0.80)
    hEkDist_log.GetYaxis().SetLabelFont(42)
    hEkDist_log.GetYaxis().SetLabelSize(0.04)
    hEkDist_log.GetYaxis().SetTitle("Number of events")
    hEkDist_log.GetYaxis().SetTitleSize(0.06)
    hEkDist_log.GetYaxis().SetTitleFont(22)
    hEkDist_log.GetYaxis().SetTitleOffset(0.70)
    hEkDist_log.Draw()

    color=[[416+1, 820-7], [632-6, 800-2], [616-2, 880+6]]
    label=[]
    hEffViewVsEk[0][0].SetTitle("Overall efficiency vs Ek")
    hEffViewVsEk[0][0].SetStats(0)
    legend = TLegend(0.66,.32,.89,.59)
    legend.SetBorderSize(1)
    legend0 = TLegend(0.66,.32,.89,.59)
    legend0.SetBorderSize(1)

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
            sTpeStrip   = "TpeStrip{0}{1}".format(v,il)
            sPEstrip    = "PEstrip{0}{1}".format(v,il)
            
            if(iv%2==0):
                icanv = iv + 2
                icanv_c = iv + 4
            else: 
                icanv = iv
                icanv_c = iv + 2

            gPad1 = canvas[il].cd(icanv)   
            #gPad1.SetLogx()              
            title = sTpeStrip + ":energy*1e-3>>hl"+str(il)+coo
            cut = ""
            tree.Draw(title, cut, "colz")
            histo = "hl"+ str(il)+coo
            exec("h_norm[il][iv]=NormalizeX(%s)" %histo)
            h_norm[il][iv].GetYaxis().SetTitleFont(22)
            h_norm[il][iv].GetYaxis().SetTitleSize(0.06)
            h_norm[il][iv].GetYaxis().SetLabelFont(42)
            h_norm[il][iv].GetYaxis().SetLabelSize(0.04)
            h_norm[il][iv].GetYaxis().SetTitleOffset(0.70)
            h_norm[il][iv].GetXaxis().SetTitleFont(22)
            h_norm[il][iv].GetXaxis().SetTitleSize(0.06)
            h_norm[il][iv].GetXaxis().SetTitleOffset(0.75)
            h_norm[il][iv].GetXaxis().SetLabelFont(42)
            h_norm[il][iv].GetXaxis().SetLabelSize(0.04)
            h_norm[il][iv].GetZaxis().SetTitleFont(22)
            h_norm[il][iv].GetZaxis().SetTitleSize(0.06)
            h_norm[il][iv].GetZaxis().SetTitleOffset(0.75)
            h_norm[il][iv].GetZaxis().SetLabelFont(42)
            h_norm[il][iv].GetZaxis().SetLabelSize(0.04)
            h_norm[il][iv].Draw("colz")
            
            gPad2 = canvas[il].cd(icanv_c)   
            #gPad2.SetLogx()            
            title = sCsiz + ":energy*1e-3>>hl"+str(il)+coo +"clu"
            cut = sPEstrip + ">=3 || " + sCsiz + " == 0"
            tree.Draw(title, cut, "colz")
            histo_clu = "hl"+ str(il)+coo+"clu"
            exec("hclu_norm[il][iv]=NormalizeX(%s)" %histo_clu)
            hclu_norm[il][iv].GetYaxis().SetTitleFont(22)
            hclu_norm[il][iv].GetYaxis().SetTitleSize(0.06)
            hclu_norm[il][iv].GetYaxis().SetLabelFont(42)
            hclu_norm[il][iv].GetYaxis().SetLabelSize(0.04)
            hclu_norm[il][iv].GetYaxis().SetTitleOffset(0.70)
            hclu_norm[il][iv].GetXaxis().SetTitleFont(22)
            hclu_norm[il][iv].GetXaxis().SetTitleSize(0.06)
            hclu_norm[il][iv].GetXaxis().SetTitleOffset(0.75)
            hclu_norm[il][iv].GetXaxis().SetLabelFont(42)
            hclu_norm[il][iv].GetXaxis().SetLabelSize(0.04)
            hclu_norm[il][iv].GetZaxis().SetTitleFont(22)
            hclu_norm[il][iv].GetZaxis().SetTitleSize(0.06)
            hclu_norm[il][iv].GetZaxis().SetTitleOffset(0.75)
            hclu_norm[il][iv].GetZaxis().SetLabelFont(42)
            hclu_norm[il][iv].GetZaxis().SetLabelSize(0.04)
            hclu_norm[il][iv].Draw("colz")

            label[il].append("Layer {0} View {1}".format(il,coo))
            cnum.cd()
            if il == 0 and iv == 0:
                hNTrigViewVsEk[il][iv].Draw("")
            else: hNTrigViewVsEk[il][iv].Draw("same")

            #if pitch == 0.25:
            cefficiency.cd()
            cefficiency.SetGrid()

            hEffViewVsEk[il][iv].Divide(hNTrigViewVsEk[il][iv],hEkDist_log, Nlogbins, Nlogbins,"B")

            hEffViewVsEk[il][iv].GetYaxis().SetTitle("Efficiency")
            hEffViewVsEk[il][iv].GetYaxis().SetRangeUser(0, 1.05)
            hEffViewVsEk[il][iv].GetYaxis().SetTitleFont(22)
            hEffViewVsEk[il][iv].GetYaxis().SetLabelFont(42)
            hEffViewVsEk[il][iv].GetYaxis().SetLabelSize(0.04)
            hEffViewVsEk[il][iv].GetYaxis().SetTitleSize(0.05)
            hEffViewVsEk[il][iv].GetYaxis().SetTitleOffset(0.70) 
            hEffViewVsEk[il][iv].GetXaxis().SetTitle("Energy [MeV]")
            hEffViewVsEk[il][iv].GetXaxis().SetTitleFont(22)
            hEffViewVsEk[il][iv].GetXaxis().SetLabelFont(42)
            hEffViewVsEk[il][iv].GetXaxis().SetLabelSize(0.04) 
            hEffViewVsEk[il][iv].GetXaxis().SetTitleSize(0.05)
            hEffViewVsEk[il][iv].GetXaxis().SetTitleOffset(0.80)      

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
            legend.Draw()

            #elif pitch == 0.5:
            heffpitch05.Divide(htrigpitch05,hEkDist_log, Nlogbins, Nlogbins,"B")

            ceffTwoPlanes.cd()
            ceffTwoPlanes.SetGrid()

            heffpitch05.SetStats(0)
            heffpitch05.GetYaxis().SetTitle("Efficiency")
            heffpitch05.GetYaxis().SetRangeUser(0, 1.05)
            heffpitch05.GetYaxis().SetTitleFont(22)
            heffpitch05.GetYaxis().SetLabelFont(42)
            heffpitch05.GetYaxis().SetLabelSize(0.04)
            heffpitch05.GetYaxis().SetTitleSize(0.05)
            heffpitch05.GetYaxis().SetTitleOffset(0.70) 
            heffpitch05.GetXaxis().SetTitle("Energy [MeV]")
            heffpitch05.GetXaxis().SetTitleFont(22)
            heffpitch05.GetXaxis().SetLabelFont(42)
            heffpitch05.GetXaxis().SetLabelSize(0.04) 
            heffpitch05.GetXaxis().SetTitleSize(0.05)
            heffpitch05.GetXaxis().SetTitleOffset(0.80)      

            heffpitch05.SetLineColor(color[il][iv])
            heffpitch05.SetLineWidth(1)
            heffpitch05.SetMarkerColor(color[il][iv])
            heffpitch05.SetMarkerStyle(20)

            heffpitch05.Draw("")
            ceffTwoPlanes.Update()

    for il in range(totLayers):
        for iv in range(nview):
            canvasposreco[il][iv].cd()
            hposreco[il][iv].SetStats(0)
            hposreco[il][iv].GetYaxis().SetTitleFont(22)
            hposreco[il][iv].GetYaxis().SetTitleSize(0.05)
            hposreco[il][iv].GetYaxis().SetLabelFont(42)
            hposreco[il][iv].GetYaxis().SetLabelSize(0.04)
            hposreco[il][iv].GetYaxis().SetTitleOffset(0.70)
            hposreco[il][iv].GetXaxis().SetTitleFont(22)
            hposreco[il][iv].GetXaxis().SetTitleSize(0.05)
            hposreco[il][iv].GetXaxis().SetTitleOffset(0.75)
            hposreco[il][iv].GetXaxis().SetLabelFont(42)
            hposreco[il][iv].GetXaxis().SetLabelSize(0.04)
            hposreco[il][iv].GetZaxis().SetTitleFont(22)
            hposreco[il][iv].GetZaxis().SetTitleSize(0.05)
            hposreco[il][iv].GetZaxis().SetTitleOffset(0.75)
            hposreco[il][iv].GetZaxis().SetLabelFont(42)
            hposreco[il][iv].GetZaxis().SetLabelSize(0.04)
            hposreco[il][iv].Draw("colz")

            canvasdelta[il][iv].cd()
            hdelta[il][iv].SetStats(0)
            hdelta[il][iv].GetYaxis().SetTitleFont(22)
            hdelta[il][iv].GetYaxis().SetTitleSize(0.05)
            hdelta[il][iv].GetYaxis().SetLabelFont(42)
            hdelta[il][iv].GetYaxis().SetLabelSize(0.04)
            hdelta[il][iv].GetYaxis().SetTitleOffset(0.70)
            hdelta[il][iv].GetXaxis().SetTitleFont(22)
            hdelta[il][iv].GetXaxis().SetTitleSize(0.05)
            hdelta[il][iv].GetXaxis().SetTitleOffset(0.75)
            hdelta[il][iv].GetXaxis().SetLabelFont(42)
            hdelta[il][iv].GetXaxis().SetLabelSize(0.04)
            hdelta[il][iv].GetZaxis().SetTitleFont(22)
            hdelta[il][iv].GetZaxis().SetTitleSize(0.05)
            hdelta[il][iv].GetZaxis().SetTitleOffset(0.75)
            hdelta[il][iv].GetZaxis().SetLabelFont(42)
            hdelta[il][iv].GetZaxis().SetLabelSize(0.04)
            hdelta[il][iv].Draw("colz")
            cresidual[il][iv].cd(1)
            gauss = TF1("gauss","gaus",-1,+1)
            gauss.SetLineColor(632+1)
            Gamma = 0.2
            #a = Gamma/(2*TMath.Pi())
            a = Gamma/2
            lorentz = TF1("lorentz","[0]*1/({0}*(pow(x-[1],2)+pow([2],2)))".format(TMath.Pi()),-1,+1)
            lorentz.SetLineColor(416+2)
            lorentz.SetParameter(0, a)
            lorentz.SetParameter(2, a)
            crystball = TF1("crystball","crystalball",-1,+1)
            #crystball = TF1("crystball", "ROOT::Math::crystalball_function(x, [0], [1], [2], [3])", -1,+1)
            crystball.SetLineColor(616+4)
            crystball.SetParameters(1, 0, 5, 1.5, 2)
            '''crystball.SetParameter(0, 1.5)#alpha
            crystball.SetParameter(1, )#N
            crystball.SetParameter(2, 0)#mean
            crystball.SetParameter(3, 0.1)#sigma'''
            hresidual[il][iv].GetYaxis().SetTitleFont(22)
            hresidual[il][iv].GetYaxis().SetTitleSize(0.06)
            hresidual[il][iv].GetYaxis().SetLabelFont(42)
            hresidual[il][iv].GetYaxis().SetLabelSize(0.04)
            hresidual[il][iv].GetYaxis().SetTitleOffset(0.70)
            hresidual[il][iv].GetXaxis().SetTitleFont(22)
            hresidual[il][iv].GetXaxis().SetTitleSize(0.06)
            hresidual[il][iv].GetXaxis().SetTitleOffset(0.75)
            hresidual[il][iv].GetXaxis().SetLabelFont(42)
            hresidual[il][iv].GetXaxis().SetLabelSize(0.04)
            hresidual[il][iv].Fit("gauss","QR")
            #hresidual[il][iv].Fit("lorentz","QR+")
            #hresidual[il][iv].Fit("crystball","QR+")
            hresidual[il][iv].Draw()
            cresidual[il][iv].cd(2)
            hresidual2[il][iv].GetYaxis().SetTitleFont(22)
            hresidual2[il][iv].GetYaxis().SetTitleSize(0.06)
            hresidual2[il][iv].GetYaxis().SetLabelFont(42)
            hresidual2[il][iv].GetYaxis().SetLabelSize(0.04)
            hresidual2[il][iv].GetYaxis().SetTitleOffset(0.70)
            hresidual2[il][iv].GetXaxis().SetTitleFont(22)
            hresidual2[il][iv].GetXaxis().SetTitleSize(0.06)
            hresidual2[il][iv].GetXaxis().SetTitleOffset(0.75)
            hresidual2[il][iv].GetXaxis().SetLabelFont(42)
            hresidual2[il][iv].GetXaxis().SetLabelSize(0.04)
            hresidual2[il][iv].Fit("gauss","QR")
            #hresidual2[il][iv].Fit("lorentz","QR+")
            #hresidual2[il][iv].Fit("crystball","QR+")
            hresidual2[il][iv].Draw()
            cresidual[il][iv].cd(3)
            hresidual10[il][iv].GetYaxis().SetTitleFont(22)
            hresidual10[il][iv].GetYaxis().SetTitleSize(0.06)
            hresidual10[il][iv].GetYaxis().SetLabelFont(42)
            hresidual10[il][iv].GetYaxis().SetLabelSize(0.04)
            hresidual10[il][iv].GetYaxis().SetTitleOffset(0.70)
            hresidual10[il][iv].GetXaxis().SetTitleFont(22)
            hresidual10[il][iv].GetXaxis().SetTitleSize(0.06)
            hresidual10[il][iv].GetXaxis().SetTitleOffset(0.75)
            hresidual10[il][iv].GetXaxis().SetLabelFont(42)
            hresidual10[il][iv].GetXaxis().SetLabelSize(0.04)
            hresidual10[il][iv].Fit("gauss","QR")
            #hresidual10[il][iv].Fit("lorentz","QR+")
            #hresidual10[il][iv].Fit("crystball","QR+")
            hresidual10[il][iv].Draw()
            canvasdelta[il][iv].Write()
            cresidual[il][iv].Write()

    cchi2Ek.cd(1)
    hchi2xVsEk.GetYaxis().SetTitleFont(22)
    hchi2xVsEk.GetYaxis().SetTitleSize(0.06)
    hchi2xVsEk.GetYaxis().SetLabelFont(42)
    hchi2xVsEk.GetYaxis().SetLabelSize(0.04)
    hchi2xVsEk.GetYaxis().SetTitleOffset(0.70)
    hchi2xVsEk.GetYaxis().SetTitle("#chi2")
    hchi2xVsEk.GetXaxis().SetTitleFont(22)
    hchi2xVsEk.GetXaxis().SetTitleSize(0.06)
    hchi2xVsEk.GetXaxis().SetTitleOffset(0.75)
    hchi2xVsEk.GetXaxis().SetLabelFont(42)
    hchi2xVsEk.GetXaxis().SetLabelSize(0.04)
    hchi2xVsEk.GetXaxis().SetTitle("Primary energy [MeV]")
    hchi2xVsEk.GetZaxis().SetTitleFont(22)
    hchi2xVsEk.GetZaxis().SetTitleSize(0.06)
    hchi2xVsEk.GetZaxis().SetTitleOffset(0.75)
    hchi2xVsEk.GetZaxis().SetLabelFont(42)
    hchi2xVsEk.GetZaxis().SetLabelSize(0.04)
    hchi2xVsEk.Draw("colz")
    cchi2Ek.cd(2)
    hchi2yVsEk.GetYaxis().SetTitleFont(22)
    hchi2yVsEk.GetYaxis().SetTitleSize(0.06)
    hchi2yVsEk.GetYaxis().SetLabelFont(42)
    hchi2yVsEk.GetYaxis().SetLabelSize(0.04)
    hchi2yVsEk.GetYaxis().SetTitleOffset(0.70)
    hchi2yVsEk.GetYaxis().SetTitle("#chi2")
    hchi2yVsEk.GetXaxis().SetTitleFont(22)
    hchi2yVsEk.GetXaxis().SetTitleSize(0.06)
    hchi2yVsEk.GetXaxis().SetTitleOffset(0.75)
    hchi2yVsEk.GetXaxis().SetLabelFont(42)
    hchi2yVsEk.GetXaxis().SetLabelSize(0.04)
    hchi2yVsEk.GetXaxis().SetTitle("Primary energy [MeV]")
    hchi2yVsEk.GetZaxis().SetTitleFont(22)
    hchi2yVsEk.GetZaxis().SetTitleSize(0.06)
    hchi2yVsEk.GetZaxis().SetTitleOffset(0.75)
    hchi2yVsEk.GetZaxis().SetLabelFont(42)
    hchi2yVsEk.GetZaxis().SetLabelSize(0.04)
    hchi2yVsEk.Draw("colz")
    
    cXchi2Ek_quant = plot_PSFvsEk (cXchi2Ek_quant, hEkDist_log, hchi2xVsEk, ekmin, ekmax, eK, root_file,pitch,GapLayers,FibRadius, "#chi^{2}", fflag = 0)
    cYchi2Ek_quant = plot_PSFvsEk (cYchi2Ek_quant, hEkDist_log, hchi2yVsEk, ekmin, ekmax, eK, root_file,pitch,GapLayers,FibRadius, "#chi^{2}", fflag = 0)
    
    cthetachi2.cd(1)
    hthetaVschi2x.GetYaxis().SetTitleFont(22)
    hthetaVschi2x.GetYaxis().SetTitleSize(0.06)
    hthetaVschi2x.GetYaxis().SetLabelFont(42)
    hthetaVschi2x.GetYaxis().SetLabelSize(0.04)
    hthetaVschi2x.GetYaxis().SetTitleOffset(0.70)
    hthetaVschi2x.GetYaxis().SetTitle("#theta [#circ]")
    hthetaVschi2x.GetXaxis().SetTitleFont(22)
    hthetaVschi2x.GetXaxis().SetTitleSize(0.06)
    hthetaVschi2x.GetXaxis().SetTitleOffset(0.75)
    hthetaVschi2x.GetXaxis().SetLabelFont(42)
    hthetaVschi2x.GetXaxis().SetLabelSize(0.04)
    hthetaVschi2x.GetXaxis().SetTitle("#chi2")
    hthetaVschi2x.GetZaxis().SetTitleFont(22)
    hthetaVschi2x.GetZaxis().SetTitleSize(0.06)
    hthetaVschi2x.GetZaxis().SetTitleOffset(0.75)
    hthetaVschi2x.GetZaxis().SetLabelFont(42)
    hthetaVschi2x.GetZaxis().SetLabelSize(0.04)
    hthetaVschi2x.Draw("colz")
    cthetachi2.cd(2)
    hthetaVschi2y.GetYaxis().SetTitleFont(22)
    hthetaVschi2y.GetYaxis().SetTitleSize(0.06)
    hthetaVschi2y.GetYaxis().SetLabelFont(42)
    hthetaVschi2y.GetYaxis().SetLabelSize(0.04)
    hthetaVschi2y.GetYaxis().SetTitleOffset(0.70)
    hthetaVschi2y.GetYaxis().SetTitle("#theta [#circ]")
    hthetaVschi2y.GetXaxis().SetTitleFont(22)
    hthetaVschi2y.GetXaxis().SetTitleSize(0.06)
    hthetaVschi2y.GetXaxis().SetTitleOffset(0.75)
    hthetaVschi2y.GetXaxis().SetLabelFont(42)
    hthetaVschi2y.GetXaxis().SetLabelSize(0.04)
    hthetaVschi2y.GetXaxis().SetTitle("#chi2")
    hthetaVschi2y.GetZaxis().SetTitleFont(22)
    hthetaVschi2y.GetZaxis().SetTitleSize(0.06)
    hthetaVschi2y.GetZaxis().SetTitleOffset(0.75)
    hthetaVschi2y.GetZaxis().SetLabelFont(42)
    hthetaVschi2y.GetZaxis().SetLabelSize(0.04)
    hthetaVschi2y.Draw("colz")    

    cefficiency.Write()
    ceffTwoPlanes.Write()
    canvas[0].Write()
    canvas[1].Write()
    canvas[2].Write()
    cchi2Ek.Write()
    cthetachi2.Write()
    
    root_file.Write()
    root_file.Close()
