#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from ROOT import TFile, TTree, TChain, TH1D, TH2D, TMath, TCanvas, TH1F,TRandom3,TGraphErrors, TGraph,TF1

import sys
import array as ary
import numpy as np
import pandas as pd
import os
import time
import glob
import math
import pdb

#parquet libraries
import pyarrow as pa
import pyarrow.parquet as pq

from sklearn.metrics import r2_score

import functions

ran = TRandom3()

if __name__ == "__main__":

    mainDir = sys.argv[3]#"/lustrehome/federica1992/Nuses/Analysis/"
    
    #InputFile
    try:
        fin = sys.argv[1] #senza .root
        #fin     = "/lustrehome/llorusso/Sim_Geant/Builds/JobsOutput/NUSESe-Pow_0.1-5000-0.375_onAxis/rootOutput/NUSESe-Pow_0.1-5000-0.375_onAxis_1000000-evt.root"
        #fgeo    = "/lustrehome/llorusso/Sim_Geant/Builds/JobsOutput/NUSESe-Pow_0.1-5000-0.375_onAxis/rootOutput/NUSESe-Pow_0.1-5000-0.375_onAxis_1000000-evt.txt"
        fgeo = fin
        
    except:
        print("ERROR MESSAGE: \n =====> Insert first argument = ROOTfile.root - second argument = Geofile.txt - third argument = pitch (mm)<=====")
        sys.exit(1)
        
    fin  = fin + ".root"
    fgeo = fgeo + ".txt"
    tfile   = TFile(fin)
    
    print('fin ',fin)
    #OutTreeFile
    pitch   = float(sys.argv[2])#0.25 #mm
    ooo = (fin.split("/")[-1]).split(".r")[0]
    print("ooo ",ooo)
    OutputName = "Results/Tree/Tree_"+str(ooo)+"_pitch"+str(pitch)+".root"
    OutputFile = mainDir+OutputName
    print("OUTPUT:", OutputFile)
    root_file = TFile(OutputFile, "RECREATE")
    
    OutTree           = TTree("ClusterTree", "ClusterTree")
    #BigClusterTree    = TTree("BigClusterTree", "BigClusterTree")
    
    calo_type= sys.argv[4]
    save_csv=0
    
    if (len(sys.argv)>5):      
        save_csv=int(sys.argv[5])
        if(save_csv==1):
             print('********************you are saving the csv files********************')
    
    save_parquet=0
    if(len(sys.argv)>6):
        save_parquet=int(sys.argv[6])
        if(save_parquet==1):
           print('********************you are saving the parquet files********************')
        
    
    save_npz=0
    if(len(sys.argv)>7):        
        save_npz=int(sys.argv[7])
        if(save_npz==1):
            print('********************you are saving the npz files********************')
            
    #OutCSV file
    folder_path_csv=mainDir+'Results/CSVFile/'
    if not os.path.exists(folder_path_csv):
        os.makedirs(folder_path_csv)
 
     #Out npz file
    folder_path_npz=mainDir+'Results/npzFile/'
    if not os.path.exists(folder_path_npz):
        os.makedirs(folder_path_npz)
        
    print('_________________________________')
    if calo_type=='HERD':
        print("HERD calo type")
        #csv files with coordinates after clustering procedure and max cl search
        OutputNameCSV_x = "TrkDataHERD_viewx_afterCL_max_"+str(ooo)+"_pitch"+str(pitch)+".csv"
        OutputNameCSV_x = folder_path_csv+OutputNameCSV_x
        print('OutputNameCSVHERD_x ',OutputNameCSV_x)
        OutputNameCSV_y = "TrkDataHERD_viewy_afterCL_max_"+str(ooo)+"_pitch"+str(pitch)+".csv"
        OutputNameCSV_y = folder_path_csv+OutputNameCSV_y
        print('OutputNameCSVHERD_y ',OutputNameCSV_y)
        print('_________________________________')
        #csv files with coordinates before clustering procedure
        OutputNameCSV_x_nocl = "TrkDataHERD_viewx_raw_"+str(ooo)+"_pitch"+str(pitch)+".csv"
        OutputNameCSV_x_nocl = folder_path_csv+OutputNameCSV_x_nocl
        print('OutputNameCSVHERD_x_nocl ',OutputNameCSV_x_nocl)
        OutputNameCSV_y_nocl = "TrkDataHERD_viewy_raw_"+str(ooo)+"_pitch"+str(pitch)+".csv"
        OutputNameCSV_y_nocl = folder_path_csv+OutputNameCSV_y_nocl
        print('OutputNameCSVHERD_y_nocl ',OutputNameCSV_y_nocl)  

        print('_________________________________')
        #csv files with coordinates after clustering procedure
        OutputNameCSV_x_cl = "TrkDataHERD_viewx_cl_"+str(ooo)+"_pitch"+str(pitch)#+".csv"
        OutputNameCSV_x_cl = folder_path_csv+OutputNameCSV_x_cl
        print('OutputNameCSVHERD_x_cl ',OutputNameCSV_x_cl)
        OutputNameCSV_y_cl = "TrkDataHERD_viewy_cl_"+str(ooo)+"_pitch"+str(pitch)#+".csv"
        OutputNameCSV_y_cl = folder_path_csv+OutputNameCSV_y_cl
        print('OutputNameCSVHERD_y_cl ',OutputNameCSV_y_cl)  
        
        
    
    if calo_type=='CsI':
        
        print('_________________________________')
        #csv files with coordinates after clustering procedure and max cl search
        OutputNameCSV_x = "TrkData_viewx_afterCL_max_"+str(ooo)+"_pitch"+str(pitch)+".csv"
        OutputNameCSV_x = folder_path_csv+OutputNameCSV_x
        print('OutputNameCSV_x ',OutputNameCSV_x)
        OutputNameCSV_y = "TrkData_viewy_afterCL_max_"+str(ooo)+"_pitch"+str(pitch)+".csv"
        OutputNameCSV_y = folder_path_csv+OutputNameCSV_y
        print('OutputNameCSV_y ',OutputNameCSV_y)
        
        print('_________________________________')
        #csv files with coordinates before clustering procedure
        OutputNameCSV_x_nocl = "TrkData_viewx_raw_"+str(ooo)+"_pitch"+str(pitch)+".csv"
        OutputNameCSV_x_nocl = folder_path_csv+OutputNameCSV_x_nocl
        print('OutputNameCSV_x_nocl ',OutputNameCSV_x_nocl)
        OutputNameCSV_y_nocl = "TrkData_viewy_raw_"+str(ooo)+"_pitch"+str(pitch)+".csv"
        OutputNameCSV_y_nocl = folder_path_csv+OutputNameCSV_y_nocl
        print('OutputNameCSV_y_nocl ',OutputNameCSV_y_nocl)  
        
        print('_________________________________')
        #csv files with coordinates after clustering procedure
        OutputNameCSV_x_cl = "TrkData_viewx_cl_"+str(ooo)+"_pitch"+str(pitch)+".csv"
        OutputNameCSV_x_cl = folder_path_csv+OutputNameCSV_x_cl
        print('OutputNameCSV_x_cl ',OutputNameCSV_x_cl)
        OutputNameCSV_y_cl = "TrkData_viewy_cl_"+str(ooo)+"_pitch"+str(pitch)+".csv"
        OutputNameCSV_y_cl = folder_path_csv+OutputNameCSV_y_cl
        print('OutputNameCSV_y_cl ',OutputNameCSV_y_cl)  

    
    #saving numpy array for HERD or Nuses
    if calo_type=='HERD':
        tr_np_name = "TrkDataHERD_raw_"+str(ooo)+"_pitch"+str(pitch)+".npz"  
        tr_np_name = folder_path_npz+tr_np_name
        tr_np_name_true = "TrkDataHERD_true_"+str(ooo)+"_pitch"+str(pitch)+".npz"  
        tr_np_name_true = folder_path_npz+tr_np_name_true
        tr_np_name_cl = "TrkDataHERD_cl_"+str(ooo)+"_pitch"+str(pitch)+".npz"  
        tr_np_name_cl = folder_path_npz+tr_np_name_cl
        print('tr_np_name ',tr_np_name)
        print('tr_np_name_true ',tr_np_name_true)
        print('tr_np_name_cl ',tr_np_name_cl)
    if calo_type=='CsI':
        tr_np_name = "TrkData_raw_"+str(ooo)+"_pitch"+str(pitch)+".npz"  
        tr_np_name = folder_path_npz+tr_np_name
        tr_np_name_true = "TrkData_true_"+str(ooo)+"_pitch"+str(pitch)+".npz"  
        tr_np_name_true = folder_path_npz+tr_np_name_true
        tr_np_name_cl = "TrkData_cl_"+str(ooo)+"_pitch"+str(pitch)+".npz"  
        tr_np_name_cl = folder_path_npz+tr_np_name_cl
        print('tr_np_name ',tr_np_name)
        print('tr_np_name_true ',tr_np_name_true)
        print('tr_np_name_cl ',tr_np_name_cl)
        
    Track_info,Calo_info,gxcfib,gycfib,gzcfib,gxcpix,gycpix,gzcpix = functions.read_GeoFile(fgeo)
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
    # print(len(gxcfib[0][0]))
    # print ("x",gxcfib[0][0])
    # print ("y",gycfib[0][1])
    #print ("z",gzcfib[0])
    

    nview = int(Views)
    totLayers = int(Layers)
    nfiber  = int(Fibers)
    rfiber  = FibRadius #mm
    xoff    = rfiber#rfiber*0.5
    yoff    = rfiber#rfiber*0.5
    cx0strip = gxcfib[0][0][0] - rfiber + xoff
    cy0strip = gycfib[0][1][0] - rfiber + yoff
    
    stripNo = int((nfiber*rfiber+rfiber)/pitch)#128
    sigma_x = pitch/2#sqrt(12)#tobeasked
    print("[INFO]FIBER RADIUS:", rfiber*1000, " um")
    print("[INFO]Strip pitch:",pitch, " mm")
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

    # Define output treec
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

    # Informazioni sui cluster nelle viste orizzontale (x=0) e verticale(y=1) per i  layer
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
    truex=[]  
    predx=[]
    predy=[]
    truey=[]
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
        
        sResiduals3    = "Residuals3{0}".format(pos)
        sResiduals4    = "sResiduals3{0}[nlayers]/F".format(pos)
        exec("%s = ary.array('f',range(totLayers))" % sResiduals3)

        OutTree.Branch(sResiduals, eval(sResiduals),sResiduals2)
        OutTree.Branch(sResiduals3, eval(sResiduals3),sResiduals4)

    tree = tfile.Get("TrackerDigi")
    ptree = tfile.Get("Primary")
    nhit = tree.GetEntries()
    print('nhit: ',nhit)
    #nhit=10
    
    
    root_file.cd()
    
    #list for the dataframe
    raw_x=[]    #cluster,max cluster search lists
    raw_y=[]    #cluster,max cluster search lists
    raw_x_cl=[]   #cluster list
    raw_y_cl=[]   #cluster list
    raw_x_noclu=[]  #raw lists
    raw_y_noclu=[]  #raw lists
    true_x=[]   #list of true (primary track) values
    true_y=[]   #list of true (primary track) values

    #start for loop on hit
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
        #print('PRIMARY TRACK POSITION ', xmc,' ',ymc)
        cx = ptree.PrimaryParticleDirectionX
        cy = ptree.PrimaryParticleDirectionY
        cz = ptree.PrimaryParticleDirectionZ
       # print('PRIMARY TRACK direction ', cx,' ',cy,' ',cz)
        #reading of the DigiTrackerTree 
        tree.GetEntry(i) 
        PrimaryEnergy = (tree.PrimaryParticleEnergy)*1e3 #in keV
        dE_vec        = tree.Energy_keV
        view_vec      = tree.ViewNo
        layer_vec     = tree.LayerNo
        fiber_vec     = tree.FiberNo

        energy[0] = PrimaryEnergy
        
        #print('i ',i ,"Layers:",list(layer_vec),'len ', len(layer_vec))
        #print("View:",list(view_vec),'len ', len(view_vec))
        #print("Fibers:",list(fiber_vec), 'len ', len(fiber_vec))
        
    #list of objects per hit
        Strip_vec   = []
        PEstrip_vec = []
        Strips      = []
        Strip_r     = []
        Strips_cl   = []
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

        x0_1=[]
        x0_2=[]
        num_0=[]
        den_0=[]
        bc_noclmax=[]
        pe_noclmax=[]
        de1_noclmax=[]
        de2_noclmax=[]
        sumde2_noclmax=[]
        deb_noclmax=[]
        
        #ccc.append([])    
        
        for il in range(totLayers):
            
            Strips.append([])
            Strip_r.append([])
            Strips_cl.append([])
            eStrips.append([])
            PEStrips.append([])
            z.append([])
            String_zs.append([])

            Nclu.append([])
            Fstrip.append([])
            Clusiz.append([])
            Clunpe.append([])

            totPE.append([])

            x0_1.append([])
            x0_2.append([])
            num_0.append([])
            den_0.append([])
            bc_noclmax.append([])
            pe_noclmax.append([])
            de1_noclmax.append([])
            de2_noclmax.append([])
            sumde2_noclmax.append([])
            deb_noclmax.append([])
            
            for iv in range(nview):
                Strips[il].append([])
                Strip_r[il].append([])
                Strips_cl[il].append([])
                eStrips[il].append([])
                PEStrips[il].append([])
                z[il].append([])
                String_zs[il].append([])

                Nclu[il].append(0)
                Fstrip[il].append([])
                Clusiz[il].append([])
                Clunpe[il].append([])

                totPE[il].append(0)
        
                x0_1[il].append([])
                x0_2[il].append([])
                num_0[il].append([])
                den_0[il].append([])
                bc_noclmax[il].append([])
                pe_noclmax[il].append([])
                de1_noclmax[il].append([])
                de2_noclmax[il].append([])
                sumde2_noclmax[il].append([])
                deb_noclmax[il].append([])
                

        for j in range(len(fiber_vec)):
            
            mpe_fiber = (dE_vec[j]*fiberyield*trapeff*PDE)#numero medio di PE nella fibra

            if view_vec[j] == 0:
                cfib = gxcfib[layer_vec[j]][view_vec[j]][fiber_vec[j]]
                strip_centres, strip_index = functions.cstripCalc(rfiber,cx0strip,xoff,pitch,fiber_vec[j],gxcfib[layer_vec[j]][view_vec[j]][fiber_vec[j]],stripNo)
                
            if view_vec[j] == 1:
                cfib = gycfib[layer_vec[j]][view_vec[j]][fiber_vec[j]]
                strip_centres, strip_index = functions.cstripCalc(rfiber,cy0strip,yoff,pitch,fiber_vec[j],gycfib[layer_vec[j]][view_vec[j]][fiber_vec[j]],stripNo)
                
            Strip_vec.extend(strip_index)
            

            for istrip in range(len(strip_index)):
                Afract = round(functions.FracArea(rfiber,pitch,cfib,strip_centres[istrip]),1)
                mpe = dE_vec[j]*Afract*fiberyield*trapeff*PDE
                npe = ran.Poisson(mpe)#numero medio di PE nella strip
                PEstrip_vec.append(npe)###########################
            
            # print("#########fiber",j,":",fiber_vec[j], "Strip_vec",Strip_vec,' ' ,PEstrip_vec)

            '''for kk in range(len(strip_index)):
                Lay.append(layer_vec[j])
                View.append(view_vec[j])'''
            il = layer_vec[j]
            iv = view_vec[j]

            if (iv == 0): 
                strip_centres, strip_index = functions.cstripCalc(rfiber,cx0strip,xoff,pitch,fiber_vec[j],gxcfib[il][iv][fiber_vec[j]],stripNo)
                cfib=gxcfib[il][iv][fiber_vec[j]]
                
            
            if (iv == 1): 
                strip_centres, strip_index = functions.cstripCalc(rfiber,cy0strip,yoff,pitch,fiber_vec[j],gycfib[il][iv][fiber_vec[j]],stripNo)
                cfib=gycfib[il][iv][fiber_vec[j]]
                
            for istrip in range(len(strip_index)):
                Afract = round(functions.FracArea(rfiber,pitch,cfib,strip_centres[istrip]),1)
                mpe = dE_vec[j]*Afract*fiberyield*trapeff*PDE
                npe = ran.Poisson(mpe)
                if npe >= 3:#3
                    if strip_index[istrip] not in Strips[il][iv]:
                        Strips[il][iv].append(strip_index[istrip])
                        Strip_r[il][iv].append(strip_centres[istrip])
                        z[il][iv].append(gzcfib[il][iv][fiber_vec[j]])
                        eStrips[il][iv].append(dE_vec[j]*Afract)
                        PEStrips[il][iv].append(npe) 
                        #print('_______________lyn. ', il, ' view ', iv ,strip_centres[istrip])
                    else:
                        a=np.where(np.array(Strips[il][iv])==strip_index[istrip])[0][0]
                        eStrips[il][iv][a]+=(dE_vec[j]*Afract)
                        PEStrips[il][iv][a]+=npe
                        #print('_____^^^^^__________lyn. ', il, ' view ', iv ,strip_centres[istrip])
                    totPE[il][iv] += npe
        #hit for loop
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

                # if(iv==0):
                #     print("STRIP:",Strips[il][iv]) 
                #     print(ordered_strip[il][iv])
                #     print("CENTRE:",Strip_r[il][iv])
                #     print(ordered_centre[il][iv])
                #     print("ENERGY:",eStrips[il][iv])
                #     print(ordered_energy[il][iv])
                #     print("PE",PEStrips[il][iv])
                #     print(ordered_PE[il][iv])
                
        
        Strips   = ordered_strip
        Strip_r  = ordered_centre
        eStrips  = ordered_energy
        PEStrips = ordered_PE 
        #####fine ciclo sulle hit dell'evento i    
        # print("strip_r0x",Strip_r[0][0])
        # print("strip_r0y",Strip_r[0][1])
        # print("strip_r1x",Strip_r[1][0])
        # print("strip_r1y",Strip_r[1][1])
        # print("fine ciclo sulle hit dell'evento",i)
        #posizioni cluster prima della ricerca del baricentro 
        for il in range(totLayers):
            for iv in range(nview):
                if(iv==0):
                    raw_x_noclu.append([i,il,iv,Strip_r[il][iv],round(np.mean(z[il][iv]),3),PEStrips[il][iv]])
                    #print('clnocluev ',i,' ly ',il,' iv ',iv)
                    # print('****PE strip ',i,' il ',il,' iv ',iv ,' ' , PEStrips[il][iv],' ',Strip_r[il][iv])
                elif (iv==1):
                   # print('clnocluev ',i,' ly ',il,' iv ',iv)
                    raw_y_noclu.append([i,il,iv,Strip_r[il][iv],round(np.mean(z[il][iv]),3),PEStrips[il][iv]])
                    
                    #print('****PE strip ',i,' il ',il,' iv ',iv ,' ' , PEStrips[il][iv],' ',Strip_r[il][iv])

        
        vtrig= np.zeros((totLayers,nview))
        Tch[0] = (len(PEstrip_vec))
        
        ####CLUSTERING
        nclu,istrip0,nstrip0,charge0 = functions.Clustering(Strip_vec,PEstrip_vec)  
        # print('*****clustering******')
        # print('ev ',i, ' nclu ',nclu,' istrip0 ',istrip0,' nstrip0 ',nstrip0,' charge0 ',charge0)
        Nclu_tot += nclu

        for kk in range(len(istrip0)):
            Fstrip_tot.append(istrip0[kk])
            Clusiz_tot.append(nstrip0[kk])
            Clunpe_tot.append(charge0[kk])
                    
        for il in range(totLayers):
            for iv in range(nview):
                if iv == 0: 
                    v = "Hor"
                    # print('+++++++++++++ev ',i,' ly ',il, ' iv ',iv,' Stripr ',Strip_r[il][iv], ' Stripidx ',Strips[il][iv])
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

                nclu,istrip0,nstrip0,charge0 = functions.Clustering(Strips[il][iv],PEStrips[il][iv])
                #print('ev ',i,' ly ',il,' iv ',iv, ' test_clustering: istrip0  ',istrip0, ' nstrip0 ',nstrip0,' charge0 ',charge0 )
                
                Nclu[il][iv] = nclu
                for kk in range(len(istrip0)):
                    Fstrip[il][iv].append(istrip0[kk])
                    Clusiz[il][iv].append(nstrip0[kk])
                    Clunpe[il][iv].append(charge0[kk])
               # print("FSTRIP:", Fstrip)
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
        idx_cl=[]
        indices_cl=[]
        Strips_cl=[]
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
        for il in range(totLayers):
            istrip_vec.append([])
            idx_cl.append([])
            indices_cl.append([])
            Strips_cl.append([])
            bc.append([])
            debc.append([])
            primary_track.append([])
            for iv in range(nview):
                bc[il].append([])
                debc[il].append([])
                primary_track[il].append([])
                #print("view",iv)
                istrip_vec[il].append([])
                idx_cl[il].append([])
                indices_cl[il].append([])
                Strips_cl[il].append([])
                
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

                    if iv == 0: 
                        Layx.append(il)
                    elif iv == 1:
                        Layy.append(il)
                   # print('clsize',Clusiz[il][iv])
                    # print('fstrip ',Fstrip[il][iv][0] )
                    #for k in range(0,Clusiz[il][iv][0]):                        
                        #idx_cl[il][iv].append(Fstrip[il][iv][0]+k) 
                        #print(' ly ',il,' iv ', iv,' index cl',Fstrip[il][iv][0], 'new index ',Fstrip[il][iv][0]+k)
                    
                    for j in range(0,len(Clusiz[il][iv])):
                        idx_cl[il][iv].append([])       
                       # indices_cl[il][iv].append([]) 
                        #Strips_cl[il][iv].append([])
                        for k in range(0,Clusiz[il][iv][j]):                        
                            idx_cl[il][iv][j].append(Fstrip[il][iv][j]+k) 
                            #print('j ',j,' k ',k,' ly ',il,' iv ', iv,' Fstrip ',Fstrip[il][iv][j], 'new index ',Fstrip[il][iv][j]+k)
                    
                   
                    #print('ev ',i,' il ',il,' iv ',iv, 'idx_cl ',idx_cl[il][iv])
                                  
                    
                    for k in range(0,len(idx_cl[il][iv])):
                        tmp=[index for index, element in enumerate(Strips[il][iv]) if element in idx_cl[il][iv][k]]
                        #print('tmp_index ',tmp)
                        indices_cl[il][iv].append(tmp)
                  
                    # print('******** indices_cl ',indices_cl[il][iv] )
                    # print('Stripr ',Strip_r[il][iv] )
                    for scl in range(0,len(indices_cl[il][iv])):
                        x0_1[il][iv].append([indices_cl[il][iv][scl][0]])
                        x0_2[il][iv].append([indices_cl[il][iv][scl][-1]])  
                        idx1=indices_cl[il][iv][scl][0]
                        idx2=indices_cl[il][iv][scl][-1]+1                  
                        #print('____________')
                        #print(i,' il ',il,' iv ',iv,' idx1 ',indices_cl[il][iv][scl][0])
                        #print(i,' il ',il,' iv ',iv,'idx2 ',indices_cl[il][iv][scl][-1])
                        # print(x0_1[il][iv][scl])
                        # print(x0_2[il][iv][scl])
                        #print('____________')
                        #print('Selected strip r ',Strip_r[il][iv][indices_cl[il][iv][scl][0]:indices_cl[il][iv][scl][-1]+1])
                        #print('Selected PE ',PEStrips[il][iv][indices_cl[il][iv][scl][0]:indices_cl[il][iv][scl][-1]+1])
                        tmp_num=sum(np.multiply(Strip_r[il][iv][idx1:idx2],PEStrips[il][iv][idx1:idx2]))
                        
                        
                        #print('num ',tmp_num)
                        tmp_den=sum(PEStrips[il][iv][idx1:idx2])
                        #print('den ',tmp_den)
                        #àprint('[DEBUG fuori if ****] ev ',i,' il ',il,' iv ',iv,' num ',tmp_num,' tmp_den ', tmp_den,' tmp_num/tmp_den ',tmp_num/tmp_den)
                        if(tmp_den>=3):

                            de1_noclmax=np.array(PEStrips[il][iv][idx1:idx2])*sigma_x/tmp_den                          
                            sumde2_noclmax=np.array(Strip_r[il][iv][idx1:idx2])*tmp_den-np.array(Strip_r[il][iv][idx1:idx2])*np.array(PEStrips[il][iv][idx1:idx2])                                               
                            de2_noclmax=np.array(sumde2_noclmax)*pow(np.array(PEStrips[il][iv][idx1:idx2]),1/2)/pow(tmp_den,2)                            
                            #print('[DEBUG****] ev ',i,' il ',il,' iv ',iv,' num ',tmp_num,' ', tmp_den,' tmp_num/tmp_den ',tmp_num/tmp_den)
                            bc_noclmax[il][iv].append(round(tmp_num/tmp_den,3))
                            pe_noclmax[il][iv].append(round(tmp_den,3))
                            #breakpoint()
                            deb_noclmax[il][iv].append(round(pow(sum(pow(de1_noclmax,3))+sum(pow(de2_noclmax,2)),0.5),2))                            
                           
                            Strips_cl[il][iv].append(Strip_r[il][iv][indices_cl[il][iv][scl][0]:indices_cl[il][iv][scl][-1]+1]) 
                        
                    # for i_bc in range(0,len(bc_noclmax[il][iv])):
                    #     print('______Baricentro di carica_______')
                    #     print('ev ',i,' ly ',il,' iv ',iv,' bc ',bc_noclmax[il][iv][i_bc],'debc ',deb_noclmax[il][iv][i_bc])
                    
                    
                    
                    
                    String_zs[il][iv] = round(np.mean(z[il][iv]),3)
                    if iv==0:
                        
                        raw_x_cl.append([i,il,iv,bc_noclmax[il][iv], String_zs[il][iv],deb_noclmax[il][iv],pe_noclmax[il][iv]])
                        print('cl_ev NOMAX',i,' ly ',il,' iv ',iv,' ',bc_noclmax[il][iv],' ',String_zs[il][iv],' ',pe_noclmax[il][iv])
                    if iv==1:
                        
                        raw_y_cl.append([i,il,iv,bc_noclmax[il][iv], String_zs[il][iv],deb_noclmax[il][iv],pe_noclmax[il][iv]])
                        print('cl_ev NOMAX',i,' ly ',il,' iv ',iv,' ',bc_noclmax[il][iv],' ',String_zs[il][iv],' ',pe_noclmax[il][iv])
                    maxClu    = max(eval(sTnpe)[:eval(sTclu)[0]])
                    indexmax  = eval(sTnpe)[:eval(sTclu)[0]].index(maxClu)
                    fstripmax = eval(sFstr)[indexmax]

                    indexfstripmax    = eval(sstrip).index(fstripmax)
                    #print('sCsiz ',eval(sCsiz))
                    for istr in range(eval(sCsiz)[indexmax]):                  
                        istrip_vec[il][iv].append(indexfstripmax+istr)
                        
                    x0 = istrip_vec[il][iv][0]
                    x1 = istrip_vec[il][iv][-1]+1
                    #print('x1 ',x1)
                    num = sum(np.multiply(Strip_r[il][iv][x0:x1],PEStrips[il][iv][x0:x1]))
                    den = sum(PEStrips[il][iv][x0:x1])
                   
                    #print('den by lop ',den)
                    if den >= 3:
                        de1 = np.array(PEStrips[il][iv][x0:x1])*sigma_x/den
                        sumde2 = np.array(Strip_r[il][iv][x0:x1])*den-np.array(Strip_r[il][iv][x0:x1])*np.array(PEStrips[il][iv][x0:x1])
                        de2 = np.array(sumde2) *pow(np.array(PEStrips[il][iv][x0:x1]),1/2)/pow(den,2)
                        bc[il][iv]   = round((num/den),3)
                        
                        debc[il][iv] = (pow(sum(pow(de1,2))+sum(pow(de2,2)),0.5))
                        
                        if iv==0:
                            primary_track[il][iv] = (functions.evaluateTrackPoint(xmc,zmc,cx,cz,String_zs[il][iv],flag))
                            
                            x.append(bc[il][iv])
                            x_vec[il] = bc[il][iv]
                            dx.append(debc[il][iv])
                            dx_vec[il] = debc[il][iv]
                            zx.append(String_zs[il][iv])
                            zx_vec[il] = String_zs[il][iv]
                            dzx.append(2*rfiber)
                            dzx_vec[il] = 2*rfiber
                                                        
                            xtrack.append(primary_track[il][iv])
                            xMC_vec[il] = primary_track[il][iv]
                            #print('xMC_vec ',xMC_vec)
                            zxMC_vec[il] = String_zs[il][iv]
                            #print('###after clustering ev ', i, ' layers ', il, ' view ', iv, ' xrec ',
                                #bc[il][iv],' zx_rec ',String_zs[il][iv],' xMC_vec ',primary_track[il][iv], 'zxtrack ',String_zs[il][iv])
                            raw_x.append([i,il,iv, bc[il][iv],String_zs[il][iv],debc[il][iv],2*rfiber,round(xMC_vec[il],3),String_zs[il][iv]])
                            true_x.append([i,il,iv,round(xMC_vec[il],3)])
                            # print('xMC_vec[il] ',xMC_vec[il])
                            #print('primary_track_point x ', primary_track[il][iv],' xmc ',xmc,' zmc ',zmc, ' xmcvec ',xMC_vec[il],' z ', String_zs[il][iv])
                        elif iv==1:
                            primary_track[il][iv] = (functions.evaluateTrackPoint(ymc,zmc,cy,cz,String_zs[il][iv],flag))
                            
                            y.append(bc[il][iv])
                            y_vec[il] = bc[il][iv]
                            dy.append(debc[il][iv])
                            dy_vec = debc[il][iv]
                            zy.append(String_zs[il][iv])
                            zy_vec[il] = String_zs[il][iv]
                            dzy.append(2*rfiber)
                            dzy_vec[il] = 2*rfiber
                            ytrack.append(primary_track[il][iv])
                            yMC_vec[il] = primary_track[il][iv]
                            zyMC_vec[il] = String_zs[il][iv]
                            #print('primary_track_point y ', primary_track[il][iv],' ymc ',ymc,' zmc ',zmc,' yMC ',yMC_vec[il],'z ',String_zs[il][iv])
                            #print('###after clustering ev ', i, ' layers ', il, ' view ', iv, ' xrec ',
                                #bc[il][iv],' zx_rec ',String_zs[il][iv],' yMC_vec ',primary_track[il][iv], 'zxtrack ',String_zs[il][iv])
                            raw_y.append([i,il,iv, bc[il][iv],String_zs[il][iv],debc[il][iv],2*rfiber,round(yMC_vec[il],3),String_zs[il][iv]])
                            true_y.append([i,il,iv,round(yMC_vec[il],3)])
                            # print('ev ', i, ' layers ', il, ' view ', iv, ' yrec ',
                            #     bc[il][iv],' zy_rec ',String_zs[il][iv],' ytrack ',primary_track[il][iv], 'zytrack ',String_zs[il][iv])


        for ihit in range(len(Layx)):
            Ilayerx[ihit] = Layx[ihit]
        #print("layx",Layx)
        for ihit in range(len(Layy)):
            Ilayery[ihit] = Layy[ihit]                                      
        #print("layy",Layy)

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
        '''print ("x",x, dx)
        print ("y",y, dy)
        print ("zx",zx, dzx)
        print ("zy",zy, dzy)'''
        
        '''if np.shape(x)[0] < totLayers and np.shape(y)[0] < totLayers: 
            nlayers[0]  = 0
            trk_flag[0] = 0'''
        #elif np.shape(x)[0] == totLayers and np.shape(y)[0] == totLayers: 
        if np.shape(x)[0] == totLayers and np.shape(y)[0] == totLayers: 
            nlayers[0]  = totLayers
            nlayersx[0]  = totLayers
            nlayersy[0]  = totLayers
            trk_flag[0] = 1
            x_reco = TGraphErrors(totLayers,zx,x,dzx,dx)
            x_reco.SetTitle("Hor")
            x_true = TGraph(totLayers,zx,xtrack)  
            y_reco = TGraphErrors(totLayers,zy,y,dzy,dy)
            y_reco.SetTitle("Ver")
            y_true = TGraph(totLayers,zy,ytrack)  
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
                    x_reco.Draw("AP")
                    x_true.Draw("*")
                    # print("xreco", x, " MC x ", xtrack)
                if iv == 1:
                    xmin = min(zy)-0.05
                    xmax = max(zy)+0.05
                    rettay = TF1("rettay","[0]*x+[1]",xmin,xmax)
                    rettaMC = TF1("rettaMC","[0]*x+[1]",min(ytrack),max(ytrack))
                    # print('ytrack',ytrack)
                    y_reco.SetMarkerStyle(11)
                    #y_reco.GetXaxis().SetLimits(-1,1)
                    y_reco.GetYaxis().SetRangeUser(min(y),max(y))
                    y_reco.Fit("rettay","Q")
                    y_true.Fit("rettaMC","Q")
                    y_reco.Draw("AP")
                    y_true.Draw("*")
                #ccc[i][iv].Write() 
                #x_reco.Write()
                #y_reco.Write()
                #x_true.Write()
                #y_true.Write()
            xchi = rettax.GetChisquare()/rettax.GetNDF() 
            xprob = rettax.GetProb()
            mxrec = (rettax.GetParameter(0))
            qxrec = rettax.GetParameter(1)
            ychi = rettay.GetChisquare()/rettay.GetNDF()
            yprob = rettay.GetProb() 
            myrec = (rettay.GetParameter(0))
            qyrec = rettay.GetParameter(1)

            mx[0] = mxrec
            qx[0] = qxrec
            chi2x[0] = xchi
            my[0] = myrec
            qy[0] = qyrec
            chi2y[0] = ychi
            
            crFx,crFy,crFz = functions.evaluateRecoCosineDirector(mxrec, myrec)
            scalarF = functions.ScalarProduct(cx,cy,cz,crFx,crFy,crFz)   
            thetaF = TMath.ACos(scalarF)*180/TMath.Pi()
            
            cosx[0] = crFx
            cosy[0] = crFy
            cosz[0] = crFz

            cosxMC[0] = cx
            cosyMC[0] = cy
            coszMC[0] = cz
            
            angle_deg[0] = thetaF
            
            rx = functions.CalcResiduals(totLayers,x,zx,mxrec,qxrec)
            ry = functions.CalcResiduals(totLayers,y,zy,myrec,qyrec)
            r2x = functions.CalcResiduals(totLayers,xtrack,zx,mxrec,qxrec)
            r2y = functions.CalcResiduals(totLayers,ytrack,zy,myrec,qyrec)
            predx.append(np.array(functions.PredValue(totLayers,zx,mxrec,qxrec)))
            truex.append(xtrack)           
            predy.append(np.array(functions.PredValue(totLayers,zy,myrec,qyrec)))
            truey.append(ytrack)           
            for il in range(totLayers):
                for iv in range(nview):
                    if iv == 0: 
                        pos = "x"
                    else:
                        pos = "y"
                    sResiduals = "Residuals{0}".format(pos)
                    exec("%s[il] = r%s[il]" % (sResiduals,pos))
                    sResiduals3 = "Residuals3{0}".format(pos)
                    exec("%s[il] = r2%s[il]" % (sResiduals3,pos))
                    #print(r2_score(xtrack,functions.PredValue(totLayers,zx,mxrec,qxrec)))
                           
        
        else:
            nlayers[0] = 0
            nlayersx[0]  = np.shape(x)[0]
            nlayersy[0]  = np.shape(y)[0]

            trk_flag[0] = 0
        
        OutTree.Fill()
        
   # print('r2_x',r2_score(truex,predx))
    #print('r2_y',r2_score(truey,predy))
    # print(len(truex),' true ',truex,' ',predx)
      
    #__________________________________fine ciclo sugli eventi_________________________________________#
    #_________________________saving information in csv and npz files__________________________________#
    
    #csv file 
    data_trackingx_cl=pd.DataFrame(raw_x,columns=['Ev','ly','view','xrec','zxrec','dx','dxz','xMC','zMC'])
    data_trackingy_cl=pd.DataFrame(raw_y,columns=['Ev','ly','view','yrec','zyrec','dy','dyz','yMC','zMC'])
    
    data_trackingx_nocl=pd.DataFrame(raw_x_noclu,columns=['Ev','ly','view','x_hit','zx_hit','PE'])   
    data_trackingx_nocl['xMC']=data_trackingx_cl['xMC']
    
    data_trackingy_nocl=pd.DataFrame(raw_y_noclu,columns=['Ev','ly','view','y_hit','zy_hit','PE'])
    data_trackingy_nocl['yMC']=data_trackingy_cl['yMC']
    
    data_trackingx_cl_nomax=pd.DataFrame(raw_x_cl, columns=['Ev','ly','view','x_hit','zx_hit','dx_hit','PE'])
    data_trackingx_cl_nomax['xMC']=data_trackingx_cl['xMC']
    
    data_trackingy_cl_nomax=pd.DataFrame(raw_y_cl, columns=['Ev','ly','view','y_hit','zy_hit','dy_hit','PE'])
    data_trackingy_cl_nomax['yMC']=data_trackingy_cl['yMC']

    if(save_parquet==1):
        tab_y_cl_nomax=pa.Table.from_pandas(data_trackingy_cl_nomax)
        pq.write_table(tab_y_cl_nomax,OutputNameCSV_y_cl+'.parquet')
        tab_x_cl_nomax=pa.Table.from_pandas(data_trackingx_cl_nomax)
        pq.write_table(tab_x_cl_nomax,OutputNameCSV_x_cl+'.parquet')
    if(save_csv==1):
        #data_trackingx_cl.to_csv(OutputNameCSV_x)
        #data_trackingy_cl.to_csv(OutputNameCSV_y)   
        data_trackingx_nocl.to_csv(OutputNameCSV_x_nocl)
        data_trackingy_nocl.to_csv(OutputNameCSV_y_nocl)
        data_trackingx_cl_nomax.to_csv(OutputNameCSV_x_cl+".csv")
        data_trackingy_cl_nomax.to_csv(OutputNameCSV_y_cl+".csv")   
    
    #raw_y_noclu.append(data_trackingy_cl['yMC'].tolist())
    #raw_x_noclu.append(data_trackingx_cl['xMC'].tolist())
    
    data_array_x=np.array(raw_x_noclu,dtype=object)
    data_array_y=np.array(raw_y_noclu,dtype=object)
    data_array_truex=np.array(true_x)
    data_array_truey=np.array(true_y)
    
    data_array_cl_x=np.array(raw_x_cl,dtype=object)
    data_array_cl_y=np.array(raw_y_cl,dtype=object)
    
    if(save_npz==1):
        np.savez_compressed(tr_np_name,x_view=data_array_x,y_view=data_array_y)
        np.savez_compressed(tr_np_name_true,x_view=data_array_truex,y_view=data_array_truey)
        np.savez_compressed(tr_np_name_cl,x_view=data_array_cl_x,y_view=data_array_cl_y)

    tfile.Close()
    OutTree.Write()
    root_file.Write()
    root_file.Close()
