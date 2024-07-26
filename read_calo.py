import pylab

import math

import matplotlib as mpl

import scipy
import csv

from matplotlib import pyplot as plt
from matplotlib import gridspec
import matplotlib.colors as colors

import matplotlib.patches as patches
from matplotlib.collections import PatchCollection

from scipy import optimize

from matplotlib.axes import Axes
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import poisson
from scipy.optimize import curve_fit


import os
import sys
import time

#import pyfits
from ROOT import TFile, TTree, TChain, TH1D, TMath, TCanvas, TH1F, TH2F

import array as ary
import numpy as np

import glob


def read_GeoFile(fgeo):
    f = open(fgeo,'r')
    lines = len(open(fgeo).readlines())
    
    overall_track_geom = np.zeros(7) # n.of Layer,Views,Fibers, FibLength, FibRadius, TrackerLength, TrackerWidth
    overall_cal_geom   = np.zeros(7) # Pixel Size, Depth, nPix_per_view, TotPixels, CaloLength, CaloWidth
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
                nPix = int(aaa[3]) #numero pixel lungo X ed Y
                totPix = nPix*nPix*25
                CalLayer = int(aaa[4]) #number of layer along z axis
                CaloLength = float(aaa[5])
                CaloWidth  = float(aaa[6])
                #CaloLength = float(aaa[4])
                #CaloWidth  = float(aaa[5])
                overall_cal_geom[0] = pixSize
                overall_cal_geom[1] = pixDepth
                overall_cal_geom[2] = nPix   
                overall_cal_geom[3] = totPix 
                overall_cal_geom[4] = CaloLength
                overall_cal_geom[5] = CaloWidth
                print("nPix: ",nPix)    
                print("pixTot= ", totPix)
                gxcpix = np.zeros(totPix)
                
                gycpix = np.zeros(totPix)
                
                gzcpix = np.zeros(totPix)
            else:
                iXPix = int(aaa[1])
                iYPix = int(aaa[2])
                iZPix = int(aaa[3])
                iPix  = int(aaa[4])
                xcPix = float(aaa[5])
                ycPix = float(aaa[6])
                zcPix = float(aaa[7])

                gxcpix[iPix] = xcPix
                gycpix[iPix] = ycPix
                gzcpix[iPix] = zcPix
                
            

            ical +=1

    f.close()
    return overall_track_geom, overall_cal_geom,gxcfib,gycfib,gzcfib,gxcpix,gycpix,gzcpix

if __name__ == '__main__':

    froo = sys.argv[1]
    fgeo = sys.argv[2]

    froo += ".root"
    fgeo += ".txt"
    
    # print "[INFO] ROOT file",froo
    # print "[INFO] Geometry file",fgeo

    ### Load Geometry from geo file
    Track_info,Calo_info,gxcfib,gycfib,gzcfib,gxcpix,gycpix,gzcpix = read_GeoFile(fgeo)
    Layers        = Track_info[0]
    Views         = Track_info[1]
    Fibers        = Track_info[2]
    FibLength     = Track_info[3]
    FibRadius     = Track_info[4]
    TrackerLength = Track_info[5]
    TrackerWidth  = Track_info[6]
    # print ""
    # print "\t --- Tracker (DIM",TrackerLength,"x",TrackerLength,"x",TrackerWidth,"mm^3) ---"
    # print "[TRACKER_GEO_INFO]: Number of Layers in the Tracker = ", Layers
    # print "[TRACKER_GEO_INFO]: Number of Views per Layer = ", Views
    # print "[TRACKER_GEO_INFO]: Number of Fibers per View = ", Fibers
    # print "[TRACKER_GEO_INFO]: Fiber Length = ",FibLength,"mm","Fiber Radius = ",FibRadius,"mm"

    PixelSize   = Calo_info[0]
    PixelDepth  = Calo_info[1]
    nPix        = Calo_info[2]
    TotPix      = Calo_info[3]
    #CalLayer =    Calo_info[4]
    CaloLength  = Calo_info[4]
    CaloWidth   = Calo_info[5]
    
    #print ("")
    #print ("\t --- Calorimeter (DIM",CaloLength,"x",CaloLength,"x",CaloWidth,"mm^3) ---")
    #print ("[CALO_GEO_INFO]: Dimension = " ,nPix,"x",nPix, "x 25","TotalPixels= ", TotPix)
    #print ("[CALO_GEO_INFO]: Pixel Size = ", PixelSize,"mm")
    #print ("[CALO_GEO_INFO]: Pixel Depth = ", PixelDepth,"mm")
    #print ("")
    
    h0H = TH1F("h0H", "h0H", 100, 0, 100)
    #h1H = TH1F("h1H", "h1H", 600, 0, 600)
    #h2H = TH2F("h2", "h2", 200, 0, 200, 600, 0, 600)
    #plus_minus = "±"
    

    ### Load DATA from root file
    fr = TFile(froo)
    tpr=fr.Get("Primary")
    tkr=fr.Get("TrackerDigi")
    #cal=fr.Get("CalorimeterHit")
    cal=fr.Get("CalorimeterDigi")
    nevts = tkr.GetEntries()
    #with open('enePix_gamma.csv', mode='a', newline='') as file:
               #writer = csv.writer(file)
    ### Read Data
    
    ngood = 0
    #dist_weigthed =[]
    #R = {}
    #posizioni_laterali = []
    energyZ = {}
    #layer_z_array=[]
    z_hit = []
    x_hit = []
    y_hit = []
    ene_hit = []
    R_hit = []
    radii = []
    primary_trace_x = []
    primary_trace_y = []
    primary_trace_z = []
    for i in range(1):        
    #for i in range(nevts):
        try:
            tpr.GetEntry(i)
            tkr.GetEntry(i)
            cal.GetEntry(i)

            #--- Primary Track ---#
            pdg = tpr.PrimaryParticlePDG
            primary_en = tpr.PrimaryParticleEnergy
            x0 = tpr.PrimaryParticlePositionX
            y0 = tpr.PrimaryParticlePositionY
            z0 = tpr.PrimaryParticlePositionZ
            cx = tpr.PrimaryParticleDirectionX
            cy = tpr.PrimaryParticleDirectionY
            cz = tpr.PrimaryParticleDirectionZ
            
            primID  = tpr.eventID
            trackID = tkr.eventID
            calID   = cal.eventID
            z_hit.append([])
            x_hit.append([])
            y_hit.append([])
            ene_hit.append([])
            R_hit.append([])
            #nt = 1000
            #xp = np.zeros(nt)
            #yp = np.zeros(nt)
            #zp = np.zeros(nt)

                # -- Create the primary track until it interact with the detector -- #
            #nnt = 0
            #for k in range(nt):
             #   tt = float(k)*10.   
                #zz = z0 + cz*tt
                    #if(zz<ztopxz or zz<ztopyz):
                    #    break
                #xp[k] = x0 + cx*tt
                #yp[k] = y0 + cy*tt
                #zp[k] = z0 + cz*tt
               # t = (zp[k]-z0)/cz
            #nnt = k
        
            #--- Calorimeter ---#
            
            #hitpix = int(cal.CalPixels_Hit)
            hitpix = int(cal.DigiCalPixels_Hit)
            cpix_x = np.empty(hitpix, dtype=float)
            cpix_y = np.empty(hitpix, dtype=float)
            cpix_z = np.empty(hitpix, dtype=float)
            epix   = np.zeros(hitpix, dtype=float)
            layer_z = np.empty(hitpix, dtype=int)
            layer_x  = np.empty(hitpix, dtype=int)
            layer_y  = np.empty(hitpix, dtype=int)
            R = np.empty(hitpix, dtype=float)
            
            #w = np.empty(hitpix, dtype=float)
            #patches_pix =  []
            nt = 50
            xp = np.zeros(nt)
            yp = np.zeros(nt)
            zp = np.zeros(nt)

                # -- Create the primary track until it interact with the detector -- #
            nnt = 0
            for k in range(nt):
                    tt = float(k)*10.
                    zz = z0 + cz*tt
                    #if(zz<ztopxz or zz<ztopyz):
                    #    break
                    xp[k] = x0 + cx*tt
                    yp[k] = y0 + cy*tt
                    zp[k] = z0 + cz*tt

            nnt = k
            for k in range(nnt):
                 primary_trace_x.append(xp[k])
                 primary_trace_y.append(yp[k])
                 primary_trace_z.append(zp[k])
            jjcal=0
            for jcal in range(hitpix):
                    
                    calUnitNo = cal.CalUnitNo[jcal]
                    calPixX   = cal.CalPixelX[jcal]
                    calPixY   = cal.CalPixelY[jcal]
                    calPixZ   = cal.CalPixelZ[jcal]
                    enePix    = cal.CalEnergy[jcal]
                    
                    #w = [cal.CalEnergy[jcal] / sum(cal.CalEnergy[jcal])]
                    #w = np.array(w)
    
                    #print("pesi= ", w)
                    if(enePix<=0):
                        continue

                    epix[jjcal] = enePix
                  
                    cpix_x[jjcal] = gxcpix[calUnitNo]
                    cpix_y[jjcal] = gycpix[calUnitNo]
                    cpix_z[jjcal] = gzcpix[calUnitNo]
                    
                    #z_hit[i].append(cpix_z[jjcal]/PixelDepth)  
                    #x_hit[i].append(cpix_x[jjcal]/PixelDepth) 
                    #y_hit[i].append(cpix_y[jjcal]/PixelDepth) 

                    z_hit[i].append(cpix_z[jjcal])  
                    x_hit[i].append(cpix_x[jjcal]) 
                    y_hit[i].append(cpix_y[jjcal]) 

                    ene_hit[i].append(epix[jjcal])
                    if calPixZ  in energyZ:
                     energyZ[calPixZ] += enePix
                    else:
                     energyZ[calPixZ] = enePix
                     
                     
                   
                    #w_z=[w[j] for j in range(len(w)) if cal.CalPixelZ[jcal] == calPixZ]

                    #w[jjcal] = (epix[jjcal]/ sum(cal.CalEnergy))
                    layer_z[jjcal] = int(cpix_z[jjcal]/PixelDepth)
                    #layer_z_array = np.array(layer_z[jjcal])
                    layer_x[jjcal] = int(cpix_x[jjcal]/PixelDepth)
                    layer_y[jjcal] = int(cpix_y[jjcal]/PixelDepth)
                    #Z_new = list(energyZ.keys())
                    #Ez_new = list(energyZ.values())
                           
                   
                   
                    
                    jjcal+=1

            
      
            #fig = plt.figure(figsize=(26,13))
            #plt.scatter(layer_z_array, R)
            #fig = plt.figure(figsize=(10, 7))
            #ax = fig.add_subplot(111, projection='3d') 
                  
           # gspec_right = gridspec.GridSpec(ncols=3, nrows=1, width_ratios=[2, 2, 2])
            #gspec_right = gridspec.GridSpec(ncols=2, nrows=1, width_ratios=[2,2])
            #ax4 = fig.add_subplot(gspec_right[:, 0])
           # ax1 = fig.add_subplot(gspec_right[:, 0])
           # ax2 = fig.add_subplot(gspec_right[:, 1])
            #ax3 = fig.add_subplot(gspec_right[:, 2])
           # ax5 = fig.add_subplot(gspec_right[:, 2], projection='3d')  
            #ax5 = fig.add_subplot( projection='3d') 
            
                
           # scatter = ax1.scatter(x_hit, z_hit, c=ene_hit, cmap='rainbow',norm=mpl.colors.LogNorm())#, s=epix)
           # scatter = ax2.scatter(y_hit, z_hit, c=ene_hit, cmap='rainbow',norm=mpl.colors.LogNorm())#, s=epix)   
            #scatter = ax3.scatter(cpix_x, cpix_y, c=epix, cmap='rainbow',norm=mpl.colors.LogNorm()) 
            #scatter1 = ax1.scatter(cpix_x, cpix_z, c=epix,  cmap='rainbow', label='Calorimeter energy')
            #scatter = ax2.scatter(cpix_y, cpix_z, c=epix,  cmap='rainbow', label='Calorimeter energy')
            #scatter = plt.scatter(cpix_z, R_w)
            #scatter = plt.scatter(cpix_z, R, c='red')
            
            #scatter1 = ax1.scatter(layer_x, layer_z, c=epix,  cmap='rainbow', label='Calorimeter energy')
            #scatter = ax2.scatter(layer_y, layer_z, c=epix,  cmap='rainbow', label='Calorimeter energy')
            #scatter3d = ax5.scatter(x_hit, y_hit, z_hit, c=ene_hit, cmap='rainbow',norm=mpl.colors.LogNorm())
           
            
            #ax.set_xlim([-(nPix*1.5),(nPix*1.5)]) 
            #ax.set_ylim([-(nPix*1.5),(nPix*1.5)]) 
            #ax.set_zlim([-(nPix*1.5),(nPix*1.5)]) 
            #ax1.set_xlim(+25,-25) 
            #ax2.set_xlim(+25,-25) 
            #ax1.set_ylim(+25,-25) 
            #ax2.set_ylim(+25,-25) 
           
            #ax4.set_xlim([-(nPix*2)+1,(nPix*2)+1.]) 
            #ax4.set_ylim([-(nPix*2)+1,(nPix*2)+1.]) 
            #plt.draw()
                    #plt.pause(0.1)
                    #plt.show()
                    #time.sleep(2)
            #plt.show()
           

             # Calcola il centro di ogni bin
            
            #plt.savefig("./DisplayPlots_gamma/PDG"+str(i)+".png")
            #
            
            #h0H.Fill(enePix)
            

            
                            
        except KeyboardInterrupt as e:
            #print ("Exit")
            raise
   
    #fig = plt.figure(figsize=(30, 15))   
    
    fig = plt.figure(figsize=(26,13))    
    #gspec_right = gridspec.GridSpec(ncols=2, nrows=1, width_ratios=[2, 2])       
    ##ax1 = fig.add_subplot(gspec_right[:, 0])
    #ax2 = fig.add_subplot(gspec_right[:, 1])
    #ax5 = fig.add_subplot(gspec_right[:, 2], projection='3d')
    ax5 = fig.add_subplot(111, projection='3d')  
    
    for i in range(1):
    #    fig = plt.figure(figsize=(30, 15))   
        #fig, axs = plt.subplots(1, 2, figsize=(26, 13))
    #    axs = fig.subplots(1, 3)
        #ax = fig.add_subplot(111, projection='3d') 
        #scatter = ax1.scatter(x_hit[i], z_hit[i], c=ene_hit[i], cmap='rainbow',norm=mpl.colors.LogNorm())
        
        #ax1.set_title(f'Event {i}_x-z view')
        #ax1.set_xlabel('X')
        #ax1.set_ylabel('Z')
        #ax1.set_title('2D Scatter Plot x-z')

        
        #scatter = ax2.scatter(y_hit[i], z_hit[i], c=ene_hit[i], cmap='rainbow',norm=mpl.colors.LogNorm())
        #ax2.set_title(f'Event {i}_y-z view')
        #ax2.set_xlabel('Y')
        #ax2.set_ylabel('Z')
        #ax2.set_title('2D Scatter Plot y-z')
        
        scatter3d = ax5.scatter(x_hit[i], y_hit[i], z_hit[i], c=ene_hit[i], cmap='rainbow')#,norm=mpl.colors.LogNorm())
        scatter3d = ax5.scatter(primary_trace_x, primary_trace_y, primary_trace_z, c='blue')
        colorbar= plt.colorbar(scatter3d, ax=ax5, label='Calorimeter energy')
       

        ax5.set_title(f'Event {i}_#3D view')

        ax5.set_xlabel('X (mm)') 
        ax5.set_ylabel('Y (mm)')  
        ax5.set_zlabel('Z (mm)')  
        
        #ax5.set_xlabel('LayerX') 
        #ax5.set_ylabel('LayerY')  
        #ax5.set_zlabel('LayerZ')  
        
       
        #colorbar = plt.colorbar(scatter3d, ax=ax5, label=f'Calorimeter energy_{i}')
        
    # Prima subplot (x-z view)
        #axs[0].scatter(x_hit[i], z_hit[i], c=ene_hit[i], cmap='rainbow',norm=mpl.colors.LogNorm())
        #axs[0].scatter(yp[:nnt], zp[:nnt], marker = "*", color="red")
        
        #axs[0].scatter(yp[:nnt], zp[:nnt], marker = "*", color="red")
       
        
        #axs[0].set_title(f'Event {i}_y-z view')
        
        ##axs[0].set_ylabel('x [mm] ')
        #axs[0].set_xlabel('z [mm]')
    # Seconda subplot (y-z view)
        
        #axs[1].scatter(primary_trace_y, primary_trace_z)
        #axs[1].scatter(y_hit[i], z_hit[i], c=ene_hit[i], cmap='rainbow',norm=mpl.colors.LogNorm())
        
        #axs[1].set_title(f'Event {i}_y-z view')
        #axs[1].set_ylim([-(nPix*1.5),(nPix*1.5)]) 
       # axs[1].set_zlim([-(nPix*1.5),(nPix*1.5)]) 
        #axs[1].set_ylabel('y [mm]')
        #axs[1].set_xlabel('z [mm]')

     # Prima subplot (3dview)
        #axs[2] = fig.add_subplot(1, 3, 3, projection='3d')
        #axs[2].scatter(x_hit[i], y_hit[i], z_hit[i], c=ene_hit[i], cmap='rainbow', marker=".")
        
    #    axs[2].set_title(f'Event {i}_#3D view')
    #    axs[2].set_xlabel('x [mm]')
    #    axs[2].set_ylabel('y [mm] ')
    #    axs[2].set_zlabel('z [mm] ')

    
    plt.savefig("./DisplayPlots/PDG_proton_inclined"+str(i)+".png")
      
   