import pylab

import math

import matplotlib as mpl

from matplotlib import pyplot as plt
from matplotlib import gridspec
import matplotlib.colors as colors

import matplotlib.patches as patches
from matplotlib.collections import PatchCollection

from matplotlib.axes import Axes
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from mpl_toolkits.mplot3d import Axes3D

import os
import sys
import time

#import pyfits
from ROOT import TFile, TTree, TChain, TH1D, TMath, TCanvas, TH1F

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
    
    print ("")
    print ("\t --- Calorimeter (DIM",CaloLength,"x",CaloLength,"x",CaloWidth,"mm^3) ---")
    print ("[CALO_GEO_INFO]: Dimension = " ,nPix,"x",nPix, "x 25","TotalPixels= ", TotPix)
    print ("[CALO_GEO_INFO]: Pixel Size = ", PixelSize,"mm")
    print ("[CALO_GEO_INFO]: Pixel Depth = ", PixelDepth,"mm")
    print ("")
    
    h0H = TH1F("h0H", "h0H", 100, 0, 100)
    ### Load DATA from root file
    fr = TFile(froo)
    tpr=fr.Get("Primary")
    tkr=fr.Get("TrackerDigi")
    cal=fr.Get("CalorimeterHit")
    nevts = tkr.GetEntries()
    
    ### Read Data
    ngood = 0
    #for i in range(10):
    for i in range(nevts):
        try:
            tpr.GetEntry(i)
            tkr.GetEntry(i)
            cal.GetEntry(i)

            #--- Primary Track ---#
            pdg = tpr.PrimaryParticlePDG
            x0 = tpr.PrimaryParticlePositionX
            y0 = tpr.PrimaryParticlePositionY
            z0 = tpr.PrimaryParticlePositionZ
            cx = tpr.PrimaryParticleDirectionX
            cy = tpr.PrimaryParticleDirectionY
            cz = tpr.PrimaryParticleDirectionZ
            
            primID  = tpr.eventID
            trackID = tkr.eventID
            calID   = cal.eventID
            
            nt = 1000
            xp = np.zeros(nt)
            yp = np.zeros(nt)
            zp = np.zeros(nt)

                # -- Create the primary track until it interact with the detector -- #
            #nnt = 0
            #for k in range(nt):
            #    tt = float(k)*10.   
                #zz = z0 + cz*tt
            #        #if(zz<ztopxz or zz<ztopyz):
                    #    break
            #    xp[k] = x0 + cx*tt
            #    yp[k] = y0 + cy*tt
            #    zp[k] = z0 + cz*tt
            #    t = (zp[k]-z0)/cz
            #nnt = k
           
            #--- Calorimeter ---#
            
            hitpix = int(cal.CalPixels_Hit)
            #print("hitpix= ", hitpix)
            #cpix_x = np.empty(hitpix, dtype=np.float)
            #cpix_y = np.empty(hitpix, dtype=np.float)
            #cpix_z = np.empty(hitpix, dtype=np.float)
            #epix   = np.zeros(hitpix, dtype=np.float)
            cpix_x = np.empty(hitpix, dtype=float)
            cpix_y = np.empty(hitpix, dtype=float)
            cpix_z = np.empty(hitpix, dtype=float)
            epix   = np.zeros(hitpix, dtype=float)
            layer_z  = np.empty(hitpix, dtype=int)
            layer_x  = np.empty(hitpix, dtype=int)
            layer_y  = np.empty(hitpix, dtype=int)
            #patches_pix =  []
            
            jjcal=0
            for jcal in range(hitpix):
                    
                    calUnitNo = cal.CalUnitNo[jcal]
                    calPixX   = cal.CalPixelX[jcal]
                    calPixY   = cal.CalPixelY[jcal]
                    calPixZ   = cal.CalPixelZ[jcal]
                    enePix    = cal.CalEnergy[jcal]
                    

                    if(enePix<=0):
                        continue

                    epix[jjcal] = enePix
                    cpix_x[jjcal] = gxcpix[calUnitNo]
                    cpix_y[jjcal] = gycpix[calUnitNo]
                    cpix_z[jjcal] = gzcpix[calUnitNo]
                    #print("cpixz= ", cpix_z[jjcal])
                    #if(cpix_x[jjcal]>200):
                    #    print(i,' jcal ',jcal,' jjcal ',jjcal,' epix[jjcal] ',epix[jjcal],' x ',cpix_x[jjcal],' y ',cpix_y[jjcal], ' z ',cpix_z[jjcal] )
                    
                    layer_z[jjcal] = int(cpix_z[jjcal]/PixelDepth)
                    layer_x[jjcal] = int(cpix_x[jjcal]/PixelDepth)
                    layer_y[jjcal] = int(cpix_y[jjcal]/PixelDepth)
                    #print("layer: ", layer)
                    
                    #h0H.Fill(enePix)
                    
                    #ax2 = plt.hist(epix, bins=100, color='blue', alpha=0.7) 
                    jjcal+=1
            
            
            #h0H.Fill(enePix)      
            #print("hit= ", hitpix)
            #fig = plt.figure(figsize=(24,12))
            #gspec_right = gridspec.GridSpec(ncols=3, nrows=1, width_ratios=[1.5, 1.5, 2])
            #ax1 = fig.add_subplot(gspec_right[:, 0])
            #ax2 = fig.add_subplot(gspec_right[:, 1])
            #ax = fig.add_subplot(gspec_right[:, 2])
            #ax5 = fig.add_subplot(gspec_right[:, 2], projection='3d')  
           
            #scatter = ax1.scatter(cpix_x, cpix_z, c=epix, cmap='rainbow')#, s=epix)
            #scatter = ax2.scatter(cpix_y, cpix_z, c=epix, cmap='rainbow')#, s=epix)          
            #scatter3d = ax5.scatter(cpix_x, cpix_y, layer_z, c=epix, marker=".", s=10, norm=mpl.colors.LogNorm())
            #scatter3d = ax5.scatter(layer_x, layer_y, layer_z, c=epix, marker=".", s=9, norm=mpl.colors.LogNorm())

            #scatter3d= ax5.scatter(xp, yp, zp, marker=".")
            #ax5.plot(xp[:nnt], zp[:nnt], c='red')#"k--")
           
            
            # Additional settings or labels for the plot
            #ax1.set_xlabel('X')
            #ax1.set_ylabel('Z')
            #ax1.set_title('2D Scatter Plot x-z')

            #ax2.set_xlabel('Y')
            #ax2.set_ylabel('Z')
            #ax2.set_title('2D Scatter Plot y-z')

            #Add a colorbar to show the correspondence between color and layer number
            #colorbar = plt.colorbar(scatter, ax=ax1, label='Calorimeter energy')
            #colorbar = plt.colorbar(scatter, ax=ax2, label='Calorimeter energy')
            #colorbar = plt.colorbar(scatter, ax=ax5, label='Calorimeter energy')


            #ax5 = fig.add_subplot(111, projection='3d')  
            #ax5.set_xlabel('X Layer') 
            #ax5.set_ylabel('Y Layer')  
            #ax5.set_zlabel('Z Layer')  
            
            #ax5.set_xlim([-(nPix*2)+1,(nPix*2)+1.]) 
            #x5.set_ylim([-(nPix*2)+1,(nPix*2)+1.]) 
            #ax5.set_zlim([-(nPix*2)+1,(nPix*2)+1.]) 
            #ax5.set_zlim([-(CaloLength*1.5)+1,0])
        
            #ax1.set_xlim(-80,+80) 
            #ax2.set_xlim(-80,+80) 
            #for spine in ax.spines.values():
              #spine.set_visible(False)
            #ax5.set_zlim([0.,CalLayer])  
            
            #plt.draw()
                    #plt.pause(0.1)
                    #plt.show()
                    #time.sleep(2)
            #plt.show()
            #plt.savefig("./DisplayPlots/PDG"+str(i)+".png")
            #plt.savefig("./DisplayPlots_proton/PDG"+str(i)+".png")
            h0H.Fill(enePix)
                                
                                
        except KeyboardInterrupt as e:
            #print ("Exit")
            raise
    c=TCanvas('c','c')
    c.cd()
    h0H.Draw()  
    c.SaveAs("energy_proton.png")
                  
    

    
        
        
