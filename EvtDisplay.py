###come argomento Rootfile e Geofile senza estensioni

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

import os
import sys
import time

#import pyfits
from ROOT import TFile, TTree, TChain, TH1D, TMath, TCanvas, TH1F
#import ROOT

import array as ary
import numpy as np

import glob

############# IMPORTANTE #############
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

    overall_track_geom = np.zeros(8) # n.of Layer,Views,Fibers, FibLength, FibRadius, TrackerLength, TrackerWidth, GapBetweenLayers
    overall_cal_geom   = np.zeros(6) # Pixel Size, Depth, nPix_per_view, TotPixels, CaloLength, CaloWidth
    itr = 0
    ical = 0
    for i in range(lines):
        line = f.readline()
        aaa = line.split()
        Detector = str(aaa[0])
        if(Detector == "TRACKER"):
            if(itr==0): #primo rigo del file
                layer  = int(aaa[1])
                view   = int(aaa[2])
                fiber  = int(aaa[3])
                lfib   = float(aaa[4])
                rfib   = float(aaa[5])
                TrackerLength = float(aaa[6])
                TrackerWidth  = float(aaa[7])
                GapLayers = float(aaa[8])

                overall_track_geom[0] = layer
                overall_track_geom[1] = view
                overall_track_geom[2] = fiber
                overall_track_geom[3] = lfib
                overall_track_geom[4] = rfib
                overall_track_geom[5] = TrackerLength
                overall_track_geom[6] = TrackerWidth
                overall_track_geom[7] = GapLayers

                gxcfib = np.zeros((layer,view,fiber))
                gycfib = np.zeros((layer,view,fiber))
                gzcfib = np.zeros((layer,view,fiber))
            else:      #tutti gli altri righi
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

def init_plots(TrackerLength,TrackerWidth,CaloLength):    
    
    fig = plt.figure(figsize=(15,7))

    # create a 2-row 2-column container as the left container
    gspec_left = gridspec.GridSpec(ncols=3, nrows=2, width_ratios=[3.4,4.2,8])    

    ax2 = fig.add_subplot(gspec_left[:, 0])
    ax3 = fig.add_subplot(gspec_left[:, 1])   

    # now the plots are on top of each other, we'll have to adjust their edges so that they won't overlap
    gspec_left.update(wspace=0.15)

    # create a 1-row 1-column grid as the right container
    gspec_right = gridspec.GridSpec(ncols=2, nrows=1, width_ratios=[3,5])

    ax4 = fig.add_subplot(gspec_right[0, 1])

    # also, we want to get rid of the horizontal spacing in the left gridspec
    gspec_right.update(wspace=1.5, top = 0.65, bottom = 0.3)
    
    '''gspec = gridspec.GridSpec(ncols=3, nrows=3, width_ratios=[1,1,2], height_ratios=[0.5,1,0.5])
    
    ax2 = fig.add_subplot(gspec[:, 0])
    ax3 = fig.add_subplot(gspec[:, 1])
    ax4 = fig.add_subplot(gspec[1, 2])
    
    move_figure(fig,20,20)
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1,wspace=0.8)'''
    
    plt.ion()
    plt.show()
    
    ax2.set_xlim([-(TrackerLength*0.5)+1.,+(TrackerLength*0.5)+1.])
    ax2.set_ylim([-(TrackerWidth*0.5) +1.,+(TrackerWidth*0.5)+1. ])
    plt.draw()
    
    ax3.set_xlim([-(TrackerLength*0.5)+1.,+(TrackerLength*0.5)+1.])
    ax3.set_ylim([-(TrackerWidth*0.5) +1.,+(TrackerWidth*0.5)+1. ])
    plt.draw()
    
    ax4.set_xlim([-(CaloLength*0.5)+1.,+(CaloLength*0.5)+1.])
    ax4.set_ylim([-(CaloLength*0.5)+1.,+(CaloLength*0.5)+1.])
    plt.pause(0.01)
    plt.draw()

    axes = [ax2,ax3,ax4]
    return fig, axes

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
    GapLayers     = Track_info[7]
    # print ""
    # print "\t --- Tracker (DIM",TrackerLength,"x",TrackerLength,"x",TrackerWidth,"mm^3) ---"
    # print "[TRACKER_GEO_INFO]: Number of Layers in the Tracker = ", Layers
    # print "[TRACKER_GEO_INFO]: Number of Views per Layer = ", Views
    # print "[TRACKER_GEO_INFO]: Number of Fibers per View = ", Fibers
    # print "[TRACKER_GEO_INFO]: Fiber Length = ",FibLength,"mm","Fiber Radius = ",FibRadius,"mm"
    # print "[TRACKER_GEO_INFO]: Fiber LengthGap between layers = ", GapLayers,"mm"

    PixelSize   = Calo_info[0]
    PixelDepth  = Calo_info[1]
    nPix        = Calo_info[2]
    TotPix      = Calo_info[3]
    CaloLength  = Calo_info[4]
    CaloWidth   = Calo_info[5]
    # print ""
    # print "\t --- Calorimeter (DIM",CaloLength,"x",CaloLength,"x",CaloWidth,"mm^3) ---"
    # print "[CALO_GEO_INFO]: Dimension = ",nPix,"x",nPix," TotalPixels = ",TotPix
    # print "[CALO_GEO_INFO]: Pixel Size = ", PixelSize,"mm"
    # print "[CALO_GEO_INFO]: Pixel Depth = ", PixelDepth,"mm"
    # print "" 
    
    rfib     = FibRadius*2  # if is magnify by 10

    ### Prepare EvtDisplay
    cmap = mpl.cm.jet
    cmap.set_under('w') #Set the color for low out-of-range values
    cmap.set_bad('w') #Set the color for masked values

    fig,axes = init_plots(TrackerLength,TrackerWidth,CaloLength)
    ax2 = axes[0]
    ax3 = axes[1]
    ax4 = axes[2]

    #Scale of ColorBar
    emin_track = 0.     #keV
    emax_track = 100.  #keV

    emin_cal = 0.       #MeV
    emax_cal = 100.     #MeV

    ### Load DATA from root file
    fr = TFile(froo)
    tpr=fr.Get("Primary")
    tkr=fr.Get("TrackerDigi")
    cal=fr.Get("CalorimeterDigi")
    
    nevts = tkr.GetEntries()    
    # print "[DATA_INFO]: Number of events = ", nevts
    
    ### Read Data
    ngood = 0
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

            hitfib = int(tkr.DigiFibers_Hit)

            efib0 = np.zeros(hitfib, dtype=np.float)
            efib1 = np.zeros(hitfib, dtype=np.float)

            cfib_0_x = np.empty(hitfib, dtype=np.float)
            cfib_0_z = np.empty(hitfib, dtype=np.float)
            cfib_1_y = np.empty(hitfib, dtype=np.float)
            cfib_1_z = np.empty(hitfib, dtype=np.float)

            hitfib0 = 0
            hitfib1 = 0
            goodfib = 0

            ztopxz = -9999.
            ztopyz = -9999.

            primID  = tpr.eventID
            trackID = tkr.eventID
            calID   = cal.eventID
            
            #--- Tracker ---#
            for j in range(hitfib):
                layer = tkr.DigiLayerNo[j]
                view = tkr.DigiViewNo[j]
                fiber = tkr.DigiFiberNo[j]
                enefib = tkr.DigiEnergy_keV[j]

                if(enefib<=0):
                    continue
                
                goodfib += 1
                if(view==1):
                    efib0[hitfib0] = enefib
                    cfib_0_x[hitfib0] = gxcfib[layer][view][fiber]
                    cfib_0_z[hitfib0] = gzcfib[layer][view][fiber]
                    if(cfib_0_z[hitfib0]>ztopxz):
                        ztopxz = cfib_0_z[hitfib0]
                    hitfib0 += 1
                else:
                    efib1[hitfib1] = enefib
                    cfib_1_y[hitfib1] = gycfib[layer][view][fiber]
                    cfib_1_z[hitfib1] = gzcfib[layer][view][fiber]
                    if(cfib_1_z[hitfib1]>ztopyz):
                        ztopyz = cfib_1_z[hitfib1]
                    hitfib1 += 1

            if(goodfib>0):
                # print "Showing_",ngood,"(ievt = ",i,") -> fibersHit = ",goodfib
                # print "**[DEBUG]: Prim_ID",primID,"Track_ID",trackID,"Cal_ID",calID
                ngood += 1
                nt = 1000
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

                # -- Create the Patch collection for the Hits: trattini colorati nel grafico a seconda dell'energia depositata -- #
                patches_fib_xz =  []
                patches_fib_yz =  []

                ix = hitfib0
                iy = hitfib1
                for j in range(ix):
                    xst = cfib_0_x[j]
                    zst = cfib_0_z[j]
                    patches_fib_xz.append(patches.Circle((xst,zst), rfib))
                for j in range(iy):
                    yst = cfib_1_y[j]
                    zst = cfib_1_z[j]
                    patches_fib_yz.append(patches.Circle((yst,zst), rfib))

                ax2.cla() #clears an axes, i.e. the currently active axes in the current figure. It leaves the other axes untouched
                ax2.set_xlabel('X (mm)')
                ax2.set_ylabel('Z (mm)')
                ax2.set_xlim([-(TrackerLength*0.5)-1.,+(TrackerLength*0.5)+1.])
                #ax2.set_xlim(-1,1)
                ax2.set_ylim([-(TrackerWidth*0.5) -1.,+(TrackerWidth*0.5)+1. ])
                ax2.set_title("Tracker XZ View")
                ax2.plot(xp[:nnt], zp[:nnt], ":", color="darkgrey")
                #ax2.grid(axis = 'y')

                
                pafib_xz = PatchCollection(patches_fib_xz, cmap=cmap, alpha=0.9, lw=0., norm=colors.Normalize(vmin=emin_track, vmax=emax_track)) #LogNorm in logaritmic
                pafib_xz.set_array(efib0)
                ax2.add_collection(pafib_xz)

                '''try:
                    cbar_xz.remove()
                except NameError:
                    pass
                cbar_xz = fig.colorbar(pafib_xz,ax=ax2) 
                cbar_xz.set_label("Energy(keV)")
                #cbar_xz.set_clim(emin_track, emax_track)'''
                

                #plt.pause(0.1)
                plt.draw()

                ax3.cla()
                ax3.set_xlabel('Y (mm)')
                #ax3.set_ylabel('Z (mm)')
                ax3.set_xlim([-(TrackerLength*0.5)-1.,+(TrackerLength*0.5)+1.])
                ax3.set_ylim([-(TrackerWidth*0.5) -1.,+(TrackerWidth*0.5)+1. ])
                ax3.set_title("Tracker YZ View")
                ax3.tick_params(labelleft=False)
                #ax3.tick_params(left=False)

                ax3.plot(yp[:nnt], zp[:nnt], ":", color="red")
                pafib_yz = PatchCollection(patches_fib_yz, cmap=cmap, alpha=0.9, lw=0., norm=colors.Normalize(vmin=emin_track, vmax=emax_track))
                pafib_yz.set_array(efib1)
                ax3.add_collection(pafib_yz)

                # try:
                #     cbar_yz.remove()
                # except NameError:
                #     pass
                # cbar_yz = fig.colorbar(pafib_yz,ax=ax3) 
                # cbar_yz.set_label("Energy(keV)")
                # cbar_yz.set_clim(emin_track,emax_track)

                ########################DISPLAY LAYERS##########################
                diam = 2*FibRadius
                h = pow(3.,0.5)*FibRadius
                half_layer= diam + h
                gap = GapLayers
                # print("@@@@@@@", gap, diam)
                gap_um = "mm"
                '''if (gap_um == "cm"):
                    gap = gap*10 #in mm'''

                for l in range(int(Layers)):
                    '''a = l*(int(Views)*diam+gap)
                    layerZ_bottom = -(diam+gap+2*diam) + a
                    layerZ_top = -(diam+gap) + a'''
                    o_top = -(half_layer+gap)
                    layerZ_top = l*(2*half_layer + gap)+ o_top
                    layerZ_bottom = layerZ_top-2*half_layer

                    ax2.plot([-(TrackerLength*0.5),(TrackerLength*0.5)],[layerZ_top,layerZ_top], color = "darkgrey")
                    ax2.plot([-(TrackerLength*0.5),(TrackerLength*0.5)],[layerZ_bottom,layerZ_bottom], color = "darkgrey")
                    ax3.plot([-(TrackerLength*0.5),(TrackerLength*0.5)],[layerZ_top,layerZ_top], color = "darkgrey")
                    ax3.plot([-(TrackerLength*0.5),(TrackerLength*0.5)],[layerZ_bottom,layerZ_bottom], color = "darkgrey")
                
                ##################################################

                #plt.pause(0.1)
                plt.draw()

                #--- Calorimeter ---#
                hitpix = int(cal.DigiCalPixels_Hit)
                cpix_x = np.empty(hitpix, dtype=np.float)
                cpix_y = np.empty(hitpix, dtype=np.float)
                cpix_z = np.empty(hitpix, dtype=np.float)
                epix   = np.zeros(hitpix, dtype=np.float)
                patches_pix =  []
                
                jjcal = 0
                for jcal in range(hitpix):
                    calUnitNo = cal.CalUnitNo[jcal]
                    calPixX   = cal.CalPixelX[jcal]
                    calPixY   = cal.CalPixelY[jcal]
                    enePix    = cal.CalEnergy[jcal]

                    if(enePix<=0):
                        continue

                    epix[jjcal] = enePix
                    cpix_x[jjcal] = gxcpix[calUnitNo]
                    cpix_y[jjcal] = gycpix[calUnitNo]
                    

                    patches_pix.append(patches.Rectangle((cpix_x[jjcal],cpix_y[jjcal]), PixelSize,PixelSize))
                    jjcal += 1
                    
                ax4.cla() 
                ax4.set_xlabel('X (mm)')
                ax4.set_ylabel('Y (mm)')
                ax4.set_xlim([-(CaloLength*0.5)+1.,+(CaloLength*0.5)+1.])
                ax4.set_ylim([-(CaloLength*0.5)+1.,+(CaloLength*0.5)+1.])
                ax4.set_title("Calorimeter XY View")
               
                paPix = PatchCollection(patches_pix, cmap=cmap, alpha=0.9, lw=0., norm=colors.Normalize(vmin=emin_cal, vmax=emax_cal))
                paPix.set_array(epix[:jjcal])
                ax4.add_collection(paPix)
                ax4.plot(x0,y0, marker='X', markersize=6, markerfacecolor='w',markeredgewidth=1., markeredgecolor='k',alpha=0.70) 
                # try:
                #     cbar_cal.remove()
                # except NameError:
                #     pass
                cbar_cal = fig.colorbar(paPix,ax=ax4) 
                cbar_cal.set_label("Energy(MeV)")
                #cbar_cal.set_clim(emin_cal,emax_cal)

                plt.draw()
                plt.pause(0.1)
                plt.show()
                time.sleep(2)
                plt.savefig("./DisplayPlots/PDG"+str(pdg)+"_"+str(ngood)+".png")
                
        except KeyboardInterrupt as e:
            # print "Exit"
            raise


    # print "Number of Tracker good events = ", ngood
    raw_input("Press to Exit")

#sys.exit()

