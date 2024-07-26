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
import pandas as pd

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
                #print("nPix: ",nPix)    
                #print("pixTot= ", totPix)
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
   
    PixelSize   = Calo_info[0]
    PixelDepth  = Calo_info[1]
    nPix        = Calo_info[2]
    TotPix      = Calo_info[3]
    #CalLayer =    Calo_info[4]
    CaloLength  = Calo_info[4]
    CaloWidth   = Calo_info[5]
   
    ### Load DATA from root file
    fr = TFile(froo)
    tpr=fr.Get("Primary")
    tkr=fr.Get("TrackerDigi")
    #cal=fr.Get("CalorimeterHit")
    cal=fr.Get("CalorimeterDigi")
    nevts = tkr.GetEntries()
    
   
    ### Read Data
    threshold = 0.9
    total_energy = 0
   
    ngood = 0
    energyZ = {}
   
    z_hit = []
    x_hit = []
    y_hit = []
    ene_hit = []
    R_hit = []
    
    #layer_z_values = []
  
    #total_energies = []
    #distances_90percent = []
    #Z =[]
    
    layer_z_values_all = []
    total_energies_all = []
    distances_90percent_all = []
    Z_all = []
   
   
   
    data_param = []
    h1_R1 = TH1F("h1_R1", "R1 Distribution", 50, 0, 100)
    h1_R2 = TH1F("h1_R2", "R2 Distribution", 50, 0, 100)
    #h1_R3 = TH1F("h1_R3", "R3 Distribution", 50, 0, 100)
    h1_R4 = TH1F("h1_R4", "R4 Distribution", 50, 0, 100)
    h1_R5 = TH1F("h1_R5", "R5 Distribution", 50, 0, 100)
    h1_R6 = TH1F("h1_R6", "R6 Distribution", 50, 0, 100)
    for i in range(10): 
        
        R1_values = []
        R2_values = []
        R3_values = []
        R4_values = []
        R5_values = []
        R6_values = []
        R7_values = []
        
        layer_z_values = []
        total_energies = []
        #distances_90percent = []
        Z = []
        
       
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
                    
                   
                    #writer.writerow([enePix])
                    cpix_x[jjcal] = gxcpix[calUnitNo]
                    cpix_y[jjcal] = gycpix[calUnitNo]
                    cpix_z[jjcal] = gzcpix[calUnitNo]
                    
                    z_hit[i].append(cpix_z[jjcal])  
                    x_hit[i].append(cpix_x[jjcal]) 
                    y_hit[i].append(cpix_y[jjcal]) 
                    ene_hit[i].append(epix[jjcal])
                    
                   
                    layer_z[jjcal] = int(cpix_z[jjcal]/PixelDepth)
                    layer_x[jjcal] = int(cpix_x[jjcal]/PixelDepth)
                    layer_y[jjcal] = int(cpix_y[jjcal]/PixelDepth)
                
                    layer_z_values.append(layer_z[jjcal])
                    
                    total_energy_event = np.sum(ene_hit[i])   
                    
                    jjcal+=1
                    
            energy90 = 0.9 * total_energy_event       
          
            for layer in range(1,26):
              total_energy_layer = 0
            
             
              if -layer in layer_z_values:
                for j, k in enumerate(layer_z_values):
                  
                  if k == -layer: 
                    total_energy_layer += epix[j]
                   
                
                total_energies.append(total_energy_layer)
               
                if total_energy_layer <= energy90 and  total_energy_layer>0 :
                       #print("energy layer = ", total_energy_layer)
                       #print(energy90)
                       R_molier = np.sqrt((cpix_x[j] - x0)**2 + (cpix_y[j] - y0)**2)
                       #print("R_moliere = ", R_molier)
                       #h1_R4.Fill(R_molier)
               
               
                Z.append(layer)
                max_energy_deposit_layer = max(total_energies)
                R3 = total_energy_layer/total_energy_event
                R2 = max_energy_deposit_layer/total_energy_event  
               
                R3_values.append(R3)
                #h1_R3.Fill(R3)
                #R_moliere_values.append(R_molier)
               
               
                #print("Layer", layer, "Energia depositata nel layer = ",  total_energy_layer) #, "Energia totale =", total_energy_event, "R3 = ", R3)
            
            #fuori loop layer   
            max_energy_index = np.argmax(total_energies)
            z_max_energy = layer_z_values[max_energy_index]
            #R_moliere_mean_value=sum( R_moliere_values)/len(R_moliere_values)
            #print("mean value moliere = ", R_moliere_mean_value)
            #print("Max energy = ",  max_energy_deposit_layer, "nel layer", z_max_energy)
            R6 = z_max_energy *PixelDepth #sono mm
            #print(R6)
            z_first_layer = Z[0] * PixelDepth
            z_last_layer = Z[-1] * PixelDepth    #sono mm
           
            #print( z_max_energy)
            layer_z_values_all.append(layer_z_values)
            total_energies_all.append(total_energies)
            #distances_90percent_all.append(distances_90percent)
            Z_all.append(Z) 
            energy_last_layer = total_energies_all[-1][-1]
            
            
            R4_values.append(R_molier)
            R1 = energy_last_layer/total_energy_event
            R1_values.append(R1)
            R2_values.append(R2)
            R5_values.append(z_last_layer)
            
            h1_R4.Fill(R_molier)
            h1_R1.Fill(R1)
            h1_R2.Fill(R2)
            h1_R5.Fill(z_last_layer)
            
            #R6 =  z_last_layer - z_first_layer #sviluppo longitudinale dello sciame (mm)
            #print("R6 = ", R6)
            #print("Sviluppo longitudinale dello sciame:", longitudinal_development)
            #R7 = (calPixX + calPixY + calPixZ)
            R6_values.append(R6)
            h1_R6.Fill(R6)
            #R7_values.append(R7)
            
            #data_param.append([R1_values, R2_values, R3_values,  R4_values,  R5_values,  R6_values,  R7_values])
            
            #data=pd.DataFrame(data_param,columns=['R1','R2', 'R3', 'R4', 'R5','R6', 'R7'])
            
            #data.to_csv('parametri.csv', index=False)

          
        except KeyboardInterrupt as e:
            #print ("Exit")
            raise
       
output_root_file = TFile("histograms.root", "RECREATE")
h1_R1.Write()
h1_R2.Write()
#h1_R3.Write()
h1_R4.Write()
h1_R5.Write()
h1_R6.Write()
output_root_file.Close()
    
histograms = [h1_R1, h1_R2, h1_R4, h1_R5, h1_R6]
titles = ["R1 Distribution e inclined", "R2 Distribution e inclined", "R4 Distribution e inclined" , "R5 Distribution e inclined", "R6 Distribution e inclined"]
for i, hist in enumerate(histograms):
        plt.figure()
        plt.hist(hist, bins=50, range=(0, 100), color= 'blue')
        plt.title(titles[i])
        #plt.xlabel('Value')
        #plt.ylabel('Frequency')
        #plt.grid(True)
        plt.savefig(f'histogram_{i + 1}.png')
        plt.show()