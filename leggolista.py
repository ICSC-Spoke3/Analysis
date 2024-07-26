
import pylab
import math
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import gridspec
import matplotlib.colors as colors
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
from matplotlib.axes import Axes
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import os
import sys
import time
from ROOT import TFile, TTree, TChain, TH1D, TMath, TCanvas, TH1F, TH2F
import array as ary
import numpy as np
import glob

if __name__ == '__main__':
    
    list = sys.argv[1]

    eventNo = []
    particle = []
    gEnergy = []
    gEnergy_um = []
    LayNo = []
    gap_mm = []

    with open(list) as fptr:
        fnamedata = fptr.readlines()
        #print ("fname: ", fnamedata)

        for jf in range(len(fnamedata)):
            rootfile = fnamedata[jf].split('\n')[0]
            print ("[INFO] ROOT file\t",rootfile)
            #print (rootfile.split('_')[5])
            eventNo.append(int((rootfile.split('_')[5]).split('-')[0]))
            particle.append((rootfile.split('_')[1]).split('/')[4])
            gEnergy.append((rootfile.split('_')[2]).split('-')[0])
            gEnergy_um.append(((rootfile.split('_')[2]).split('-')[1]))
            LayNo.append(int((rootfile.split('_')[3]).split('-')[0]))
            gap_mm.append((rootfile.split('_')[4]).split('-')[0])
            #print (gEnergy)
            #print (gEnergy_um)
            #print (LayNo)
            #print (gap_mm)
            #print (particle)
            #print (eventNo)


    

    
    
