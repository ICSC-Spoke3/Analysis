#useful function for CluTrack
import numpy as np

import os
import time
import glob
import math

def evaluateTrackPoint(x0,z0,cx,cz,zP,flag):
    t = (zP-z0)/cz
    x = x0 + cx*t
    x=round(x,2)
    if(flag==1):
        print( "t = ",t,"xTr = ",x)
    # print('evaluate_track_point ',x)
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
        StripIndex = round(NStripPerFiber*(iFib/2-1)+ stripFirstfib)
    #print("StripIndex",StripIndex)

    cstrip = c0 + pitch*(0.5+StripIndex)    
    #print("CENTRO STRIP:",cstrip)

    xminstrip = cstrip - pitch/2
    xmaxstrip = cstrip + pitch/2 
    #print ("xminstrip",xminstrip)
    #print ("xmaxstrip",xmaxstrip)
    
    strip_centres.append(round(cstrip,2))
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
        
        strip_centres.append(round(cstrip1,2))
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
        
        strip_centres.append(round(cstrip1,2))
        strip_index.append(StripIndex2)

        # print("cstrip1",round(cstrip1,2))
        #print("xmaxstrip2",xmaxstrip2)
        #print("strip index:",StripIndex2)

    return strip_centres, strip_index

def getStripCenter(Rfib, c0, offset, pitch, iFib, cfib, StripNo, stripIndex):
    xminfib = cfib - Rfib
    xmaxfib = cfib + Rfib

    NStripPerFiber = 2 * Rfib / pitch
    stripFirstfib = (2 * Rfib - offset) / pitch

    if iFib == 0:
        StripIndex = 0
    else:
        StripIndex = round(NStripPerFiber * (iFib / 2 - 1) + stripFirstfib)

    cstrip = c0 + pitch * (0.5 + StripIndex)

    xminstrip = cstrip - pitch / 2
    xmaxstrip = cstrip + pitch / 2

    while round(xminstrip, 2) > round(xminfib, 2) and xminstrip <= xmaxfib and xmaxstrip >= xminfib and StripIndex > 0:
        cstrip = cstrip - pitch
        xminstrip = cstrip - pitch / 2
        xmaxstrip = cstrip + pitch / 2
        StripIndex -= 1

    stripIndex += StripIndex

    cstrip = c0 + pitch * (0.5 + stripIndex)
    return cstrip

def FracArea(Rfib,pitch,cfib,cstrip):

    AFib = math.pi*pow(Rfib,2)
        
    #print("cfib",cfib)
    #print("cstrip",cstrip)
    #print("pitch",pitch)
    #print("Rfib",Rfib)
    
    cfib = cfib
    xminfib = cfib - Rfib
    xmaxfib = cfib + Rfib
    xminstrip = cstrip - pitch/2
    xmaxstrip = cstrip + pitch/2 
    
    #print("xminfib-xmaxfib",xminfib,xmaxfib)
    #print("xminstrip-xmaxstrip",xminstrip,xmaxstrip)
    
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
    #print("xmin - cfib",xmin - cfib)
    #print("xmax - cfib",xmax - cfib)
    #print("costh1",costh1)
    #print("costh2",costh2)


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

def PredValue (n,z,m,q):
    r = []
    for i in range(n):
        x_fit = m*z[i]+q
        r.append(x_fit)
    return r




