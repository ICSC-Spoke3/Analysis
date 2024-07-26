#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 17:26:26 2023

@author: federica
"""
import pandas as pd
import os
import numpy as np
import seaborn as sns
import sys
import matplotlib.pyplot as plt

#read the csv files produced by CluTrack.py macro 
#re-arrange data in suitable shape for ML algorithm

if __name__ == "__main__":
    
    try:
        mainDir = sys.argv[1]
    except:
        print("Please insert the folder in which the csv data are stored") 
        sys.exit(1)
    #InputFile
    fin = sys.argv[2]  
    if('csv' in fin ):
        print('Please insert the file name without <<.csv>>')
        sys.exit()
    
    
    data_file=mainDir+fin+'.csv'
    data_trk=pd.read_csv(data_file,index_col=0)
    data_trk=data_trk.drop('zMC',axis=1)
    if 'viewx' in fin:
        column_order=['xrec1','zxrec1','xrec2','zxrec2','xrec3','zxrec3','dx1','dxz1','dx2','dxz2','dx3','dxz3','xMC1','xMC2','xMC3']
        # Genera un numero sequenziale per ogni ripetizione di un evento
        data_trk['count'] = data_trk.groupby('Ev').cumcount() + 1
    
        # Crea un dataframe pivotato utilizzando il metodo pivot_table
        data_trk_processed = data_trk.pivot_table(index='Ev', columns='count', values=['xrec', 'zxrec','dx','dxz','xMC'])
    
        # Rinomina le colonne del dataframe pivotato
        data_trk_processed.columns = [f'{col[0]}{col[1]}' for col in data_trk_processed.columns]
        data_trk_processed=data_trk_processed.reindex(columns=column_order)
        data_trk_processed.reset_index(inplace=True)
        data_trk_processed.insert(1,'view',0)
        data_trk_processed=data_trk_processed.drop('view',axis=1)
        # data_trk_processed=data_trk_processed.drop('xMC2',axis=1)
        # data_trk_processed=data_trk_processed.drop('xMC3',axis=1)
    else:       
        column_order=['yrec1','zyrec1','yrec2','zyrec2','yrec3','zyrec3','dy1','dyz1','dy2','dyz2','dy3','dyz3','yMC1','yMC2','yMC3']
        # Genera un numero sequenziale per ogni ripetizione di un evento
        data_trk['count'] = data_trk.groupby('Ev').cumcount() + 1
    
        # Crea un dataframe pivotato utilizzando il metodo pivot_table
        data_trk_processed = data_trk.pivot_table(index='Ev', columns='count', values=['yrec', 'zyrec','dy','dyz','yMC'])
    
        # Rinomina le colonne del dataframe pivotato
        data_trk_processed.columns = [f'{col[0]}{col[1]}' for col in data_trk_processed.columns]
        data_trk_processed=data_trk_processed.reindex(columns=column_order)
        data_trk_processed.reset_index(inplace=True)
        data_trk_processed.insert(1,'view',0)
        data_trk_processed=data_trk_processed.drop('view',axis=1)
        # data_trk_processed=data_trk_processed.drop('xMC2',axis=1)
        # data_trk_processed=data_trk_processed.drop('xMC3',axis=1)
       
        
        
    data_trk_processed=data_trk_processed.drop('Ev',axis=1)
    data_trk_processed.to_csv(mainDir+fin+'processed.csv')
    print(fin)
    print(data_trk_processed)
       

    