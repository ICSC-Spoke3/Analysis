#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 15:58:44 2023

@author: federica
"""
#test of the possible regression algorithm with ML
import pandas as pd
import os
import numpy as np
import seaborn as sns
import sys
import matplotlib.pyplot as plt

from sklearn.model_selection import StratifiedKFold,RepeatedStratifiedKFold,RepeatedKFold
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression,LinearRegression,Ridge,Lasso,ElasticNet, RANSACRegressor
from sklearn.metrics import recall_score,roc_auc_score,f1_score,RocCurveDisplay,roc_curve,auc,accuracy_score
from sklearn.model_selection import cross_val_score,KFold,learning_curve
from sklearn.metrics import mean_absolute_error, r2_score, mean_squared_error

from sklearn.ensemble import RandomForestRegressor
from sklearn.multioutput import MultiOutputRegressor
import xgboost as xg

import ROOT

#%% read the preprocessed dataset

if __name__ == "__main__":

    mainDir = sys.argv[1]
    
    #InputFile
    try:
        fin = sys.argv[2]         
    except:
        print("ERROR MESSAGE: \n =====> Insert first argument = the folder which contains the data; Insert second argument = the name of your csv data file")
        sys.exit(1)
    best_model=LinearRegression()
    ndata=1000
    if sys.argv[3] !='':
        best_model =sys.argv[3]
    if sys.argv[4]!='':
        ndata=sys.argv[4]
    
    trk_data_x=pd.read_csv(mainDir+fin,index_col=0)       
    trk_data_x=trk_data_x.dropna()
    trk_data_x=trk_data_x.iloc[:int(ndata),:]
    #select a 10% of dataset to be a validation independent dataset
    trk_data_x_cv=trk_data_x.sample(frac=0.9,random_state=0)
    #trk_data_x_cv=trk_data_x_cv.reset_index()
    trk_data_x_val=trk_data_x.drop(trk_data_x_cv.index)
    #trk_data_x_val=trk_data_x_val.reset_index()
    col_name=trk_data_x.columns
    tot_ncol=len(col_name)
    targ_col=[]
   
    for i in range(0,len(col_name)):
        if 'MC' in col_name[i]:
            targ_col.append(i)
            
            
#from here we set the valriables for ML
    targ=trk_data_x_cv.iloc[:,targ_col[0]:] #target value=MC truth
    X=trk_data_x_cv.iloc[:,:len(col_name)-len(targ_col)]   #data features
    
    
    
    n_plits=5
    n_repeats=100
    
    kf=RepeatedKFold(n_splits=n_plits,n_repeats=n_repeats,random_state=0)

    size=n_plits*n_repeats
    r2=np.zeros(size)
    mae=np.zeros(size)
    mse=np.zeros(size)
    r2_train=np.zeros(size)
    mae_train=np.zeros(size)
    mse_train=np.zeros(size)
    #list to memorize some predicted and true data
    y_true=list()
    y_pred=list()
    ytrue_train=list()
    y_train_pred=list()
    print('cross val ')

    for i, (train_index, test_index) in enumerate(kf.split(X,targ)):
        
        X_train = X.iloc[train_index]
        X_test = X.iloc[test_index]
        y_train = targ.iloc[train_index].values
        y_test = targ.iloc[test_index].values
              
        if(best_model=='LinearRegressor'):
            print('Linear regression selected')
            multioutput_regressor=MultiOutputRegressor(LinearRegression(fit_intercept=(True)))
            best_model=multioutput_regressor
            best_model.fit(X_train,y_train)
                
        if(best_model=='Ridge'):
            print('ridge model selected')
            multioutput_regressor=MultiOutputRegressor(Ridge())
            space = dict()
            space['estimator__alpha']=[0.20, 0.5]
            cv_inner = KFold(n_splits=3, shuffle=True, random_state=1)
            search = GridSearchCV(multioutput_regressor, space, scoring='r2', cv=cv_inner, refit=True)
            # # execute search
            result = search.fit(X_train, y_train)
            ## get the best performing model fit on the whole training set
            best_model = result.best_estimator_
            best_model.fit(X_train,y_train)

        if(best_model=='Lasso'):
            print('Lasso model selected')
            multioutput_regressor=MultiOutputRegressor(Lasso())
            space = dict()
            space['estimator__alpha']=[0.2, 0.5,1]
            cv_inner = KFold(n_splits=3, shuffle=True, random_state=1)
            search = GridSearchCV(multioutput_regressor, space, scoring='r2', cv=cv_inner, refit=True)
            # # execute search
            result = search.fit(X_train, y_train)
            ## get the best performing model fit on the whole training set
            best_model = result.best_estimator_
            best_model.fit(X_train,y_train)
            
        if(best_model=='ElasticNet'):
            print('ElasticNet model selected')
            multioutput_regressor=MultiOutputRegressor(ElasticNet())
            space = dict()
            space['estimator__alpha']=[0.2, 0.5,1]
            space['estimator__l1_ratio']=[0.1,0.5]
            cv_inner = KFold(n_splits=3, shuffle=True, random_state=1)
            search = GridSearchCV(multioutput_regressor, space, scoring='r2', cv=cv_inner, refit=True)
            # # execute search
            result = search.fit(X_train, y_train)
            ## get the best performing model fit on the whole training set
            best_model = result.best_estimator_
            best_model.fit(X_train,y_train)

        if(best_model=='RandomForest'):
            print('Random Forest model selected')
            multioutput_regressor=MultiOutputRegressor(RandomForestRegressor())
            space = dict()
            space['estimator__n_estimators'] = [50,80,100]
            space['estimator__max_depth'] = [None,3]
            cv_inner = KFold(n_splits=3, shuffle=True, random_state=1)
            search = GridSearchCV(multioutput_regressor, space, scoring='r2', cv=cv_inner, refit=True)
            # # execute search
            result = search.fit(X_train, y_train)
            ## get the best performing model fit on the whole training set
            best_model = result.best_estimator_
            best_model.fit(X_train,y_train)
            
        if(best_model=='XGBoost'):
            print('XGBoost model selected')
            multioutput_regressor=MultiOutputRegressor(xg.XGBRegressor())
            space = dict()
            space['estimator__n_estimators'] = [50,80,100]
            space['estimator__learning_rate'] = [0.1,0.2,0.3]
            space['estimator__reg_alpha']=[0, 0.5, 1, 5]
            space['estimator__reg_lambda']=[0, 0.5, 1, 5]
            cv_inner = KFold(n_splits=3, shuffle=True, random_state=1)
            search = GridSearchCV(multioutput_regressor, space, scoring='r2', cv=cv_inner, refit=True)
            # # execute search
            result = search.fit(X_train, y_train)
            ## get the best performing model fit on the whole training set
            best_model = result.best_estimator_
            best_model.fit(X_train,y_train)
        
        p=best_model.predict(X_test)
        p_train=best_model.predict(X_train)
        y_true.append(y_test)
        y_pred.append(p)
        ytrue_train.append(y_train)
        y_train_pred.append(p_train)
        
        r2[i]=r2_score(y_test, p)
        mae[i]=mean_absolute_error(y_test, p)
        mse[i]=mean_squared_error(y_test,p)

        r2_train[i]=r2_score(y_train, p_train)
        mae_train[i]=mean_absolute_error(y_train, p_train)
        mse_train[i]=mean_absolute_error(y_train,p_train)
        
    print('R2  : %.3f +/- (%.3f)'%(np.mean(r2),np.std(r2)))
    print('MAE : %.3f +/- (%.3f)'%(np.mean(mae),np.std(mae)))
    print('MSE : %.3f +/- (%.3f)'%(np.mean(mse),np.std(mse)))
    
    print('MAE on train : %.3f  +/- (%.3f): ' % (np.mean(mae_train),np.std(mae_train)))
    print('R2 on train : %.3f  +/- (%.3f): ' % (np.mean(r2_train),np.std(r2_train)))
    

#validation test
    print('validation start')
    X_train_val=trk_data_x_cv.iloc[:,:len(col_name)-len(targ_col)]
    y_train_val=trk_data_x_cv.iloc[:,targ_col[0]:] 
    X_test_val=trk_data_x_val.iloc[:,:len(col_name)-len(targ_col)]
    y_test_val=trk_data_x_val.iloc[:,targ_col[0]:]
    
   
    if(best_model=='LinearRegressor'):
        print('Linear regression model')
        multioutput_regressor=MultiOutputRegressor(LinearRegression(fit_intercept=(True)))
        best_model=multioutput_regressor
        best_model.fit(X_train_val,y_train_val)
        
    if(best_model=='Ridge'):
        multioutput_regressor=MultiOutputRegressor(Ridge())
        print('ridge model selected')
        space = dict()
        space['estimator__alpha']=[0.20, 0.5]
        cv_inner = KFold(n_splits=3, shuffle=True, random_state=1)
        search = GridSearchCV(multioutput_regressor, space, scoring='r2', cv=cv_inner, refit=True)
        # # execute search
        result = search.fit(X_train_val, y_train_val)
        ## get the best performing model fit on the whole training set
        best_model = result.best_estimator_
        best_model.fit(X_train_val,y_train_val)

    if(best_model=='Lasso'):
        print('Lasso model selected')
        multioutput_regressor=MultiOutputRegressor(Lasso())
        space = dict()
        space['estimator__alpha']=[0.2, 0.5,1]
        cv_inner = KFold(n_splits=3, shuffle=True, random_state=1)
        search = GridSearchCV(multioutput_regressor, space, scoring='r2', cv=cv_inner, refit=True)
        # # execute search
        result = search.fit(X_train_val, y_train_val)
        ## get the best performing model fit on the whole training set
        best_model = result.best_estimator_
        best_model.fit(X_train_val,y_train_val)
        
    if(best_model=='ElasticNet'):
        print('ElasticNet model selected')
        multioutput_regressor=MultiOutputRegressor(ElasticNet())
        space = dict()
        space['estimator__alpha']=[0.2, 0.5,1]
        space['estimator__l1_ratio']=[0.1,0.5]
        cv_inner = KFold(n_splits=3, shuffle=True, random_state=1)
        search = GridSearchCV(multioutput_regressor, space, scoring='r2', cv=cv_inner, refit=True)
        # # execute search
        result = search.fit(X_train_val, y_train_val)
        ## get the best performing model fit on the whole training set
        best_model = result.best_estimator_
        best_model.fit(X_train_val,y_train_val)
        
    if(best_model=='RandomForest'):
        print('Random Forest model selected')
        multioutput_regressor=MultiOutputRegressor(RandomForestRegressor())
        space = dict()
        space['estimator__n_estimators'] = [50,80,100]
        space['estimator__max_depth'] = [None,3]
        cv_inner = KFold(n_splits=3, shuffle=True, random_state=1)
        search = GridSearchCV(multioutput_regressor, space, scoring='r2', cv=cv_inner, refit=True)
        # # execute search
        result = search.fit(X_train_val, y_train_val)
        ## get the best performing model fit on the whole training set
        best_model = result.best_estimator_
        best_model.fit(X_train_val,y_train_val)
        
    if(best_model=='XGBoost'):
        print('XGBoost model selected')
        multioutput_regressor=MultiOutputRegressor(xg.XGBRegressor())
        space = dict()
        space['estimator__n_estimators'] = [50,80,100]
        space['estimator__learning_rate'] = [0.1,0.2,0.3]
        space['estimator__reg_alpha']=[0, 0.5, 1, 5]
        space['estimator__reg_lambda']=[0, 0.5, 1, 5]
        cv_inner = KFold(n_splits=3, shuffle=True, random_state=1)
        search = GridSearchCV(multioutput_regressor, space, scoring='r2', cv=cv_inner, refit=True)
        # # execute search
        result = search.fit(X_train_val, y_train_val)
        ## get the best performing model fit on the whole training set
        best_model = result.best_estimator_
        best_model.fit(X_train_val,y_train_val)

    print(best_model)
    y_train_pred=best_model.predict(X_train_val)
    y_test_pred=best_model.predict(X_test_val)
    
    r2_val=r2_score(y_test_val, y_test_pred)
    mae_val=mean_absolute_error(y_test_val, y_test_pred)
    mse_val=mean_squared_error(y_test_val, y_test_pred)
    print('on validation r2 : %.3f '%r2_val)
    print('on validation mae: %.3f '%mae_val)
    print('on validation mse: %.3f '%mse_val)
    
    view=''
    if 'viewx' in fin:
        view='viewx'
    else:
        view='viewy'
        
        #learning curve
    print('learningcurve evaluation')
    X_tr_ln,X_test_ln,y_train_ln,y_test_ln=train_test_split(X,targ,test_size=0.20,random_state=0)
    train_sizes = np.linspace(0.1, 1.0, 10)
    
    # Calcola i punteggi di training e di validazione tramite le learning curve
    train_sizes, train_scores, test_scores = learning_curve(
        best_model, X_train_val, y_train_val, train_sizes=train_sizes, cv=5)
    
    # Calcola i punteggi medi e le deviazioni standard per i punteggi di training e di validazione
    train_scores_mean = np.mean(train_scores, axis=1)
    train_scores_std = np.std(train_scores, axis=1)
    test_scores_mean = np.mean(test_scores, axis=1)
    test_scores_std = np.std(test_scores, axis=1)
    
    # Plotta i risultati delle learning curve
    fig3,pl3=plt.subplots(figsize=(8,6))
    pl3.set_title("Learning Curve")
    pl3.set_xlabel("Training Examples")
    pl3.set_ylabel("Score")
    pl3.grid()
    
    pl3.fill_between(train_sizes, train_scores_mean - train_scores_std,
                     train_scores_mean + train_scores_std, alpha=0.1, color="r")
    pl3.fill_between(train_sizes, test_scores_mean - test_scores_std,
                     test_scores_mean + test_scores_std, alpha=0.1, color="g")
    pl3.plot(train_sizes, train_scores_mean, 'o-', color="r", label="Training Score")
    pl3.plot(train_sizes, test_scores_mean, 'o-', color="g", label="Cross-validation Score")
    
    pl3.legend(loc="best")
    
    fig3.savefig('../Results/plot/'+best_model.estimator.__class__.__name__+'_'+view+'_ln_curves.png')
    
    
    #___________visualization plots_______________#
    # fig2,pl2=plt.subplots()
    # resisuals=[]
    # for i in range(0,len(y_test_pred)):
    #     resisuals.append(float(y_test_pred[i][0])-float(y_test[i][0]))
    # pl2.hist(resisuals,100)
    # fig2.savefig('./'+'residuals.png')
        

    resisuals_x1=[]
    resisuals_x2=[]
    resisuals_x3=[]
    

    if view=='viewx':
        h_Res_viewx_x1=ROOT.TH1F('Residual first x1','Residual first x1',100,-10,10)
        h_Res_viewx_x2=ROOT.TH1F('Residual first x2','Residual first x2',100,-10,10)
        h_Res_viewx_x3=ROOT.TH1F('Residual first x3','Residual first x3',100,-10,10)

    if view=='viewy':
        h_Res_viewx_x1=ROOT.TH1F('Residual first y1','Residual first y1',100,-10,10)
        h_Res_viewx_x2=ROOT.TH1F('Residual first y2','Residual first y2',100,-10,10)
        h_Res_viewx_x3=ROOT.TH1F('Residual first y3','Residual first y3',100,-10,10)
        
      
    if view=='viewx':
        for i in range(0,len(y_test_pred)):
            resisuals_x1.append(float(y_test_pred[i][0])-float(np.array(y_test_val.values)[i][0]))
            resisuals_x2.append(float(y_test_pred[i][1])-float(np.array(y_test_val.values)[i][1]))
            resisuals_x3.append(float(y_test_pred[i][2])-float(np.array(y_test_val.values)[i][2]))
            h_Res_viewx_x1.Fill(float(y_test_pred[i][0])-float(np.array(y_test_val.values)[i][0]))
            h_Res_viewx_x2.Fill(float(y_test_pred[i][1])-float(np.array(y_test_val.values)[i][1]))
            h_Res_viewx_x3.Fill(float(y_test_pred[i][2])-float(np.array(y_test_val.values)[i][2]))
    else:
        for i in range(0,len(y_test_pred)):
            resisuals_x1.append(float(y_test_pred[i][0])-float(np.array(y_test_val.values)[i][0]))
            resisuals_x2.append(float(y_test_pred[i][1])-float(np.array(y_test_val.values)[i][1]))
            resisuals_x3.append(float(y_test_pred[i][2])-float(np.array(y_test_val.values)[i][2]))
            h_Res_viewx_x1.Fill(float(y_test_pred[i][0])-float(np.array(y_test_val.values)[i][0]))
            h_Res_viewx_x2.Fill(float(y_test_pred[i][1])-float(np.array(y_test_val.values)[i][1]))
            h_Res_viewx_x3.Fill(float(y_test_pred[i][2])-float(np.array(y_test_val.values)[i][2]))        
   

    myfile = ROOT.TFile('../Results/plot/tracking_fiber_ai'+view+best_model.estimator.__class__.__name__+'.root', 'RECREATE' )
   
    
    h_Res_viewx_x1.Write()
    h_Res_viewx_x2.Write()
    h_Res_viewx_x3.Write()

      
    myfile.Close()
    

