# -*- coding: utf-8 -*-
"""
Created on Sun Nov 17 22:21:09 2019

Written by Zhenqi (Pete) Shi @ 
"""
import pandas as pd
import numpy as np
import pyphi as phi
import pyphi_plots as pp
import scipy.io as spio
import random

#NIRData=spio.loadmat('AbemaData.MAT')
#X = np.array(NIRData['SPECTRA_B'])
#Y = np.array(NIRData['Y'])
#
#X = phi.snv(X)
##X,M = phi.savgol(5,1,2,X)
#
#Y = Y
#
#BLOCK=5
#ITE=20
#LV=10

def CV (X,Y,BLOCK,ITER,LV):
    RMSECV=np.zeros((ITER,LV))
    
    for i in range(ITER):
        A=range(X.shape[0])
        B=list(A)
        random.shuffle(B)
        X=X[B,:]
        Y=Y[B,:]
        C=int(X.shape[0]/BLOCK)
        Yval_pred=np.zeros((X.shape[0],LV)) 
        
        for r in range(BLOCK):
            if r==0:
                BGN_v=0
                END_v=C
                D=range(X.shape[0])
                E=list(D)
                c_=E[C+1:X.shape[0]:1]
    #            print('A')
    
            elif r==BLOCK-1:
                BGN_v=C*r
                END_v=X.shape[0]
                D=range(X.shape[0])
                E=list(D)
                c_=E[0:C*r:1]
    #            print('hehe')
            else:
                BGN_v=C*r
                END_v=C*(r+1)
                D=range(X.shape[0])
                E=list(D)
                temp_1=E[0:C*r:1]
                temp_2=E[C*(r+1)+1:X.shape[0]:1]
                c_=[*temp_1,*temp_2]
    #            print('haha')
    #            
            print(BGN_v)
            print(END_v)
    #        print(c_)
            Xcal=X[c_,:]
            Ycal=Y[c_,:]        
            Xval=X[BGN_v:END_v,:]
            Yval=Y[BGN_v:END_v,:]               
            
            for k in range(LV):
                pls_calibration=phi.pls(Xcal,Ycal,k+1,mcsX='center',mcsY=True)
                yhatval_pls = phi.pls_pred(Xval,pls_calibration)
                Yval_pred[BGN_v:END_v,k]=yhatval_pls['Yhat'].ravel()          
        
        for k in range(LV):
            RMSECV[i,k]=np.sqrt(np.sum((Y.ravel()-Yval_pred[:,k])**2)/Y.shape[0])
    
    from matplotlib.pylab import plt
    for i in range(RMSECV.shape[0]):
        plt.plot(np.arange(1,LV+1,1),RMSECV[i,:])
        plt.title("Scree plot", fontsize=16, fontweight='bold')
        plt.xlabel("# of LV")
        plt.ylabel("RMSECV")
        plt.show()
    
    return RMSECV