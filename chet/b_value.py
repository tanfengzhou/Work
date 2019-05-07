#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 16:30:21 2017

@author: Jesse
"""
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import stats



def exponential_func(x, a, b, c):
    return a * np.exp(-b * x) + c

def bvalue_func(event_MLmag,mc,binwidth):

# Input:
# - event_MLmag: an array of  event ML magnitudes.
# - mc: magnitude of completeness
# - binwidth: desired width of histogram bins

#  The number of bins is calculated based on the minimum and maximum magnitude
#  of the data set
    nbin = np.arange(math.floor(min(event_MLmag)),math.ceil(max(event_MLmag)),binwidth)
    ML_hist = plt.hist(event_MLmag, nbin);
    ML   = np.array(ML_hist[1])
    ML   = ML[:len(nbin)-1] + np.diff(ML) / 2.0
#    ML   = ML[:nbin] + np.diff(ML) / 2.0
    N    = np.array(ML_hist[0]);
    
    Ncum = np.zeros(ML.shape)
    
    for iml in range(len(ML)):
      # cumulative number of events for M greater than a given value:
      Ncum[iml] = np.sum(np.sum(event_MLmag >=ML[iml]))
    
    # estimate a- and b-values for the magnitude of completeness given above:
    imc   = (ML >= mc) & (ML <= max(event_MLmag))
    ML_   = ML[imc]
    Ncum_ = Ncum[imc]
    slope,intercept,r_value,p_value,std_err = stats.linregress(ML_,np.log10(Ncum_))
    a = intercept
    b = abs(slope)
    iNcum_ = 10**(a - b * ML_)   # straight line using a- and b-values.
    R2   = r_value**2  # goodness of fit
    fig = plt.figure(figsize=(8.27,8.27))
    ax1 = fig.add_subplot(1,1,1)
    ax1.bar(ML,N,edgecolor="k",color=[0.25,0.25,0.25],width=binwidth)
    ax1.set_xlabel("Magnitude, ML")
    ax1.set_ylabel("Number of Events")
    ax1.set_ylim(0,max(N)*1.2)
    ax1.grid(True,linestyle='--')
    
    ax2 = ax1.twinx()
    line1, = ax2.plot(ML,np.log10(Ncum),'bo',label='Input Data')
    line2, = ax2.plot(ML_,np.log10(iNcum_),'r--',label='Fitted Line')
    ax2.set_ylabel("log10(N > ML)",color='blue')
    ax2.spines['right'].set_color('blue')
    ax2.tick_params(axis='y',colors='blue')
    plt.legend(handles=[line1,line2],loc=1)
    textbox = 'log10(N) = '+'%.3f' % (a)+' - '+'%0.3f' % (b)+'ML \nR2 = '+'%0.2f' % (R2)+'\nMc = '+'%0.2f' % (mc)
    #ax2.text(2,2.5,textbox,bbox={'facecolor':'white', 'alpha':0.5, 'pad':10},fontsize=12)
    ax2.text(0.03,0.97,textbox,bbox={'facecolor':'white', 'alpha':0.5, 'pad':10},fontsize=12,verticalalignment='top',horizontalalignment='left',transform=ax2.transAxes)
    plt.title('Entire Dataset')
    plt.show()
    fig.savefig('test.png')
    return a,b,R2

def find_best_fit(event_MLmag,low_mc,high_mc,interval):
#  This code will cycle through a range of magnitude of completeness' and find
#  the best fit for a set of data.
#  Input:
#  - event_MLmag: an array of  event ML magnitudes.
#  - low_mc: the lowest magnitude of completeness to begin calculations
#  - high_mc: the highest magnitude of completeness to finish calculations
#  - interval: the interval between mc values, and also the binwidth for histograms
    b_all = []
    R2_all = []
    mc_all = []
#  Calculate the R2, b-value, a-value, and magnitude of completeness for a range of data
    for i in range(int(low_mc/interval),int(high_mc/interval)):
        a,b,R2 = bvalue_func(event_MLmag,i*interval,interval)
        b_all.append(b)
        R2_all.append(R2)
        mc_all.append(i*interval)
#  Show R2 values on a curve plot
    plt.plot(mc_all,R2_all)
    plt.show()
#  Calculate best fit
    nPoints = len(R2_all)
    allCoord = np.vstack((range(nPoints), R2_all)).T
    np.array([range(nPoints), R2_all])
    firstPoint = allCoord[0]
    lineVec = allCoord[-1] - allCoord[0]
    lineVecNorm = lineVec / np.sqrt(np.sum(lineVec**2))
    vecFromFirst = allCoord - firstPoint
    scalarProduct = np.sum(vecFromFirst * np.tile(lineVecNorm, (nPoints, 1)), axis=1)
    vecFromFirstParallel = np.outer(scalarProduct, lineVecNorm)
    vecToLine = vecFromFirst - vecFromFirstParallel
    distToLine = np.sqrt(np.sum(vecToLine ** 2, axis=1))
    idxOfBestPoint = np.argmax(distToLine)
#  Plot the histogram for the best-fit data
    bvalue_func(event_MLmag,mc_all[idxOfBestPoint],interval)
    return R2_all[idxOfBestPoint],mc_all[idxOfBestPoint],b_all[idxOfBestPoint]
   
if __name__=='__main__':
    event_MLmag = np.random.normal(1.25,1,1000)
    interval = 0.1
    low_mc = 0
    high_mc = 3
    R2,mc,b = find_best_fit(event_MLmag,low_mc,high_mc,interval)
   
  