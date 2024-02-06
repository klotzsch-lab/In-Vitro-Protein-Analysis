# -*- coding: utf-8 -*-
"""
Python script written by Willi Weber
Humboldt-Universität zu Berlin
Mechanobiology (Prof. Klotzsch)
10.05.2023

last updated: 16.01.2024

input: .csv file with intensities of FRAP trajectories

recommended to use with SpyderIDE
"""

from PyQt5.QtWidgets import*
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

from scipy.optimize import curve_fit

#function for fitting
def neg_exp(x,a,b):
    y = a*(1-np.exp(-x*b))
   
    return y

#Data import
dlg = QFileDialog()
dlg.setFileMode(QFileDialog.ExistingFiles)
paths= []
files = {}
if dlg.exec_():
    paths = dlg.selectedFiles()
    
    for path in paths:
        path_alias = path.split("/")[-1].replace(".csv","")
        files[path_alias] = pd.read_csv(path, sep = ";")

#plot setup
fig, ax = plt.subplots()
colors = {"PLD CP29A": "#3ac966", "DPLD CP29A": "#3ebae6", "FL CP29A": "#edc142"}
colors2 = {"PLD CP29A": "#8BCD9F50", "DPLD CP29A": "#A6D7E850", "FL CP29A": "#F7DD9350"}
ax.set_xlabel("Time (s)")
ax.set_ylabel("Normalized Fluorescence Intensity (%)")
ax.set_xlim(0,50)
ax.set_ylim(0,120)

#data processing
value_dict = []
linecollection = []
linelabels = []
for file in files:
    tempdata = files[file]
    tempdata.reset_index()
    
    #data normalization
    for column in tempdata:
        normvalue = tempdata[column][0]
        tempdata[column] = tempdata[column].apply(lambda x: 100*x/normvalue)
    
    #calculation of mean and standard deviation
    stdData = tempdata.std(axis= 1).rename("std")
    meanData = tempdata.mean(axis = 1).rename("mean")
    meanPostBleach = meanData[15]
    
    #data transformation for fit
    fitMeanData = meanData.apply(lambda x: x-meanPostBleach)
    meanData = pd.concat([meanData, stdData ], axis = 1)
    fitdata = pd.concat([fitMeanData, stdData ], axis = 1)
    fitdata = fitdata[15:].reset_index()

    #fitting
    parameters, covariance = curve_fit(neg_exp, fitdata.index.values.tolist(), fitdata["mean"], bounds=(0, [100., np.inf])) 
    perr = np.sqrt(np.diag(covariance))
    a = parameters[0] 
    b = parameters[1] 
    
    #saving parameters for export
    corA = a + meanPostBleach
    t_half = np.log(0.5)/-b
    value_dict.append({"File": file, "Mobile Fraction": corA, "Mobile Fraction Error": perr[0], "Tau": b, "Tau Error": perr[1], "t1/2": t_half, "t1/2err_upper":t_half - np.log(0.5)/-(b+perr[1]), "t1/2err_lower":np.abs(t_half - np.log(0.5)/-(b-perr[1])), "n_trajectorys": len(tempdata.columns)})
    
    #preparing plot data from fit
    xs = range(0,32)
    ys = [neg_exp(x,a,b) for x in xs]
    ys = [y+meanPostBleach for y in ys]
    plotxs = [x+15 for x in xs]
    dataColor = colors2[file]
    fitColor = colors[file]
        
    #halflife drawing:
    ax.axhline(y=meanPostBleach, xmin = 15/50, xmax = (15+t_half)/50,linestyle = "--", color= "black",alpha =0.2)
    ax.axvline(x = 15+t_half, ymin = meanPostBleach/120, ymax = (neg_exp(t_half,a,b)+meanPostBleach)/120, linestyle = "--", color= "black", alpha =0.2)

    #data drawing
    meanData[:5].plot( y = "mean", ax = ax, color = dataColor)
    meanData[15:].plot( y = "mean", ax = ax,color = dataColor)
    
    #std drawing
    md = meanData
    errorFillTopStart = [md["mean"][i]+md["std"][i] for i in range(5)]
    errorFillBottomStart = [md["mean"][i]-md["std"][i] for i in range(5)]
    ax.fill_between(xs[:5],errorFillTopStart, errorFillBottomStart, color = dataColor)

    errorFillTop = [md["mean"][15:][i+15]+md["std"][15:][i+15] for i in range(len(md[15:]))]
    errorFillBottom = [md["mean"][15:][i+15]-md["std"][15:][i+15] for i in range(len(md[15:]))]
    ax.fill_between(plotxs,errorFillTop, errorFillBottom, color = dataColor)
        
    #fit drawing
    line = ax.plot(plotxs,ys,c = fitColor)
    linecollection.append(line[0])
    linelabels.append(str(file))
     

ax.legend(linecollection, linelabels, loc="lower right")
ax.grid(which = "minor")
valueFrame = pd.DataFrame(value_dict)