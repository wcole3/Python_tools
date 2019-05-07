#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 11:33:12 2018

Fitting routine for Langmuir isotherm for use with Malchite green SHS data
The isotherm is modified to account for adsorbate bulk concentration change 
resulting from adsorption

Fitting of displacement studies to the 2 species, modified Langmuir isotherm is
now supported

User can choose functionality and generated plots in main()

@author: wcole
"""

import numpy as np;
import matplotlib.pyplot as plt;
import scipy.optimize as sp;

#import the dataset
#format will be two columns of [Conc., Int.]
#calibration data
caliDataName = 'nPSB_MG_sp';
#displacement exp data
dispDataName = 'nPSB_PF_MG_SP';
#the fixed dye concentration for the displacement data (in the same units as dispData)
fixed_dye_conc = 1000;

def LangmuirCurve(x, B,a,N,K):
    """Calculate the Langmuir curve for the calibration data set
    Parameters
    ----------
    x : The data x range
    
    B : the B fit constant
    
    a : the a fit constant
    
    N: the maximum surface density of dye molecules
    
    K : the equillibrium constant for the adsorption equation
    
    Returns
    -------
    The value of the Langmuir curve at x
    
    
    """
    return B+np.square(a*(((x+N+(55.5/K))-np.sqrt(np.square(x+N+(55.5/K))-(4*x*N)))/(2*N)));

def Frac_Cov(c, N, K):
    """The curve defining the fractional coverage of the dye
    Parameters
    ---------
    c : the concentration of the dye
    
    N: the maximum surface density of dye molecules
    
    K : the equillibrium constant for the adsorption equation
    
    Returns
    -------
    The value of the fractional coverage the concentration
    
    """
    return (((c+N+(55.5/K))-np.sqrt(np.square(c+N+(55.5/K))-(4*c*N)))/(2*N));

def DispCurve(c,x,B,a,N,K):
    """The displacement Langmuir model
    Parameters
    ---------
    x : The value of x
    
    c : The concentration of the displacing molecule
    
    B : the B fit constant
    
    a : the a fit constant
    
    N: the maximum surface density of dye molecules
    
    K : the equillibrium constant for the adsorption equation
    
    Returns
    -------
    The expected SHS intensity resulting from competitive adsorption
    """
    return B+np.square(a*(x/(1+x+(K*((c-N)/55.5)))));


def DispDyeCoverage(c_1, c_2, N_d, K_1, K_2):
    """The displacing or dye molecule's coverage.  Enter desired molecule as
    species 1.
    Parameters
    ---------
    c_1 : The concentration of species 1
        
    c_2 : The concentration of species 2
        
    N_d : The maximum surface number denisty of the Dye
        
    K_1 : The equillibirum constant of species 1 adsorption
        
    K_2 : The equillibirum constant of species 1 adsorption
    
    Returns
    ------
    The surface coverage of species 1
    """
    x = (K_1 * ((c_1 - N_d)/55.5));
    y = (K_2 * ((c_2 - N_d)/55.5));
    return (x / (1 + x + y));

def IntToCoverage(intensity, caliConst):
    """
    Convert SHS intensity to fractional coverage
    Parameters
    ----------
    intensity : the SHS intensity
    
    caliConst : the constants obtained from fitting
    
    Returns
    -------
    The fractonal coverage of the dye on the molecule
    """
    return (np.sqrt((intensity - caliConst[0])))/caliConst[1];
    

def fitCaliCurve(filename):
    """Fit the calibration curve data the to Langmuir model
    Parameters
    ---------
    
    filename : the name of the file containing the calibration data
    
    Returns
    -------
    constants : the obtained fit constants
        
    corr : The covariance matrix for the fit
    
    """
    data = np.loadtxt(filename + '.txt')
    SHSint= lambda x, B, a, N, K: B+np.square(a*(((x+N+(55.5/K))-np.sqrt(np.square(x+N+(55.5/K))\
                                                   -(4*x*N)))/(2*N)))
    
    constants,corr=sp.curve_fit(SHSint,data[:,0],data[:,1], p0=(1,1,1,1)\
                                ,maxfev=1000000,bounds=([0,0,0,0],[np.inf,np.inf,np.inf,np.inf]))
    return constants, corr, data


def calcCaliCorrandR(constants, corr, data, outName):
    """Calculates the correlation matrix and the R^2 value for the 
    Parameters
    ----------
    constants : the fit constants obtained from fitCaliCurve
    
    corr : the covariance matrix obtained from fitCaliCurve
    
    data : the measured experimental data
    
    outName : the header for the outputted text file
    
    """
    print(constants)
    perr=np.sqrt(np.diag(corr))
    print(perr)
    corrmat=np.zeros([len(constants),len(constants)])
    for i in range(len(corr)):
        for j in range(len(corr)):
            
            ele=corr[i,j]
            diele=ele/(perr[i]*perr[j])
            corrmat[i,j]=round(diele,3)
    print(corrmat)
    #calculate the r^2 value
    ss_res = 0
    ss_total = 0
    residuals = np.zeros([len(data[:,0]), 1])
    for i in range(len(data[:,0])):
        residuals[i] = (LangmuirCurve(data[i,0],constants[0],constants[1],constants[2],constants[3]) - data[i,1])
        ss_res += np.square(residuals[i])
        ss_total += np.square((data[i,1] - np.average(data[:,1])))
    print(ss_res)
    print(ss_total)
    r_sq = 1 - (ss_res/ss_total)
    print(r_sq)
    #write out the fit results
    f = open(outName + "_cali_constants.txt", 'w')
    f.write("B\ta\tN\tK\n")
    for i in range(len(constants)):
        f.write('%.9f' %constants[i] + "\t")
    f.write("\n\n")
    for i in range(len(corr)):
        f.write('%.9f' %perr[i] + "\t")
    f.write("\n\n")
    f.write("Correlation matrix :\n\n")
    for i in range(len(corr)):
        for j in range(len(corr)):
            f.write('%.9f' %corrmat[i,j]+'\t')
        f.write("\n\n")
    f.write("R^2 value : \t" + '%.9f' %r_sq)
    f.close()

                  
def plotCaliCurve(constants, data, outName):
    """Plots the calibration curve data along with the obtained fit curve
    Parameters
    ---------
    constants : the fit constants obtained from fitCaliCurve
        
    data : the measured experimental data
        
    outName : the name of the output file
    
    """
    x=np.linspace(min(data[:,0]),max(data[:,0]),1000)
    plt.figure()
    plt.rcParams.update({'font.size' : 16})
    plt.scatter(data[:,0],data[:,1])
    plt.plot(x,LangmuirCurve(x,constants[0],constants[1],constants[2],constants[3]))
    #plt.xlabel("MG Concentration (nM)")
    #plt.ylabel("Relative SHS signal (Arb. Units)")
    plt.savefig(outName + "_cali_model_plot.png")
    plt.show()



def plotCaliCoverage(constants, data, outName):
    """Plot the fractional coverage of the dye
    Parameters
    ----------
    constants : the fit constants obtained from fitCaliCurve
    
    data : the measured experimental data
    
    outName : the name of the output file
    
    """
    #need to invert the data points and plot them
    expCoverages = IntToCoverage(data[:,1], constants);
    x=np.linspace(min(data[:,0]),max(data[:,0]),1000)
    plt.figure()
    plt.rcParams.update({'font.size' : 16})
    plt.plot(x, Frac_Cov(x,constants[2], constants[3]))
    plt.scatter(data[:,0], expCoverages);
    plt.ylim([0,1])
    plt.savefig(outName + "_cali_coverage_plot.png")
    plt.show()


####Start of methods for displacement curve fitting
def fitDispCurve(fileName, caliConst):
    """Fit the displacement data to the displacement curve
    Parameters
    ---------
    fileName : the file name for the experimental data
    
    caliConst : The constants obtained from fitCaliCurve
    
    Returns
    -------
    aaconst : the constants obtained from a fit to the displacement model
        
    aacorr : the covariance matrix obtained from the fit to the displacement model
    
    """
    AAdata = np.loadtxt(fileName + '.txt')
    k=caliConst[3]
    n=caliConst[2]
    #Be sure to change this appropriately to the fixed dye conc
    x=k*((fixed_dye_conc-n)/55.5)
    n = n
    #here is the code with N fixed to MG's N value
    AAshg=lambda c, B,a,K: B+np.square(a*(x/(1+x+(K*((c-n)/55.5)))))
    '''
    #Not used currently, could be useful in the future
    #More advanced fitting
    C = (fixed_dye_conc * k) / 55.5;
    D = (k * n) / 55.5;
    AAshgadv = lambda x, B, a, K : B + np.square(a * \
                                                    (((((K * x) / 55.5) + C + D) - \
                                                      np.sqrt(np.square((((K * x) / 55.5) + C + D)) -\
                                                              (4 * (((K * n) / 55.5) + D) * C))) \
                                                                /(2 * ((((K * n) / 55.5) + D)))))
    '''
    aaconst, aacorr = sp.curve_fit(AAshg, AAdata[:,0],AAdata[:,1],p0=(0,10,0.0001)\
                                 ,maxfev=1000000,bounds=([0,0,0],[min(AAdata[:,1]), np.inf, np.inf]))
    
    return aaconst, aacorr, AAdata

def calcDispCorrandR(aaconst, aacorr, caliConst, AAdata, outName):
    """Calculate the displacement curve and the R^2 value for the fit
    Parameters
    ----------
    aaconst : the constants obtained from fitDispCurve
        
    aacorr : the covariance matrix obtained from fitDispCurve
        
    caliConst : the constants obtained from fitCaliCurve
    
    AAdata : the experimental data for the displacement data
    
    outName : the output file name header
    
    """
    k=caliConst[3]
    n=caliConst[2]
    #Be sure to change this appropriately to the fixed dye conc
    x=k*((fixed_dye_conc-n)/55.5)
    n = n
    print(aaconst)
    perr2=np.sqrt(np.diag(aacorr))
    print(perr2)
    corrmat2=np.zeros([len(aaconst),len(aaconst)])
    for i in range(len(aacorr)):
        for j in range(len(aacorr)):
            ele=aacorr[i,j]
            diele=ele/(perr2[i]*perr2[j])
            corrmat2[i,j]=round(diele,3)
    print(corrmat2)
    #calculate the r^2 value
    AAss_res = 0
    AAss_total = 0
    residuals = np.zeros([len(AAdata[:,0]), 1])
    for i in range(len(AAdata[:,0])):
        residuals[i] = (DispCurve(AAdata[i,0],x,aaconst[0],aaconst[1],n,aaconst[2]) - AAdata[i,1])
        AAss_res += np.square(residuals[i])
        AAss_total += np.square((AAdata[i,1] - np.average(AAdata[:,1])))
    print(AAss_res)
    print(AAss_total)
    AAr_sq = 1 - (AAss_res/AAss_total)
    print(AAr_sq)
    #write out the fit results
    f = open(outName + "_disp_constants.txt", 'w')
    f.write("B\ta\tN\tK\n")
    for i in range(len(aaconst)):
        f.write('%.9f' %aaconst[i] + "\t")
    f.write("\n")
    for i in range(len(aacorr)):
        f.write('%.9f' %perr2[i] + "\t")
    f.write("\n\n")
    f.write("Correlation matrix :\n\n")
    for i in range(len(aacorr)):
        for j in range(len(aacorr)):
            f.write('%.9f' %corrmat2[i,j]+'\t')
        f.write("\n\n")
    f.write("R^2 value : \t" + '%.9f' %AAr_sq)
    f.close()


def plotDispCurve(aaconst, caliConst, data, outName):
    """Plot the displacement curve data and the model fit
    Parameter
    --------
    aaconst : the constants obtained from fitDispCurve
    
    caliConst : the constants obtained from fitCaliCurve
    
    data : the experimental data for the displacement data
    
    outName : the name of the saved plot
    
    
    """
    k=caliConst[3]
    n=caliConst[2]
    #Be sure to change this appropriately to the fixed dye conc
    x=k*((fixed_dye_conc - n)/55.5)
    n = n 
    z=np.linspace(min(data[:,0]),max(data[:,0]),1000)
    plt.figure()
    plt.rcParams.update({'font.size' : 16})
    plt.scatter(data[:,0]/1000,data[:,1])
    #Change the inputs below depending on whether N is floating or fixed
    plt.plot(z/1000,DispCurve(z,x,aaconst[0],aaconst[1],n,aaconst[2]))
    #Did you change above?
    #plt.xlabel("CTAB Concentration (uM)")
    #plt.ylabel("Relative SHS signal (Arb. Units)")
    plt.ylim([0,max(data[:,1])])
    plt.savefig(outName + '_disp_model_plot.png')
    plt.show()


#want to make a plot of fractional coverage wrt concentration
def plotDispCoverage(dispData, dye_conc, caliConst, dispConst, outName):
    """Plot the surface coverage of the dye and displacing molecule
    Parameters
    ---------
    dispData : The concentrations of the displacing molecule
        
    dye_conc : The fixed dye concentration used
        
    caliConst : The fit constants from the dye calibration from fitCaliCurve()
        
    dispConst : The fit constants from the displacement fit in fitDispCurve()
        
    outName : The output file name header
    
    """
    z = np.linspace(min(dispData[:,0]), max(dispData[:,0]), 1000)
    dispCoverage = DispDyeCoverage(z, dye_conc, caliConst[2], dispConst[2], caliConst[3])
    knownDyeCov = Frac_Cov(dye_conc, caliConst[2], caliConst[3])
    caliDyeCoverage = DispDyeCoverage(dye_conc, z, caliConst[2], caliConst[3], dispConst[2])
    #Covert the exp intensities to coverages
    expCoverages = IntToCoverage(dispData[:,1], dispConst);
    expCoverages = (expCoverages / expCoverages[0]) * knownDyeCov
    plt.figure()
    plt.rcParams.update({'font.size' : 16})
    plt.plot(z/1000, ((caliDyeCoverage) / caliDyeCoverage[0]) * knownDyeCov)
    plt.scatter(dispData[:,0]/1000, expCoverages);
    #plt.ylim([0,1])
    plt.savefig(outName + "_disp_coverage_plot.png")
    plt.show()
    #might as well get R^2 here
    ss_res = 0
    ss_total = 0
    residuals = np.zeros([len(dispData[:,0]), 1])
    for i in range(len(dispData[:,0])):
        predictedDyeCov = DispDyeCoverage(dye_conc, dispData[i,0], caliConst[2], caliConst[3], dispConst[2])
        scaledPredictedDyeCov = (predictedDyeCov / caliDyeCoverage[0]) * knownDyeCov
        residuals[i] = scaledPredictedDyeCov - expCoverages[i]
        ss_res += np.square(residuals[i])
        ss_total += np.square((expCoverages[i] - np.average(expCoverages[:])))
    print(ss_res)
    print(ss_total)
    r_sq = 1 - (ss_res/ss_total)
    print(r_sq)
    


#function to be run externally to batch run files
def extMain(caliDataName, dispDataName):
    #TODO
    return 0

if __name__ == "__main__":
    #fit the calibration data and calculation the quality of the fit
    caliConstants, caliCorr, caliData = fitCaliCurve(caliDataName)
    calcCaliCorrandR(caliConstants, caliCorr, caliData, caliDataName)
    plotCaliCurve(caliConstants, caliData, caliDataName)
    
    #plot the coverage of the dye versus conc
    plotCaliCoverage(caliConstants, caliData, caliDataName)
    
    #fit the displacement curve and calculate the quality of the fit
    dispConstants, dispCorr, dispData = fitDispCurve(dispDataName, caliConstants)
    calcDispCorrandR(dispConstants, dispCorr, caliConstants, dispData, dispDataName)
    plotDispCurve(dispConstants, caliConstants, dispData, dispDataName)
    #plotDispCoverage(dispData, fixed_dye_conc, caliConstants, dispConstants, dispDataName)
    
