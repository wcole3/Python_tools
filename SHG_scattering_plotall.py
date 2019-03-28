#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 11:45:04 2017
Python script to import all txt files from 'batch' folder in directory
Then fit the spectra in each to a gaussian before plotting (w/gaussian fit)
and outputing data about the files including the area of the peak and the data
file to which it contains in a txt document called "Fit_results_"+name (set by 
user).

                                                                        
Most recent change Date: March 5 2018
Altered the removeCR function to simply remove points that are affected by 
cosmic rays.
                                                                        
Author email: will.cole3@gmail.com
@author: wcole
"""

import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
import glob
import os

#A function to fill the fit result folder with useful info
def popHeaderFit(f):
    """Generate the header for the output file
        
    Parameters
    ----------
    f : name of the output file
        
    """
    f.write("Results for the Gaussian fits of SHG scattering experiments\n")
    f.write("Results are organized as follows:\n")
    f.write("File name of data\n"+
            "Norm (a.u.)  \t####\tSD\n"
            "Mean (nm)    \t####\tSD\n"+
            "Sigma (nm)   \t####\tSD\n"+
            "Slope (a.u.) \t####\tSD\n"+
            "Offset (a.u.)\t####\tSD\n\n"+
            "Correlation Matrix\n\n"+
            "Additional Notes are listed below\n\n\n")
#Another funtion to populate teh header of the peak area file
def popHeaderPA(f):
    """Generate the header for the output file
        
    Parameters
    ----------
    f : output file name
        
    """
    f.write("Peaks areas of the 400nm SH signal formatted as:\n"+
            "Peak Area\t\tFilename\n")
#Final write out function to spew out the fit results info
def writeOutFit(dataName, valueArr, covarArr, notes):
    """Write the results of the fits to a file
        
    Parameters
    ----------
    dataName : name of the data set, this will be appended to a file name header
    
    valueArr : the vector representing the fit constants obtained from gaussFit
    
    covarArr : the covarviance matrix from gaussFit
    
    notes : any user defined notes
        
    """
    # create sigmas
    sigma=np.sqrt(np.diag(covarArr))
    #This adheres to the format in popHeaderFit
    lineFormat=["Norm (a.u.)  \t","Mean (nm)    \t","Sigma (nm)   \t",
                "Slope (a.u.) \t","Offset (a.u.)\t"]
    outString="Data filename: "+dataName+'\n\n'
    for i in range(len(valueArr)):
        outString+=lineFormat[i]+'%s'%valueArr[i]+'\t'+'%s'%sigma[i]+'\n'
    outString+='\n'
    #now for the Correlation matrix
    for i in range(len(covarArr)):
        for j in range(len(covarArr)):
            if i<=j:
                ele=covarArr[i,j]
                diele=ele/(sigma[i]*sigma[j])
                corrMatElem=round(diele,3)
                outString+='%s'%corrMatElem+'\t'
        outString+='\n'
        for n in range(i+1):
            outString+='\t'
    outString+='\n'
    outString+="Additional Notes: \n"
    outString+=notes
    outString+='\n\n'
    return outString
#ALSO add in functionallity to remove all data points that fall outside the 
#spectral range of interest 385-415 nm
def sliceArr(inArr):
    """Limit the range fed to gaussFit
        
    Parameters
    ---------
    inArr : the input data array
    
    Returns
    ------
    retArr : the array with the desired range
        
    """
    retArr=np.zeros([0,2])
    for i in range(len(inArr)):
        if 394<=inArr[i,0]<=406 and 0 < inArr[i,1]:
            retArr=np.append(retArr,[inArr[i,:]],axis=0)
    return retArr        
#A function to automate removal of points from cosmic ray impacts on CCD
def removeCR(inArr):
    """Remove any cosmic rays that snuck into the data
        
    Parameters
    ----------
    inArr : the data to search for cosmic rays, generally this is any data whose value is 100x
        greater that the previous entry (spikes in the data)
        
    Returns
    -------
    outArr : the input array with any cosmic rays removed
    
    note : a string describing the points removed
        
    """
    note="none"
    tally=0
    outArr=np.zeros([0,2])
    
    #print(avgDiff, stDiff)
    
    for i in range(1,len(inArr)):
        #this loop doesnt work if the CR is at the first data point
        #There is a special case after this loop to attempt to take care of it
        check=(inArr[i,1]-inArr[i-1,1])
        if 0.5*inArr[i-1,1] <= check and 100 <= inArr[i,1]:
            tally+=1
            inArr[i,1]=inArr[i-2,1]
        if check<0 and 0.5*inArr[i-1,1] <= np.abs(check) and 100 <= np.abs(check):
            tally+=1
            inArr[i,1]=inArr[i-2,1]
        else:
            outArr=np.append(outArr,[inArr[i,:]],axis=0)
    
    if 0 < tally:
        note="Points removed by removeCR = "+ str(tally)
    return outArr,note
#A funrtion to fit the data to a gaussian and return relevant values
def gaussFit(inArr):
    """A method to fit the input data to a typical Gaussian with a slope background:
           y = a * exp( -1 *(x - b)^2 / (2 * c^2)) + (d * (x -e))
        
    Parameters
    ---------
    inArr : the data to be fit to the Gaussian; column 0 is the x values in the equation above,
        column 1 is the y values
        
    Returns
    -------
    const : the fit constants: (a, b, c, d, and e) from the equation above
    
    error : the covariance matrix for the fit
    
    peakArea : the peak area of the Gaussian obtain by integration of the above equation
    
    fitWork : bool; true if the fit worked
        
    """
    fitWork=True
    x=inArr[:,0]
    y=inArr[:,1]
    gaussian= lambda x,a,b,c,d,e: (a*np.exp(-((x-b)**2)/(2*(c**2))))+(d*(x-e))
    try:
        const, error=opt.curve_fit(gaussian, x,y,p0=(1,400,1,1,1),max_nfev=1000000,
                                   bounds=((0,398,0,0,0),(np.inf,402,10,100,1000)))
    except RuntimeError:
        print("Good fit not found")
        const=np.zeros([5,1])
        error=np.zeros([5,5])
        fitWork=False
    
    peakArea=const[0]*const[2]*np.sqrt(2*np.pi)
    return const, error, peakArea, fitWork

#plot the data with the Gaussian fit overlaid
def plotOverlay(inArr, value, filename):
    """Plot the data along with the fit gaussian
        
    Parameters
    ---------
    inArr : the data to plot, x values in column 0 and y values in column 1
    
    value : the constants obtained from gaussFit
    
    filename : the name of the output file
        
    """
    #want an info string for the Gaussian fit constants
    boxString="Mean: "+'%s'%round(value[1],2)+" nm"+"\nSigma: "\
                        +'%s'%round(value[2],2)+" nm"+"\nNorm: "+'%s'%round(value[0],2)
    x=np.linspace(inArr[0,0],inArr[len(inArr)-1,0],1000)
    gauss=(value[0]*np.exp(-((x-value[1])**2)/(2*(value[2]**2))))+\
          (value[3]*(x-value[4]))
    plt.figure(dpi=200)
    plt.scatter(inArr[:,0],inArr[:,1],label="Raw Data", marker='.', s=4, color='k')
    plt.plot(x,gauss,color='r',label="Gaussian Fit")
    plt.xlabel("Wavenumber (nm)")
    plt.ylabel("Intensity (Astronomical Units)")
    plt.text(inArr[0,0],max(inArr[:,1]),boxString, bbox=dict(facecolor='white'))
    plt.legend(loc='upper right')
    plt.savefig(filename+'.png')
    plt.close()

def main():
    #start by getting the working directory, and setting the batch and bin folder
    work_dir=os.path.dirname(os.path.realpath(__file__))
    batch_dir=work_dir+'/batch'
    bin_dir=work_dir+'/bin'
    dtype=[('Filename',list),('FitResults',list),('PeakArea',float)]
    fitList=[]
    
    #open the file used to collect results
    f=open(bin_dir+"/Fit_results.txt", 'w')
    popHeaderFit(f)
    g=open(bin_dir+"/Peak_areas.txt",'w')
    popHeaderPA(g)
    #first start pulling in files, THIS IS A RANDOM PULL
    all_files=glob.glob(batch_dir+'/*.txt')
    for file_path in all_files:
        filename=os.path.splitext(os.path.basename(file_path))[0]
        inArr=np.loadtxt(file_path)
        #for each file we want to filter for range and cosmic rays and then fit
        data=sliceArr(inArr)
        data,crNote=removeCR(data)
        value, covar, PA, test= gaussFit(data)
        #Then plot the overlay
        if test:
            plotOverlay(data,value,bin_dir+'/'+filename)
        #Collect results here and then sort then writeout
        fitList.append((filename,writeOutFit(filename,value,covar,crNote),PA))
    ftList=np.array(fitList,dtype=dtype)
    ftList=np.sort(ftList,order='Filename')
    #writeout results for data file
    for i in range(len(ftList)):
        f.write('%s'%ftList['FitResults'][i])
        g.write('%s'%ftList['PeakArea'][i]+'\t\t'+'%s' %ftList['Filename'][i]+'\n')
    
    f.close()
    g.close()
        


if __name__ == "__main__":
    main()
