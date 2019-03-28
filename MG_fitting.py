#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 11:33:12 2018

Fitting routine for Langmuir isotherm for use with Malchite green SHS data
The isotherm is modified to account for adsorbate bulk concentration change 
resulting from adsorption



@author: wcole
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sp

#import the dataset
#format will be two columns of [Conc., Int.]
data=np.loadtxt('pPSB_NYS_Cali_PP.txt')
#AAdata=np.loadtxt('SNP_AA_disp_60uMMG_pH_8.txt')

xcol=data[:,0]
ycol=data[:,1]
SHSint= lambda x, B, a, N, K: B+np.square(a*(((x+N+(55.5/K))-np.sqrt(np.square(x+N+(55.5/K))\
                                               -(4*x*N)))/(2*N)))

constants,corr=sp.curve_fit(SHSint,data[:,0],data[:,1], p0=(1,1,1,1)\
                            ,maxfev=1000000,bounds=([0,0,0,0],[np.inf,np.inf,np.inf,np.inf]))

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

                  
def LangmuirCurve(x, B,a,N,K):
    return B+np.square(a*(((x+N+(55.5/K))-np.sqrt(np.square(x+N+(55.5/K))-(4*x*N)))/(2*N)))

x=np.linspace(min(data[:,0]),max(data[:,0]),1000)

plt.figure()
plt.rcParams.update({'font.size' : 16})
plt.scatter(data[:,0]*1000,data[:,1])
plt.plot(x*1000,LangmuirCurve(x,constants[0],constants[1],constants[2],constants[3]))
#plt.xlabel("MG Concentration (nM)")
#plt.ylabel("Relative SHS signal (Arb. Units)")
plt.savefig("pPSB_NYS_cali_large.png")
plt.show()





k=constants[3]
n=constants[2]
#Be sure to change this appropriately to the fixed dye conc
x=k*((60-n)/55.5)




#Code for N floating
#AAshg=lambda c, B,a,N,K: B+np.square(a*(x/(1+x+(K*((c-N)/55.5)))))

#aaconst, aacorr=sp.curve_fit(AAshg, AAdata[:,0],AAdata[:,1],p0=(1,1,1,1)\
#                             ,maxfev=100000,bounds=([0,0,0,0],[np.inf,np.inf,np.inf,np.inf]))


#here is the code with N fixed to MG's N value
AAshg=lambda c, B,a,K: B+np.square(a*(x/(1+x+(K*((c-n)/55.5)))))

aaconst, aacorr=sp.curve_fit(AAshg, AAdata[:,0],AAdata[:,1],p0=(1,1,0.1)\
                             ,maxfev=1000000,bounds=([0,0,0],[np.inf,np.inf,np.inf]))


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


def curve2(y,B,a,N,K):
    return B+np.square(a*(x/(1+x+(K*((y-N)/55.5)))))

z=np.linspace(min(AAdata[:,0]),max(AAdata[:,0]),1000)
plt.figure()
plt.scatter(AAdata[:,0],AAdata[:,1])
#Change the inputs below depending on whether N is floating or fixed
plt.plot(z,curve2(z,aaconst[0],aaconst[1],constants[2],aaconst[2]))
#Did you change above?
plt.xlabel("CTAB Concentration (uM)")
plt.ylabel("Relative SHS signal (Arb. Units)")
plt.ylim([0,3500])
#plt.savefig("SNP_1um_60uMMG_AA_disp_pH8_plot.png")
plt.show()
