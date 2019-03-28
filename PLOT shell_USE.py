#Will Cole, 2015, will.cole3@gmail.com

#generic plotting shell

import numpy as np
import matplotlib.pyplot as plt
#self is the filename with extension
#xdata is column of x values, ydata is column of yaxis
#xlim scales xaxis, pt is type of plot, lab is lab, titl is title, 


def plot(self, xdata,ydata,xlim, pt,lab,titl,signal,bck,sd):
    """Plot a given set of data with some optional data manipulation
        
    Parameters
    ----------
    self : name of the input txt file and name of the output file
    
    xdata : the x data axis of the input txt file to be plotted
    
    ydata : the y data axis of the input txt file to be plotted
    
    xlim : tuple, the x range limit
    
    pt : the type of plot to perform, the options are:
    
            "all" : plot all of the data regardless of overlap
            
            "single" : plot only the first scan
            
            "average" : plot the average of all of the scans
            
            "errorbar" : plot the average with errorbars
    
    lab : the label of the data
    
    titl : the title of the plot
    
    signal : bool; if yes plots the signal seperately from the background
    
    bck : bool; if yes plots the background seperately from the signal
    
    sd : bool; if yes plots the SD as a function of x
    
        
    """
    data=np.loadtxt(self+'.iv')
    slic=len(data)
    for i in range(len(data)):
        if 0 <i:
            if data[i,0] == data[0,0]:
                slic=i
                break
    run=0
    for i in range(len(data)):
        if data[i,0]==data[0,0]:
            run+=1
        
    if pt is "all":
        plt.figure()
        plt.plot((((data[:,xdata]))),(data[:,ydata]), label=lab)
        plt.xlabel("Voltage applied to QCL (V)")
        plt.ylabel("Intensity (a.u.)")
        plt.title(titl,size=10)
        upper=max(data[:,xdata])
        lower=min(data[:,xdata])
        plt.xlim([lower,upper])
        if xlim is not "":
            plt.xlim(xlim)
        if lab is not "":
            plt.legend(loc="best")
        plt.show()
        plt.savefig(self+'.png')
    if pt is "single":
        plt.figure()
        #for i in range(len(data)):
        plt.plot(data[0:slic,xdata],data[0:slic,ydata],label=lab)
        plt.xlabel("Voltage (V)")
        plt.ylabel("Intensity (AU)")
        plt.title(titl,size=10)
        upper=max(data[:,xdata])
        lower=min(data[:,xdata])
        plt.xlim([lower,upper])
        if xlim is not "":
            plt.xlim(xlim)
        if lab is not "":
            plt.legend(loc="best")
        plt.show()
        plt.savefig(self+'.png')
    if pt is "average":
        plt.figure()
        runave=[]
        for i in range(0,slic):
            su=0
            for j in range(0,run):
                su+=data[(i+(j*slic)),ydata]
            runave.append(su)
        array=np.asarray(runave)/3
        #print len(array)
        plt.plot(data[0:slic,xdata],array, label=lab)
        plt.xlabel("Voltage (V)")
        plt.ylabel("Intensity (AU)")
        plt.title(titl,size=10)
        upper=max(data[:,xdata])
        lower=min(data[:,xdata])
        plt.xlim([lower,upper])
        if xlim is not "":
            plt.xlim(xlim)
        if lab is not "":
            plt.legend(loc="best")
        plt.show()
        plt.savefig(self+'.png')
    if pt is "errorbar":
        plt.figure()
        runave=[]
        for i in range(0,slic):
            su=0
            for j in range(0,run):
                su+=data[(i+(j*slic)),ydata]
            runave.append(su)
        array=np.asarray(runave)/3
        sdlist=[]
        for i in range(0, slic):
            su=0
            for j in range(0,run):
                su+=np.power((data[i+(j*slic),3]-array[i]),2)
            sdlist.append(su)
        sarray=np.power(np.asarray(sdlist)/(run-1),0.5)
        plt.errorbar(data[0:slic,xdata],array, yerr=sarray, label=lab)
        plt.xlabel("Voltage (V)")
        plt.ylabel("Intensity (AU)")
        plt.title(titl,size=10)
        upper=max(data[:,xdata])
        lower=min(data[:,xdata])
        plt.xlim([lower,upper])
        if xlim is not "":
            plt.xlim(xlim)
        if lab is not "":
            plt.legend(loc="best")
        plt.show()
        plt.savefig(self+'.png')
    if signal is 'yes':
        plt.figure()
        
        plt.plot(data[0:slic,xdata],data[0:slic,(ydata+2)],label="Signal:"+lab)
        plt.xlabel("Voltage (V)")
        plt.ylabel("Intensity (AU)")
        plt.title("Signal:"+titl,size=10)
        upper=max(data[:,xdata])
        lower=min(data[:,xdata])
        plt.xlim([lower,upper])
        if xlim is not "":
            plt.xlim(xlim)
        if lab is not "":
            plt.legend(loc="best")
        plt.show()
        plt.savefig(self+'_signal'+'.png')
   
    if bck is 'yes':
        plt.figure()
        
        plt.plot(data[0:slic,xdata],(data[0:slic,(ydata+2)]-data[0:slic,ydata]),label="Bck:"+lab)
        plt.xlabel("Voltage (V)")
        plt.ylabel("Intensity (AU)")
        plt.title("Bck:"+titl,size=10)
        upper=max(data[:,xdata])
        lower=min(data[:,xdata])
        plt.xlim([lower,upper])
        if xlim is not "":
            plt.xlim(xlim)
        if lab is not "":
            plt.legend(loc="best")
        plt.show()
        plt.savefig(self+'_bck'+'.png')
    if sd is 'yes':
        plt.figure()
        plt.plot(data[0:slic,xdata],(data[0:slic,(ydata+1)]),label="SD:"+lab)
        plt.xlabel("Voltage (V)")
        plt.ylabel("Intensity (AU)")
        plt.title("SD:"+titl,size=10)
        upper=max(data[:,xdata])
        lower=min(data[:,xdata])
        plt.xlim([lower,upper])
        if xlim is not "":
            plt.xlim(xlim)
        if lab is not "":
            plt.legend(loc="best")
        plt.show()
        plt.savefig(self+'_sd'+'.png')
#self(remove .iv), xdata, ydata, xlim, type, label, title, plt signal?, plt bck?
#the final two are just for bcksub'd scans to pull out the individual signal and bck scans

plot('2f_bck_sub_test_blank_051517', 1, 3,"", "all", "",'background 2f, test', 'yes', 'yes','yes')







        
    
