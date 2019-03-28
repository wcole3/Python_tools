#will cole, 2015, will.cole3@gmail.com
#least square fitter

import numpy as np
import scipy.optimize as opt
name='help'
#format is exp trans,int, Jup, Kup, Jlow, Klow
fitmat=np.loadtxt(name+'.txt')

ycol=fitmat[:,0]
xmat=np.zeros([len(fitmat),4])
xmat[:,0]=fitmat[:,2]
xmat[:,1]=fitmat[:,3]-1
xmat[:,2]=fitmat[:,4]
xmat[:,3]=fitmat[:,5]
#xmat[:,4]=fitmat[:,9]

#print ycol
#print xmat

xmatfit=np.matrix.transpose(xmat)
ycolfit=np.matrix.transpose(ycol)
#for mw
#ergE=lambda x, Bu, Cu,Dju,Del:((Bu*x[0]*(x[0]+1))+((Cu-Bu)*(x[1]**2))-(Dju*((x[0]*(x[0]+1))**2))+(Del*((x[0]*(x[0]+1))**3)))-((Bu*x[2]*(x[2]+1))+((Cu-Bu)*(x[3]**2))-(Dju*((x[2]*(x[2]+1))**2))+(Del*((x[2]*(x[2]+1))**3)))

#save this for perp
#ergE=lambda x, v,Bu, Cu,Dju,Djku,Dku,Bl,Cl,Djl,Djkl,Dkl:(v+((Bu*x[0]*(x[0]+1))+((Cu-Bu)*(x[1]**2))-(Dju*((x[0]*(x[0]+1))**2))-(Djku*x[0]*(x[0]+1)*(x[1]**2))-(Dku*(x[1]**4))))-((Bl*x[2]*(x[2]+1))+((Cl-Bl)*(x[3]**2))-(Djl*((x[2]*(x[2]+1))**2))-(Djkl*x[2]*(x[2]+1)*(x[3]**2))-(Dkl*(x[3]**4)))
#for parallel
ergE=lambda x, v,Bu, Cu,Dju,Djku,Bl,Djl,Djkl:(v+((Bu*x[0]*(x[0]+1))+((Cu)*x[1]**2)+(Dju*(x[0]*(x[0]+1))**2))+(Djku*((x[0]*(x[0]+1))*(x[1]**2))))-((Bl*x[2]*(x[2]+1))+(Djl*(x[2]*(x[2]+1))**2)+(Djkl*((x[2]*(x[2]+1))*(x[3]**2))))
#for asym parallel
#ergE=lambda x, v,Bu, Cu,Dju,Djku, Hju,Hjku,Hkju,Bl,Djl,Djkl,Hjl,Hjkl,Hkjl:(v+((Bu*x[0]*(x[0]+1))+((Cu)*x[1]**2)-(Dju*(x[0]*(x[0]+1))**2))-(Djku*((x[0]*(x[0]+1))*(x[1]**2)))+((Hju*((x[0])*((x[0]+1))))**3)+(Hjku*((((x[0])*((x[0]+1)))**2)*(x[1]**2)))+(Hkju*(x[0]*(x[0]+1))*(x[1]**4)))-((Bl*x[2]*(x[2]+1))-(Djl*(x[2]*(x[2]+1))**2)-(Djkl*((x[2]*(x[2]+1))*(x[3]**2)))+((Hjl*((x[2])*((x[2]+1))))**3)+(Hjkl*((((x[2])*((x[2]+1)))**2)*(x[3]**2)))+(Hkjl*(x[2]*(x[2]+1))*(x[3]**4)))
solve,suc=opt.curve_fit(ergE,xmatfit,ycol[:])
#solve,suc=opt.curve_fit(ergE,xmatfit,ycol[:])
print(solve)
perr=np.sqrt(np.diag(suc))
print(perr)


difflist=[]
for i in range(len(xmat)):
    #need to calculate the predicted transitions
    uppJ=xmat[i,0]
    uppK=xmat[i,1]
    lowJ=xmat[i,2]
    lowK=xmat[i,3]
    #uppl=xmat[i,4]
    #lowl=xmat[i,4]
    
    #below for mw
    '''
    Bu=solve[0]
    Cu=solve[1]
    Dju=solve[2]
    Del=solve[3]
    predicted=((Bu*uppJ*(uppJ+1))+((Cu-Bu)*uppK**2)-(Dju*(uppJ*(uppJ+1))**2)+(Del*((uppJ*(uppJ+1))**3)))-((Bu*lowJ*(lowJ+1))+((Cu-Bu)*lowK**2)-(Dju*(lowJ*(lowJ+1))**2)+(Del*((lowJ*(lowJ+1))**3)))
    '''
    #below is for perpendicular
    '''
    v=solve[0]
    Bu=solve[1]
    Cu=solve[2]
    Dju=solve[3]
    Djku=solve[4]
    Dku=solve[5]
    #Coru=solve[6]
    Bl=solve[6]
    Cl=solve[7]
    Djl=solve[8]
    Djkl=solve[9]
    Dkl=solve[10]
    #Corl=solve[12]
    predicted=((v+((Bu*uppJ*(uppJ+1))+((Cu-Bu)*uppK**2)-(Dju*(uppJ*(uppJ+1))**2)-(Djku*uppJ*(uppJ+1)*(uppK**2))-(Dku*(uppK**4))))-((Bl*lowJ*(lowJ+1))+((Cl-Bl)*lowK**2)-(Djl*(lowJ*(lowJ+1))**2)-(Djkl*lowJ*(lowJ+1)*(lowK**2))-(Dkl*(lowK**4))))
    '''
    #below is for parallel type
    
    v=solve[0]
    Bu=solve[1]
    C=solve[2]
    Dju=solve[3]
    Djku=solve[4]
    Bl=solve[5]
    Djl=solve[6]
    Djkl=solve[7]
    predicted=(v+(Bu*uppJ*(uppJ+1))+(C*(uppK**2))+(Dju*(uppJ*(uppJ+1))**2)+(Djku*(uppJ*(uppJ+1))*(uppK**2)))-((Bl*lowJ*(lowJ+1))+(Djl*(lowJ*(lowJ+1))**2)+(Djkl*(lowJ*(lowJ+1))*(lowK**2)))
    
    #below for asym parallel
    '''
    v=solve[0]
    Bu=solve[1]
    C=solve[2]
    Dju=solve[3]
    Djku=solve[4]
    Hju=solve[5]
    Hjku=solve[6]
    Bl=solve[7]
    Djl=solve[8]
    Djkl=solve[9]
    Hjl=solve[10]
    Hjkl=solve[11]
    predicted=(v+(Bu*uppJ*(uppJ+1))+(C*(uppK**2))-(Dju*(uppJ*(uppJ+1))**2)-(Djku*(uppJ*(uppJ+1))*(uppK**2))+(Hju*(uppJ**3)*((uppJ+1)**3))+(Hjku*((uppJ**2)*((uppJ+1)**2)*(uppK**2))))-((Bl*lowJ*(lowJ+1))-(Djl*(lowJ*(lowJ+1))**2)-(Djkl*(lowJ*(lowJ+1))*(lowK**2))+(Hjl*(lowJ**3)*((lowJ+1)**3))+(Hjkl*((lowJ**2)*((lowJ+1)**2)*(lowK**2))))
    '''
    diff=round(((fitmat[i,0])-predicted),5)
    difflist.append(diff)
print(np.std(difflist[:]))


    




corrmat=np.zeros([len(solve),len(solve)])
for i in range(len(suc)):
    for j in range(len(suc)):
        
        ele=suc[i,j]
        diele=ele/(perr[i]*perr[j])
        corrmat[i,j]=round(diele,3)
#print corrmat

f=open('corr_mat_'+name+'.txt','w')
for i in range(len(solve)):
    f.write('%s' %solve[i]+'\t'+'\t')
    f.write('%s' %perr[i])
    f.write('\n')
f.write('\n')
f.write('\n')
for i in range(len(suc)):
    for j in range(len(suc)):
        f.write('%s' %corrmat[i,j]+'\t')
    f.write('\n')
f.close()

f=open('fit_results_'+name+'.txt','w')
f.write('Exp Line'+'\t'+'Diff'+'\t'+'Int'+'\t'+'UppJ'+'\t'+'UppK'+'\t'+'LowJ'+'\t'+'LowK'+'\t\n')
for i in range(len(fitmat)):
    f.write('%s' %fitmat[i,0]+'\t'+'%s'%difflist[i]+'\t'+'%s'% fitmat[i,1]+'\t')
    for j in range(0,4):
        f.write('%s'%xmat[i,j]+'\t')
    f.write('\n')
f.close()

