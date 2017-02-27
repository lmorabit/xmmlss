#! /usr/bin/env python

import numpy
from pylab import *
import math
import pyfits
import scipy.signal.signaltools

import matplotlib
pgf_with_rc_fonts = {"pgf.texsystem": "pdflatex"}
matplotlib.rcParams.update(pgf_with_rc_fonts)

rc('text',usetex=True)
rc('font',family='serif')

rmax=50.0
magmax=19.8

fwhm=15.0 # fwhm of beam
#read in Video catalogue

hdulist=pyfits.open('/Users/Matt1/Work/GMRTdataNEW/InputCatG12.fits')
cat=hdulist[1].data
mags=cat.field('PETROMAG_R')-cat.field('EXTINCTION_R')
ids=cat.field('CATAID')
ra=cat.field('RA')
dec=cat.field('DEC')

#remove stars different fields
#keep=numpy.where((mags<magmax)&(mags>12.2)&(ra < 160))
keep=numpy.where((mags<magmax)&(mags>12.2)&(ra > 160)&(ra <200))
#keep=numpy.where((z>0.00)&(mags<magmax)&(mags>12.2)&(ra >200))

ids=cat.field('CATAID')[keep]
mags=mags[keep]
ra =cat.field('RA')[keep] 
dec =cat.field('DEC')[keep]


#stellar locus colors

#calculate n(m)

bins=numpy.linspace(12.2,magmax,39)
plotbins=numpy.linspace(12.3,magmax-0.1,38)

print 'bins'
print bins
print 'plotbins'
print plotbins


nm=hist(mags,bins=bins,log=True)
cumnm=hist(mags,bins=bins,cumulative=True)

totalsources=cumnm[0][-1]
print totalsources
rarange=ra.max()-ra.min()
decrange=dec.max()-dec.min()

print rarange
print decrange
area=rarange*decrange
print area

asecarea=rarange*decrange*3600.0*3600.0
print asecarea
#asecarea=35.0*3600.0*3600.0
sdensity=totalsources/asecarea
sdensitycirc=sdensity*math.pi*rmax*rmax

print 'Number of sources per square arcsec: \t'+'%.5f'%sdensity
print 'Number of sources per '+'%.2f'%rmax+'\" search radius:'+'\t' +'%.3f'%sdensitycirc

#read in radio source catalogue
radiofits=pyfits.open('/Users/Matt1/Work/GMRTdataNEW/G12MultiFree.fits')
#radiofits=pyfits.open('12hr_cat.fits')
radio=radiofits[1].data

snr=radio.field('S')/radio.field('ErrS')
snrkeep=numpy.where(snr>=5.0)

figure(1)
hist(snr[snrkeep], bins=100,range=(0,400),normed=True)


print 'Number of Source read in ' +numpy.str_(len(snr))
print numpy.str_(len(snrkeep[0]))+' sources match selection criteria'

#radioname=radio.field('col1')[snrkeep]
#print radioname

radioid=radio.field('INDEX')[snrkeep]
radiora=radio.field('RA')[snrkeep]
radiodec=radio.field('DEC')[snrkeep]
radioflux=radio.field('S')[snrkeep]
radiofluxerr=radio.field('ErrS')[snrkeep]


#remove astrometric offsets
#radiora=radiora-6.17421e-5
#radiodec=radiodec-0.0001579109


matchid=numpy.array([])
matchra=numpy.array([])
matchdec=numpy.array([])
matchmag=numpy.array([])
mradioid=numpy.array([])
mradiora=numpy.array([])
mradiodec=numpy.array([])
mradioflux=numpy.array([])
mradiofluxerr=numpy.array([])
mradiospecI=numpy.array([])
nmatchra=numpy.array([])
nmatchdec=numpy.array([])

nomatch=0

#match sources
for i in range(0,len(radioid)):
    #for i in range(112,115):
    radiff=3600.0*numpy.abs(radiora[i]-ra)

    decdiff=3600.0*numpy.abs(radiodec[i]-dec)
    
    rdiff1=(numpy.power(radiff,2)+numpy.power(decdiff,2))
    rdiff=numpy.power(rdiff1,0.5)
    
    roi=numpy.where(rdiff<rmax)
    
    if roi[0].size>0: 
       for p in range(0,roi[0].size):
           matchid=numpy.append(matchid,[ids[roi][p]],axis=0)
           matchra=numpy.append(matchra,[ra[roi][p]],axis=0)
           matchdec=numpy.append(matchdec,[dec[roi][p]],axis=0)
           matchmag=numpy.append(matchmag,[mags[roi][p]],axis=0)
           mradioid=numpy.append(mradioid,[radioid[i]],axis=0)           
           mradiora=numpy.append(mradiora,[radiora[i]],axis=0)
           mradiodec=numpy.append(mradiodec,[radiodec[i]],axis=0)
           mradioflux=numpy.append(mradioflux,[radioflux[i]],axis=0)
           mradiofluxerr=numpy.append(mradiofluxerr,[radiofluxerr[i]],axis=0)
           


    else:
     nomatch=nomatch+1
     nmatchra=numpy.append(nmatchra,[radiora[i]],axis=0)
     nmatchdec=numpy.append(nmatchdec,[radiodec[i]],axis=0)



print numpy.str_(len(matchid)) + ' counterparts within the search radius'
nmatch=len(matchid)
ncentre=len(radioid)
numpy.savetxt('radioid',mradioid)

print numpy.str_(nomatch) + ' sources with no counterparts'

total=hist(matchmag,bins=bins,log=True,cumulative=False)
print 'total', total[0]
#calculate real distribution, qm distribution
nmp=(nm[0]/asecarea)*math.pi*numpy.power(rmax,2)*ncentre
print 'nmp', nmp

real=total[0]-nmp
print 'real', real

numpy.savetxt('realdist',real)
realtotal=numpy.cumsum(real)[-1]
print 'real total', realtotal

qm=real/realtotal
print qm

#q0=(nmatch-sdensity*pi*numpy.power(rmax,2)*ncentre)/ncentre

#from fleuren

q0 = 0.171154

print nmatch
print ncentre
print sdensity

#print 'qcil '+numpy.str_(qcil)
#q0=1-numpy.float_(nomatch)/numpy.float_(ncentre)
#add in q0 value
#q0=0.90
print 'Q0 = '+numpy.str_(q0)
qm1=qm*q0


nm1=(nm[0]/asecarea)
print nm1
#calculate qm/nm
qm_nm=numpy.divide(qm1,nm1)
print 'qm_nm', qm_nm

for i in range(0,21):
  qm_nm[i]=500#qm_nm[18] 
  
print qm_nm
#box=scipy.signal.signaltools.boxcar(3)
#qmsm=scipy.signal.signaltools.convolve(qm1,1/3.0*box,mode='same')
#qmsm[0:18]=qm_nm[0]*nm1[0:18]

#define f(r)

def gaus(r,sigma):
    y=0.5/(math.pi*numpy.power(sigma,2))*numpy.exp(-0.5*numpy.power(numpy.divide(r,sigma),2))
    rtest=numpy.linspace(0,100,2000)
    y1=numpy.multiply(rtest,y)
    norm=numpy.trapz(y1,x=rtest)
    fr=y/(2*math.pi*norm)
    return y
tg=gaus(0.5,0.3)

print 'tg' 
print tg

## create an array of zeros the same size as the number of matches
LR=numpy.zeros([nmatch])
## find the separation
matchsep=numpy.power(numpy.power((matchra-mradiora)*3600.0,2)+numpy.power((matchdec-mradiodec)*3600.0,2),0.5)

## signal to noise ratio
msnr=numpy.divide(mradioflux,mradiofluxerr)


#calculate LR for each match pair
for i in range(0,nmatch):
    
    ## find the magnitude bin where the match is less than nm
    magbin=numpy.where(matchmag[i]<nm[1])
    
    
    #print 'magbin'
    #print magbin
    #print qm_nm
    #print 'nm1', nm[1]
    #print matchmag[i]
    
    qmsmi=qm_nm[magbin[0][0]-1]
    #print matchmag[i], 'qmsmi', qmsmi
    sigma0=max(numpy.power(numpy.power(1.0,2)+numpy.power(0.6*fwhm/msnr[i],2),0.5),0.5)
    fr=max(gaus(r=matchsep[i],sigma=sigma0),0.000001)
    #print i, matchsep[i], msnr[i], sigma0
    LR[i]=qmsmi*fr
    #print 'LR', LR[i]
    #print LR[i]
rel=numpy.zeros([nmatch])


#calculate Rel for each pair
for i in range(0,nmatch):
    idmatch=numpy.where(mradioid[i]==mradioid)
    #print mradioid[i]
    LRarray=LR[idmatch]
    LRsum=numpy.sum(LRarray)
    #print i
    #print 'LR', LR[i]
    #print 'LRsum', LRsum
    rel[i]=LR[i]/(LRsum+(1-q0))
    #print 'REL', rel[i]
    #print q0
    #print LRsum

noreliable=numpy.where(rel>0.8)
xreliable=numpy.where(rel<0.8)
print 'no of reliablesource '+ numpy.str_(len(noreliable[0]))
minus=1.0-rel[noreliable]
print 'minus'
minusum=numpy.sum(minus)
print 'No of cont'
print minusum

#print 'shit percent', minusum/(noreliable[0])

#mfluxrel=mradioflux[noreliable]
#fluxrelhist=hist(mfluxrel,bins=[0.07,0.08,0.09,0.1,0.12,0.14,0.16,0.18,0.2,0.4,0.6,0.8,1,2,4,10,40])
#fluxhist=hist(radioflux,bins=[0.07,0.08,0.09,0.1,0.12,0.14,0.16,0.18,0.2,0.4,0.6,0.8,1,2,4,10,40])

#complete=numpy.float_(mfluxrelhist[0])/numpy.float_(fluxhist[0])
#print mfluxrelhist[0]
#print fluxhist[0]
filename="LikelihoodratioG12FINAL.csv"
file=open(filename,'w')
s='#No \t GAMA_RA \t\t GAMA_DEC \t radio_ra \t radio_dec \t mag \t S \t ErrS \t sep \t radioid \t idno \t LR \t rel \t snr \n'
file.write(s)
for i in range(0,nmatch):
    s=numpy.str_(i+1)+'\t'+'%.6f'%matchra[i]+'\t '+'%.6f'%matchdec[i]+'\t'+'%.6f'%mradiora[i]+'\t '+'%.6f'%mradiodec[i]+'\t'+'%4.6f'%matchmag[i]+'\t'+'%4.6f'%mradioflux[i]+'\t'+'%4.6f'%mradiofluxerr[i]+'\t'+'%.6f'%matchsep[i]+'\t'+'%i'%mradioid[i]+'\t'+'%i'%matchid[i]+'\t\t'+'%4.4f'%LR[i]+'\t'+'%.4f'%rel[i]+'\t'+'%.4f'%msnr[i]+ '\n'
    file.write(s)
file.close

file=open('reliableG12FINAL.reg','w')
for i in range(0,len(noreliable[0])):
    s='fk5; circle('+'%.6f'%mradiora[noreliable[0][i]]+','+'%.6f'%mradiodec[noreliable[0][i]]+',2.8")# color=green'+'\n'
    file.write(s)
file.close    

file=open('nonreliableG12FINAL.reg','w')
for i in range(0,len(xreliable[0])):
    s='fk5; circle('+'%.6f'%mradiora[xreliable[0][i]]+','+'%.6f'%mradiodec[xreliable[0][i]]+',2.75")# color=magenta'+'\n'
    file.write(s)
file.close  

file=open('nomatchG12FINAL.reg','w')
for i in range(0,len(nmatchra)):
    s='fk5; circle('+'%.6f'%nmatchra[i]+','+'%.6f'%nmatchdec[i]+',2.8")# color=red'+'\n'
    file.write(s)
file.close


figure(2)
ax1=subplot(311)
semilogy(plotbins,nmp,ls='steps-',color='k',label='Background')
semilogy(plotbins,total[0],ls='steps-',color='0.55',label='Total')
semilogy(plotbins,real,ls='steps--',color='k',label='Real')


leg=legend(loc='lower right')
frame=leg.get_frame()
frame.set_edgecolor('1.0')
ylim(ymin=0.1,ymax=400)
ylabel('N(counterparts)')
#xlabel('VIDEO Kmag')
xlim(xmin=12,xmax=21.0)
#ax1.set_xticklabels([])
ax1.set_yticklabels(['','$10^{0}$','$10^{1}$','$10^{2}$'])
ax2=subplot(312)
semilogy(plotbins,qm_nm,ls='steps-',color='k')
ylabel('P(m)=q(m)/n(m)')
#xlabel('VIDEO Kmag')
xlim(xmin=12,xmax=21.0)
#ax2.set_xticklabels([])
ax2.set_yticklabels(['','$10^2$','$10^3$','$10^4$'])
ax3=subplot(313)
#ax3.set_xticklabels(['$12$','$14$','$16$','$18$','$20$'])
ax3.set_yticklabels(['$10^{-4}$','$10^{-3}$','$10^{-2}$','$10^{-1}$','$10^0$'])
semilogy(plotbins,qm,ls='steps-',color='k',label='raw')
#semilogy(plotbins,qmsm,ls='steps-',label='best guess')
#ax3.legend(loc='lower left')
xlim(xmin=12,xmax=21.0)
xlabel('R mag')
ylabel('q(m)')
savefig('qmdistributionG12FINAL.pdf',orientation='portrait', papertype='a4')
show()
figure(4)
plot(3600.0*(matchra-mradiora),3600.0*(matchdec-mradiodec),'+')
savefig('match_positionsG12FINAL.pdf',orientation='portrait')


import scipy.stats

print 'meanra'
print numpy.mean(3600.0*(matchra[noreliable]-mradiora[noreliable]))
print numpy.mean(matchra[noreliable]-mradiora[noreliable])
print numpy.std(matchra-mradiora)/numpy.power(len(matchra),0.5)

print 'meandec'
print numpy.mean(3600.0*(matchdec-mradiodec)[noreliable[0]])
print numpy.mean((matchdec-mradiodec)[noreliable])
print numpy.std(matchdec-mradiodec)/numpy.power(len(matchdec),0.5)

figure(5)
subplot(311)
semilogx(LR,rel,'.',color='k',ms=3)
xlim(xmin=0.00001,xmax=1000)
ylim(ymin=-0.1,ymax=1.1)
xlabel('Likelihood Ratio')
ylabel('Reliability')

subplot(312, xscale='log')
hist(LR,bins=(100),range=(1,1000),log=True,histtype='step',color='k',linewidth=2)
xlim(xmin=1.0,xmax=1000)
ylim(ymin=0.9,ymax=500)
xlabel('Likelihood Ratio')
ylabel('N(counterparts)')

subplot(313)
hist(rel,bins=100,log=True,histtype='step',color='k',linewidth=2)
xlabel('Reliability')
xlim(xmin=-0.1,xmax=1.1)
ylabel('N(counterparts)')
savefig('liklihoodG12FINAL.pdf',orientation='portrait')


figure(6)
subplot(211)
hist(3600.0*(matchra-mradiora)[noreliable],bins=10)

subplot(212)
hist(3600.0*(matchdec-mradiodec)[noreliable],bins=10)

show()
