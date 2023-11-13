"""
Created 21. November 2022 by Daniel Van Opdenbosch, Technical University of Munich

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. It is distributed without any warranty or implied warranty of merchantability or fitness for a particular purpose. See the GNU general public license for more details: <http://www.gnu.org/licenses/>
"""

import os
import glob
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import stats

files=glob.glob('*.csv')

for f in files:
	filename=os.path.splitext(f)[0]
	Po,P0o,Va_ccpergo=np.genfromtxt(f,unpack=True,delimiter=';',skip_header=1,usecols=(1,5,15),encoding='iso-8859-1')

	Adsargs=np.arange(0,np.where(np.isnan(Po)==True)[0][0])
	Bothargs=np.where(np.isnan(Po)==False)
	Desargs=np.arange(np.where(np.isnan(Po)==True)[0][0],len(Po))

	args=Adsargs
	P,P0,Va_ccperg=Po[args],P0o[args],Va_ccpergo[args]

	BETargs=np.where((P/P0>=0.05)&(P/P0<=0.3))
	slope,intercept,r,p,se=stats.linregress((P/P0)[BETargs],(P/(Va_ccperg*(P0-P)))[BETargs])

	BETsurface=4.35/(slope+intercept)

	plt.close('all')
	mpl.rc('text',usetex=True)
	mpl.rc('text.latex',preamble=r'\usepackage[helvet]{sfmath}')
	fig,ax1=plt.subplots(figsize=(7.5/2.54,5.3/2.54))

	ax1.plot((P/P0)[BETargs],(P/(Va_ccperg*(P0-P)))[BETargs],c='k',lw=0  ,marker='s',ms=2,mew=0)
	ax1.plot((P/P0)[BETargs],(P/P0*slope+intercept)[BETargs],c='k',lw=0.5,marker='s',ms=0,mew=0)

	ax1.text(ax1.get_xlim()[0]+(ax1.get_xlim()[1]-ax1.get_xlim()[0])/20,ax1.get_ylim()[1]-(ax1.get_ylim()[1]-ax1.get_ylim()[0])*0.2,r'$\rm{'+filename.replace('_','-')+'}$'+'\n'+r'$\rm{BET\ surface:\ }'+str(BETsurface.round(2))+r'\ \rm{m}^2/\rm{g}$',fontsize=8)

	ax1.set_xlabel(r'$P/P_0/1$',fontsize=10)
	ax1.set_ylabel(r'$P/(V_{\rm{a}}(P_0-P))/(\rm{g}\,\rm{cm}^{-3})$',fontsize=10)
	ax1.tick_params(axis='both',pad=2,labelsize=8)
	plt.tight_layout(pad=0.1)
	plt.savefig(filename+'_BETplot.pdf',transparent=True)
	plt.savefig(filename+'_BETplot.png',dpi=300)

	args=Adsargs
	P,P0,Va_ccperg=Po[args],P0o[args],Va_ccpergo[args]

	t=0.88*(P/P0)**2+6.45*(P/P0)+2.98

	STSAargs=np.where((P/P0>=0.2)&(P/P0<=0.5))
	slope,intercept,r,p,se=stats.linregress(t[STSAargs],Va_ccperg[STSAargs])

	STSA=slope*15.47

	plt.close('all')
	mpl.rc('text',usetex=True)
	mpl.rc('text.latex',preamble=r'\usepackage[helvet]{sfmath}')
	fig,ax1=plt.subplots(figsize=(7.5/2.54,5.3/2.54))

	ax1.plot(t[STSAargs],Va_ccperg[STSAargs]		  ,c='k',lw=0  ,marker='s',ms=2,mew=0)
	ax1.plot(t[STSAargs],(t*slope+intercept)[STSAargs],c='k',lw=0.5,marker='s',ms=0,mew=0)

	ax1.text(ax1.get_xlim()[0]+(ax1.get_xlim()[1]-ax1.get_xlim()[0])/20,ax1.get_ylim()[1]-(ax1.get_ylim()[1]-ax1.get_ylim()[0])*0.2,r'$\rm{'+filename.replace('_','-')+'}$'+'\n'+r'$\rm{STSA:\ }'+str(STSA.round(2))+r'\ \rm{m}^2/\rm{g}$',fontsize=8)

	ax1.set_xlabel(r'$t/\rm{\AA}$',fontsize=10)
	ax1.set_ylabel(r'$V_{\rm{a}}/(\rm{cm}^3\,\rm{g}^{-1})$',fontsize=10)
	ax1.tick_params(axis='both',pad=2,labelsize=8)
	plt.tight_layout(pad=0.1)
	plt.savefig(filename+'_STSAplot.pdf',transparent=True)
	plt.savefig(filename+'_STSAplot.png',dpi=300)

	args=Adsargs
	P,P0,Va_ccperg=Po[args],P0o[args],Va_ccpergo[args]

	rk=-4.14/np.log(P/P0)
	Vp=np.gradient(Va_ccperg)/np.gradient(np.log(rk))*0.02883/22.71

	BJHargs=np.where((np.isnan(Vp)==False)&(Vp!=0)&(rk>0))

	plt.close('all')
	mpl.rc('text',usetex=True)
	mpl.rc('text.latex',preamble=r'\usepackage[helvet]{sfmath}')
	fig,ax1=plt.subplots(figsize=(7.5/2.54,5.3/2.54))

	ax1.plot(rk[BJHargs],Vp[BJHargs],c='k',lw=0.5,marker='s',ms=2,mew=0)

	ax1.set_xscale('log')

	ax1.text(ax1.get_xlim()[0]*2,ax1.get_ylim()[1]-(ax1.get_ylim()[1]-ax1.get_ylim()[0])*0.1,r'$\rm{'+filename.replace('_','-')+'}$',fontsize=8)

	ax1.set_xlabel(r'$r_{\rm{p}}/\rm{\AA}$',fontsize=10)
	ax1.set_ylabel(r'$\rm{d}V_{\rm{a}}/\rm{d}\log{r_{\rm{p}}}/(\rm{cm}^3\,\rm{g}^{-1})$',fontsize=10)
	ax1.tick_params(axis='both',pad=2,labelsize=8)
	plt.tight_layout(pad=0.1)
	plt.savefig(filename+'_BJHplot.pdf',transparent=True)
	plt.savefig(filename+'_BJHplot.png',dpi=300)
