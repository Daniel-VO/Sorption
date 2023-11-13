"""
Created 18. August 2023 by Daniel Van Opdenbosch, Technical University of Munich

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. It is distributed without any warranty or implied warranty of merchantability or fitness for a particular purpose. See the GNU general public license for more details: <http://www.gnu.org/licenses/>
"""

import os
import sys
import glob
import numpy as np
import pygaps
import matplotlib as mpl
import matplotlib.pyplot as plt
import pygaps.characterisation as pgc
import pygaps.modelling as pgm
from pygaps.graphing.calc_graphs import psd_plot

def carbon_model(relative_p):
	return 0.88*(relative_p**2)+6.45*relative_p+2.98

files=glob.glob('*.csv')
results=[]
models=[]

for f in files:
	filename=os.path.splitext(f)[0]
	print(filename)
	PP0,Va_ccperg=np.genfromtxt(f,unpack=True,delimiter=';',skip_header=1,usecols=(10,15),encoding='iso-8859-1')

	isotherm=pygaps.PointIsotherm(material=filename,temperature=77,adsorbate='nitrogen',pressure=PP0,loading=Va_ccperg,pressure_mode='relative',pressure_unit='torr',loading_basis='molar',loading_unit='cm3(STP)',material_basis='mass',material_unit='g',temperature_unit='K')

	plt.close('all')
	isotherm.plot()
	plt.savefig(filename+'_pG_iso.png')
	plt.xscale('log'),plt.xlim([min(PP0)/2,2])
	plt.savefig(filename+'_pG_iso_log.png')

	plt.close('all')
	results.append((filename,pgc.area_BET(isotherm,verbose=True),pgc.area_langmuir(isotherm,verbose=True),pgc.t_plot(isotherm,thickness_model=carbon_model,verbose=True)))
	for i in plt.get_fignums()[0:]:
		plt.figure(i),plt.tight_layout(pad=0.1)
		plt.savefig(filename+'_pG_'+str(i-1)+'.png')

	plt.close('all')
	result_dict_meso=pgc.psd_mesoporous(isotherm);ax=psd_plot(result_dict_meso['pore_widths'],result_dict_meso['pore_distribution'],labeldiff='mesoporous',labelcum=None,method='comparison')
	result_dict_micro=pgc.psd_microporous(isotherm,psd_model='RY-CY');ax.plot(result_dict_micro['pore_widths'],result_dict_micro['pore_distribution'],label='microporous')
	# ~ result_dict_dft=pgc.psd_dft(isotherm);ax.plot(result_dict_dft['pore_widths'],result_dict_dft['pore_distribution'],label='DFT')
	ax.legend(loc='best')
	plt.savefig(filename+'_pG_PSD.png')

	plt.close('all')
	model_iso=pgm.model_iso(isotherm,model=['Henry','Langmuir','bet'],verbose=True)
	models.append(model_iso.model)
	plt.savefig(filename+'_pG_mod.png')

sys.stdout=open('results.log','w')
for i,(f,b,l,t) in enumerate(results):
	print(f+' BET: '+f"{b['area']:.2f}"+' Langmuir:'+f"{l['area']:.2f}"+' t:'+f"{t['results'][0]['area']:.2f}")
	print(models[i])
