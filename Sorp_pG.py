"""
Created 21. November 2022 by Daniel Van Opdenbosch, Technical University of Munich

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. It is distributed without any warranty or implied warranty of merchantability or fitness for a particular purpose. See the GNU general public license for more details: <http://www.gnu.org/licenses/>
"""

import os
import sys
import glob
import numpy
import pygaps
import matplotlib as mpl
import matplotlib.pyplot as plt
import pygaps.characterisation as pgc
from pygaps.graphing.calc_graphs import psd_plot

def carbon_model(relative_p):
	return 0.88*(relative_p**2)+6.45*relative_p+2.98

files=glob.glob('*.csv')
results=[]

for f in files:
	filename=os.path.splitext(f)[0]
	print(filename)
	PP0,Va_ccperg=numpy.genfromtxt(f,unpack=True,delimiter=';',skip_header=1,usecols=(10,15),encoding='iso-8859-1')

	isotherm=pygaps.core.pointisotherm.PointIsotherm(material=filename,temperature=77,adsorbate='nitrogen',pressure=PP0,loading=Va_ccperg,pressure_mode='relative',pressure_unit='torr',loading_basis='molar',loading_unit='cm3(STP)',material_basis='mass',material_unit='g',temperature_unit='K')

	plt.clf()
	isotherm.plot()
	# ~ plt.xscale('log'),plt.xlim([min(PP0)/2,2])
	plt.savefig(filename+'_pG_iso.png')
	plt.close('all')

	plt.clf()
	results.append((filename+'_BET',pgc.area_BET(isotherm,verbose=True)))
	for i in plt.get_fignums()[1:]:
		plt.figure(i),plt.tight_layout(pad=0.1)
		plt.savefig(filename+'_pG_BET_'+str(i-1)+'.png')
	plt.close('all')

	plt.clf()
	results.append((filename+'_t-plot',pgc.t_plot(isotherm,thickness_model=carbon_model,verbose=True)['results'][0]))
	plt.figure(2),plt.tight_layout(pad=0.1)
	plt.savefig(filename+'_pG_tplot.png')
	plt.close('all')

	plt.clf()
	result_dict_meso=pgc.psd_mesoporous(isotherm);ax=psd_plot(result_dict_meso['pore_widths'],result_dict_meso['pore_distribution'],labeldiff='mesoporous',labelcum=None,method='comparison')
	result_dict_micro=pgc.psd_microporous(isotherm);ax.plot(result_dict_micro['pore_widths'],result_dict_micro['pore_distribution'],label='microporous')
	# ~ result_dict_dft=pgc.psd_dft(isotherm);ax.plot(result_dict_dft['pore_widths'],result_dict_dft['pore_distribution'],label='DFT')
	ax.legend(loc='best')
	plt.savefig(filename+'_pG_PSD.png')
	plt.close('all')

sys.stdout=open('results.log','w')
print([(x,f"{y['area']:.2f}") for (x,y) in results])
