#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 14 16:51:02 2021

@author: aboubakari
"""
#------------------DECLARATION DES LIBRAIRIES---------------------------

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import sys
sys.path.insert(0, '/home/scola/routines/')
import routines_AFM as rt

Fscale = 1e-6
zscale = 1e-9
PLOTFIGLOG = False
SAVEFIGLOG = False
PLOTFIGLIN = True
SAVEFIGLIN = False

DEBUG = False

#stiffness value of cantilever
k = 180
R = 17e-9

#-------------------------MAIN------------------------------------------
    
rep = '/home/scola/RECHERCHE/ZnO/Nanovibe/STAGE_Aboubakari/experiences/2021-07-07_to_22_AFM-ramps/2021-07-07/extractions/'
rt.MakeRampTable(rep, fit_crit = 1e-2, ifit_offset = 4)
samples = rt.ImportRampTable(rep)

samples_list = []
for rampname in samples.keys():
    if 'MV969_roi07-2_0' in rampname and 'retract' in rampname and int (rampname[19:21]) in range (2,17):
        samples_list.append (rampname)
samples_list = ['MV969_roi07-2_0-00016_retract', 'MV969_roi07-2_0-00013_retract', 'MV969_roi07-2_0-00004_retract']
samples_list = ['MV969_roi07-2_0-00016_retract']

if DEBUG: print('Files to process are : \n', samples_list)


DATA = {}
for sample_to_plot in samples_list:
    #------------------IMPORTATION---------------------
    samp=samples[sample_to_plot]

    file = samp['file']
    rep = samp['rep']
    f = open(rep+file, 'r', encoding='latin1')
    content = f.read()
    f.close()
    content = content.replace(',', '.') 
    lines = content.split('\n')
    
    header_len = 1
    header = ''
    for i in range(header_len):
        header += lines[i] + '\n'
    NBlines = len(lines)
    z = []
    ramp = []
    
    for i in range(NBlines-header_len-1):
        zz, ra = lines[i+header_len].split('\t')
        z.append(float(zz))
        ramp.append(float(ra))
    
    NB_points = len(z)
    
    # reverse the retract arrays stored in chronological order
    # to start at z = 0
    if 'retract' in sample_to_plot:
        zr = []
        rr = []
        for i in range(NB_points):
            zr.append(z[NB_points-i-1])
            rr.append(ramp[NB_points-i-1])
        z = zr
        ramp = rr
    
    z = np.asarray(z) * zscale
    ramp = np.asarray(ramp) * Fscale
    
    
    ifit_max=int (samp['ifit_max'] * NB_points) - 1
    fit_crit=samp['fit_crit']
    ifit_offset=samp['ifit_offset']
    
    #-------------------------MAIN------------------------------------------
    
    # Correction1 = centrer à 0 (baseline)
    
    imin = 0
    imax = int(0.61 * NB_points)
    slope, intercept = np.polyfit(z[imin:imax], ramp[imin:imax], 1)
    
    ramp1 = ramp - (intercept + slope * z)
    
    
    
    l = len(ramp1)
    m = max(ramp1)
    rmin = min(ramp1)
    # icontact : zone de contact
    contact_crit = 1e-3
    i = 0
    while ramp1[l-i-1] > contact_crit*m:
        i += 1
    i0=l-i
    icontact = l - i - 250
    zmin = z[icontact]
    
    # ifit : zone d'ajustement des données
    
    
    i = 0
    while ramp1[l-i-1] > fit_crit*m:
        i += 1
    ifit = l - i + ifit_offset
    
    d0 = z[i0]
    d = z[i0:ifit_max] - d0
    h = (z[i0:ifit_max] - d0) - ramp1[i0:ifit_max]/k
    hfit = (z[ifit:ifit_max] - d0) - ramp1[ifit:ifit_max]/k
    
    #fit Sneddon sphere
    popt, pcov = curve_fit(rt.sneddon_sphere, hfit, ramp1[ifit:ifit_max] - rmin)
    perr = 2*np.sqrt(np.diag(pcov))
    alpha = popt[0]
    
    if PLOTFIGLIN:
        fig1 = plt.figure(figsize=(10,14))
        ax1 = fig1.add_subplot(211)
        
        ax1.plot((z[icontact:]-d0 - ramp1[icontact:]/k)/zscale, (ramp1[icontact:] - rmin)/Fscale, c='k', label='correction raideur')
        ax1.plot((z[icontact:]-d0)/zscale, (ramp1[icontact:] - rmin)/Fscale, c='grey', label='exp. data')
        ax1.plot(hfit/zscale, rt.sneddon_sphere(hfit, alpha)/Fscale, ls=':', c='r', label='Sneddon model')
        ax1.set_title(file)
        ax1.set_xlabel(f'd ({zscale:.0e} m)')
        ax1.set_ylabel(f'Force ({Fscale:.0e} N)')
        ax1.set_xlim((zmin-d0)/zscale, (max(z)-d0)/zscale)
        fitparamS = '$F_\mathrm{Sneddon}$ = '+ f'{alpha:.2e}'+f' $h^{3/2}$\nRelative fit error : {perr[0]/popt[0]:.2e}\n'
    #    fitparamH = '$F_\mathrm{Hooke}$ = '+ f'{a:.2e}'+' $ d + $ '+f'{b:.2e}'
    #    ax1.text(0.05, 0.7, fitparamH, ha = 'left', va = 'center', transform = ax1.transAxes)
        ax1.text(0.05, 0.25, fitparamS, ha = 'left', va = 'center', transform = ax1.transAxes)
        ax1.legend()    
        
        ax2 = fig1.add_subplot(212)
        ax2.plot((z[icontact:]-d0 -ramp1[icontact:]/k)/zscale, (ramp1[icontact:] - rmin)/Fscale, c='k', label='correction raideur')
        ax2.plot(hfit/zscale, rt.sneddon_sphere(hfit, alpha)/Fscale, ls=':', c='r', label='Sneddon model')
        ax2.set_xlabel(f'h ({zscale:.0e} m)')
        ax2.set_ylabel(f'Force ({Fscale:.0e} N)')
        fitparamS = '$F_\mathrm{Sneddon}$ = '+ f'{alpha:.2e}'
    #    fitparamH = '$F_\mathrm{Hooke}$ = '+ f'{a:.2e}'+' $ d + $ '+f'{b:.2e}'
    #    ax1.text(0.05, 0.7, fitparamH, ha = 'left', va = 'center', transform = ax1.transAxes)
        ax2.text(0.05, 0.25, fitparamS, ha = 'left', va = 'center', transform = ax2.transAxes)
        ax2.legend()

        fig4 = plt.figure()
        ax4 = fig4.add_subplot(111)
        ax4.plot((z[icontact:]-d0 - ramp1[icontact:]/k)/zscale, (ramp1[icontact:] - rmin)/Fscale, c='k', label='Tip indentation')
        ax4.plot((z[icontact:]-d0)/zscale, (ramp1[icontact:] - rmin)/Fscale, c='grey', label='AFM signal (tip indentation + cantelever flexion')
        ax4.plot(hfit/zscale, rt.sneddon_sphere(hfit, alpha)/Fscale, ls=':', c='r', label='Indentation Hertz model')
        ax4.set_xlabel(f'd (nm)')
        ax4.set_ylabel(f'Force ($\mu$N)')
        ax4.set_xlim(-10e-9/zscale, (max(z)-d0)/zscale)
        ax4.set_ylim(-0.1e-6/Fscale, 3.8e-6/Fscale)
        ax4_label = r'$F = \frac{4}{3} \sqrt{R} \frac{E_r}{1-\nu^2} \times h^{3/2}$'+'\n'
        ax4_label += r'$\frac{1}{E_r} = \frac{(1-\nu^2)}{E_{ZnO}} + \frac{(1-\nu_\mathrm{tip}^2)}{E_\mathrm{tip}}$'+'\n\n'
        ax4_label += 'Experimental tip radius\n$R = 17$ nm\nExtracted Young\'s modulus\n$E_{ZnO}$ = 129 GPa'
        ax4.text(0.02, 0.6, ax4_label, ha = 'left', va = 'center', transform = ax4.transAxes)
        ax4.legend()    
                
    
    if SAVEFIGLIN: plt.savefig(rep+'figures/'+sample_to_plot+'_lin.png')
    
    if PLOTFIGLOG:
        fig2 = plt.figure(figsize = (5,4))
        ax3 = fig2.add_subplot(111)
        
        ax3.plot((z[icontact:]-d0 -ramp1[icontact:]/k)/zscale, (ramp1[icontact:] - rmin)/Fscale, c='k', label='correction raideur')
        hlog = (z[i0:NB_points-1] - d0) - ramp1[i0:NB_points-1]/k
        ax3.loglog(hlog/zscale, rt.sneddon_sphere(hlog, alpha)/Fscale, ls=':', c='r', label='Sneddon model')
        ax3.set_xlabel(f'h ({zscale:.0e} m)')
        ax3.set_ylabel(f'Force ({Fscale:.0e} N)')
        ax3.set_title(file)
        ax3.set_xlim(1e-11/zscale, (max(z)-d0)/zscale)
        ax3.set_ylim(1e-7/Fscale, 2e-5/Fscale)
        fitparamS = '$F_\mathrm{Sneddon}$ = '+ f'{alpha:.2e}'+' $h^{3/2}$\n'+f'Relative fit error : {perr[0]/popt[0]:.2e}\n'
    #    fitparamH = '$F_\mathrm{Hooke}$ = '+ f'{a:.2e}'+' $ d + $ '+f'{b:.2e}'
    #    ax1.text(0.05, 0.7, fitparamH, ha = 'left', va = 'center', transform = ax1.transAxes)
        ax3.text(0.05, 0.25, fitparamS, ha = 'left', va = 'center', transform = ax3.transAxes)
        ax3.legend()    
    

    if SAVEFIGLOG: plt.savefig(rep+'figures/'+file[0:-4]+'_log.png')
    
    data = {}
    nu = rt.nu_ZnO
    data['NW number'] = int (sample_to_plot[19:21]) - 1
    data['alpha'] = alpha
    Er = 3 * alpha / 4 / np.sqrt(R)
    data['Er'] = Er
    E = (1 - nu*nu) / ( 1/Er - (1 - rt.nu_Si*rt.nu_Si) / rt.E_Si)
    data['E'] = E
    data['E_th'] = rt.E_ZnO
    data['Fit error'] = perr[0]/popt[0]
    data['Elasticity ratio'] = np.abs (E/rt.E_ZnO)
    data['Fit error'] = perr[0] / popt[0]
    DATA[sample_to_plot] = data
#------------------------------OUTPUT------------------------------------------

out  = '|NW number  | $\\alpha$ | Fit error | $E$ | $E_\mathrm{theo}^{ZnO}$ | $E_\mathrm{theo}^{SOG}$ | $E/E_\mathrm{theo}^{ZnO}$ |\n'
out += '| :- | :-: | :-: | :-: | :-: | :-: | :-: |'

for meas in DATA:
    nNW = DATA[meas]['NW number']
    alpha = DATA[meas]['alpha']
    Er = DATA[meas]['Er']
    E = DATA[meas]['E']
    E_th = DATA[meas]['E_th']
    efrac = DATA[meas]['Elasticity ratio']
    err = DATA[meas]['Fit error']
    out += f'\n| {nNW:d} ' + f'| {alpha:.2e} | {err*100:.1f} \% | {E:.2e} | {E_th:.2e} | {rt.E_SOG:.2e} | {efrac*100:.1f} \% |'
print(out)  