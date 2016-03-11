# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 11:50:59 2015

@author: jlothrin
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import rc
from pylab import *
from scipy import optimize
import emcee
import transit
import phasecurves as pc
#import time_series
import pdb
from matplotlib.backends.backend_pdf import PdfPages
import pickle
import tools
import triangle
#import transit_test
import analysis
from copy import deepcopy
from astropy.time import Time
import csv
#This is for WASP-31

def ldc_calc(star=None,filt=None,teff=None,logg=None,metal=None,binstart=500,binend=550,model='non-linear',compare=False,rescaled=True,plt=True):
    import numpy as np
    from scipy.interpolate import interp1d
    import ldtk
    from ldtk import (LDPSetCreator, BoxcarFilter,TabulatedFilter)
    from scipy.optimize import curve_fit
    #Import filters
    if filt == 'g430l':
        fil = np.loadtxt('/Users/jlothrin/Desktop/G430L.csv',delimiter=',',skiprows=1)
    if filt == 'g750l':
        fil = np.loadtxt('/Users/jlothrin/Desktop/G750L.csv',delimiter=',',skiprows=1)
    if filt == 'wfc3141':
        fil = np.loadtxt('/Users/jlothrin/Desktop/wfc3141.csv',delimiter=',',skiprows=1)
    if filt == 'IRAC_4_5':
        fil = np.loadtxt('/Users/jlothrin/Desktop/IRAC_4_5.csv',delimiter=',',skiprows=1)
    waves = fil[:,0]
    if filt[0] is 'g':
        waves = fil[:,0] / 10.
    if filt == 'IRAC_4_5':
        waves = fil[:,0]*1000.
    throughs = fil[:,1]
    f = interp1d(waves,throughs)
    waves_hires = np.linspace(min(waves),max(waves),500,endpoint=True)
    throughs_hires = f(waves_hires)
    
    
    #EVERYTHING MUST BE IN NM
    lower_wave = binstart
    upper_wave = binend
    
    w = [waves_hires[i] for i in range(0,len(waves_hires)) if waves_hires[i] > lower_wave if waves_hires[i] < upper_wave]
    t = [throughs_hires[i] for i in range(0,len(waves_hires)) if waves_hires[i] > lower_wave if waves_hires[i] < upper_wave]
        
    #filters = [BoxcarFilter('a',450,550), BoxcarFilter('b',650,750), BoxcarFilter('c',850,950)]
    #filters = [BoxcarFilter('a',650,750)]
    filters = [ldtk.TabulatedFilter('stis',w,t)]
    
    if star == '55cnc':
        sc = LDPSetCreator([5250, 100],[4.50, 0.10],[0.25, 0.05],filters)
        teff = 5250
        logg = 4.50
        metal = 0.25
    #sc = LDPSetCreator([6443,75],[4.76,0.09],[-0.08,0.05],filters)
    if star == 'wasp31':
        sc = LDPSetCreator([6250,50],[4.5,0.09],[-0.2,0.05],filters)
        teff = 6250
        logg = 4.5
        metal = -0.2
    if star == 'gj436':
        sc = LDPSetCreator([3416,100],[4.843,0.018],[0.02,0.2],filters)
        teff = 3416
        logg = 4.843
        metal = 0.02
    if star == 'gj3470':
        sc = LDPSetCreator([3652,50],[4.78,0.12],[0.18,0.08],filters)
        teff = 3652
        logg = 4.78
        metal = 0.18
    if star == 'hd97658':
        sc = LDPSetCreator([5217,33],[4.583,0.054],[-0.26,0.03],filters)
        teff = 5217
        logg = 4.583
        metal = -0.26
    if star == 'k218':
        sc = LDPSetCreator([3503,60],[4.863,0.13],[0.09,0.09],filters)
        teff = 3503
        logg = 4.863
        metal = 0.09
    if star == 'k138':
        sc = LDPSetCreator([3841,50],[4.886,0.01],[-0.280,0.1],filters)
        teff = 3841
        logg = 4.886
        metal = -0.280
    if star == 'wasp17':
        sc = LDPSetCreator([6666,30],[4.26,0.06],[-0.04,0.03],filters)
        teff = 6666
        logg = 4.26
        metal = -0.04
    if star is None:
        sc = LDPSetCreator([teff,50],[logg,0.1],[metal,0.1],filters)
                        
                      
    ps = sc.create_profiles(nsamples=100)
     
    qc,qe = ps.coeffs_qd(do_mc=False)
    nlc,nle = ps.coeffs_nl(do_mc=False)
    lnc,lne = ps.coeffs_ln(do_mc=False)
    print nlc
    print nle

    def nl_func(x,y1,y2,y3,y4):
        out = 1.-(y1*(1.-(x**(1./2.)))+y2*(1.-(x**(2./2.)))+y3*(1.-(x**(3./2.)))+y4*(1.-(x**(4./2.))))
        params = [y1,y2,y3,y4]
        bad = [1 for i in range(0,len(params)) if abs(params[i]) > 1.]
#        if np.sum(bad) > 0:
#            out = out+1000.
        return out
        
    fit= curve_fit(nl_func,ps._mu,ps.profile_averages[0],[0.2,0.9,-0.2,-0.03])
    
    if plt is True:
        figure()
        plot(ps._mu,ps.profile_averages[0],label='Model Intensity')
        plot(ps._mu,ps.profile_averages[0],'bo',label='Sampled Mu')
        #plot(ps._mu,1-qc[0,0]*(1-ps._mu)-qc[0,1]*(1-ps._mu)**2,'g--',label = 'Quad Fit')
        #plot(ps._mu,1-nlc[0,0]*(1-((ps._mu**(1./2.)))+nlc[0,1]*(1-ps._mu)+nlc[0,2]*(1-(ps._mu**(3./2.)))+nlc[0,3]*(1-(ps._mu**(2.)))),'r',label = 'Non-Linear Fit')
        plot(ps._mu,1-lnc[0,0]*(1-ps._mu),'y',label='Linear Fit')
        plot(ps._mu,nl_func(ps._mu,fit[0][0],fit[0][1],fit[0][2],fit[0][3]),'g',label = 'My non-linear fit')
        plot(ps._mu,nl_func(ps._mu,nlc[0][0],nlc[0][1],nlc[0][2],nlc[0][3]),'--k',label='Non-Linear Fit')

        
        #sing 450-460
        #plot(ps._mu,nl_func(ps._mu,0.1778,1.1859,-0.7235,0.1887),'k',label='Sing+2013')    
        #sing 460-470
        #plot(ps._mu,nl_func(ps._mu,0.2136,0.8688,-0.1927,-0.0317),'k',label='Sing+2013')
    
    
    new_mu = [(x-min(ps._mu))/(max(ps._mu)-min(ps._mu)) for x in ps._mu]
    
    refit= curve_fit(nl_func,new_mu,ps.profile_averages[0],[0.2,0.7,-0.2,-0.03])
    if plt is True:    
        plot(new_mu,ps.profile_averages[0],label='Rescaled Model Intensity')
        plot(new_mu,nl_func(np.asarray(new_mu),refit[0][0],refit[0][1],refit[0][2],refit[0][3]),'b',label = 'My rescaled non-linear fit')
    
    #my model for 55 cancri
    if compare is True:
        good_mu_old = [0.1313,0.19335,0.2554,0.31745,0.3795,0.44155,0.5036,0.56565,0.6277,0.68975,0.7518,0.8135,0.876,0.93795,1.0]
        one_sum_norm_old =      [ 0.37541731740693157,      0.46238495959398135,      0.52611395484286205,      0.57920999647415661,      0.62658004547351620,      0.67049257831603259,      0.71208564917444228,      0.75195835479518414,
              0.79045298237117534  ,    0.82776743852276846 ,     0.86401545614819741,      0.89927622441607569 ,     0.93360250587690652 ,     0.96703568174637822,       1.0000000000000000]
        good_mu = [0.1313,0.19335,0.2554,0.31745,0.3795,0.44155,0.5036,0.56565,0.6277,0.68975,0.7518,0.8135,0.876,0.93795,1.0]
        one_sum_norm = [0.38759735798258849, 0.46854567224818039,      0.52910322497992557,      0.58049544322307434,      0.62695066358655005,      0.67038264597838637,      0.71174173136349628,      0.75153440209321387,
              0.79003806657404141,      0.82741678251146478,      0.86376146640146356,      0.89913775107728544 ,     0.93359035990132300,      0.96715566457928104 ,      1.0000000000000000]
        
        if plt is True:
            plot(good_mu,one_sum_norm,'r--',label='Josh model')
            plot(good_mu,one_sum_norm,'ro',label='Josh Sampled Mu')
        
            plot(good_mu_old,one_sum_norm_old,'k--',label='Josh model old')
            plot(good_mu_old,one_sum_norm_old,'ko',label='Josh Sampled Mu old')
        
    #ylim([0.,1.0])
    #legend(loc='lower right')
    #ylabel('Normalized Intensity')
    #xlabel('Mu')
    #title('LDTK - Teff = 5250 logg = 4.33, z = 0.25 - 650-750 nm')
    #text(0.4,0.4,'Q1 = '+str(qc[0,0]))
    #text(0.4,0.36,'Q2 = '+str(qc[0,1]))
    if plt is True:
        legend(loc='lower right')
        ylabel('Normalized Intensity')
        xlabel('Mu')
        title('LDTK - Teff = '+str(teff)+', logg = '+str(logg)+', z = '+str(metal)+' - STIS - '+str(lower_wave)+'-'+str(upper_wave)+' nm')
    
    if rescaled is True:
        out = refit[0]
    if rescaled is False:
        out = fit[0]
       
    #pdb.set_trace()
        
    return out,nle
    
def ldc_calc_wrapper(star=None,filt=None,teff=None,logg=None,metal=None,model='non-linear',rescaled=True,plt=False,binstart=5315,binend=10100,nbins=11,csv_name=None,add_white=False,bins=None):
    import ldc_calc
    import numpy as np
    #bins for g750l?    
    if bins is None:
        bins = np.linspace(binstart,binend,nbins)/10.
    out=[]
    for i in range(0,len(bins)-1):
        out_tmp,nle=ldc_calc.ldc_calc(star=star,filt=filt,teff=teff,logg=logg,metal=metal,binstart=bins[i],binend=bins[i+1],model=model,compare=False,rescaled=rescaled,plt=plt)
        tmp_list = list(out_tmp)
        tmp_list.insert(0,bins[i+1])
        tmp_list.insert(0,bins[i])
        for i in range(0,4):
            tmp_list.append(nle[0][i])
        out.append(list(tmp_list))
    if csv_name is not None:
        b = open(csv_name,'w')
        a = csv.writer(b)
        a.writerows(out)
        b.close()
        out
    return out
    
def ldc_calc_k2(star=None,teff=None,tu=None,logg=None,lu=None,metal=None,metu=None,out_file=None,pickle_file=None,make_txt=True,make_fig=True):
    import numpy as np
    from scipy.interpolate import interp1d
    import ldtk
    from ldtk import (LDPSetCreator, BoxcarFilter,TabulatedFilter)
    from scipy.optimize import curve_fit
    kep_trans = np.loadtxt('/Users/jlothrin/LDCs/kepler_response_lowres1.csv',delimiter=',',skiprows=9)

    waves = kep_trans[:,0] * 1000.
    
    throughs = kep_trans[:,1]
    
    filters = [ldtk.TabulatedFilter('kepler',waves,throughs)]
    #filters = [BoxcarFilter('a',650,750)]
    sc = LDPSetCreator(teff=[teff,tu],logg=[logg,lu],z=[metal,metu],filters=filters)
    
    ps = sc.create_profiles(nsamples=100)
    
    #plot(ps._mu,ps.profile_averages[0],'k',label='Model Intensity')
    #plot(ps._mu,ps.profile_averages[0],'ko',label='Sampled Mu')
    oldmu = ps._mu
    oldprof = ps.profile_averages[0]
    
    ps.set_limb_z_try(sqrt(1-np.min(ps._mu)**2))
    
    qc,qe = ps.coeffs_qd(do_mc=True)
    qc,qm = ps.coeffs_qd(do_mc=True,return_cm=True)
    nlc,nle = ps.coeffs_nl(do_mc=True)
    lnc,lne = ps.coeffs_ln(do_mc=True)
    
    if make_txt is True:    
        text_file = open(out_file,'w')
        text_file.write(star + '\n \n')
        text_file.write('LINEAR: 1 - coeff*(1-mu) \n')
        text_file.write('Linear Coefficient: %s \n' % lnc[0,0])
        text_file.write('Linear Coefficient Uncertainty: %s \n \n' % lne[0])
        text_file.write('QUADRATIC: 1 - coeff1*(1-mu)-coeff2*(1-mu)**2 \n')
        text_file.write('Quadratic Coefficients: %s %f \n' % (qc[0,0],qc[0,1]))
        text_file.write('Quadratic Coefficients Uncertainty: %s %f \n' % (qe[0,0],qe[0,1]))
        text_file.close()
    
    #good_mu = ps_mu[2:]
    lw = 2.0
    ms = 3
    fig0 = figure()
    title(star)
    xlabel('mu')
    ylabel('I')
    ylim([0,1.0])
    plot(oldmu,oldprof,'r',label='Model Intensity',linewidth=lw)
    plot(oldmu,oldprof,'ro',label='Sampled Mu')
    plot(ps._mu,ps.profile_averages[0],'k',label='Rescaled Model Intensity',linewidth=lw)
    plot(ps._mu,ps.profile_averages[0],'ko',label='Rescaled Sampled Mu')
    plot(ps._mu,1-lnc[0,0]*(1-ps._mu),'y',label='Linear Fit',linewidth=lw)
    plot(ps._mu,1-qc[0,0]*(1-ps._mu)-qc[0,1]*(1-ps._mu)**2,'g--',label = 'Quad Fit',linewidth=lw)
    plot(ps._mu,ps.profile_averages[0]+ps.profile_uncertainties[0],'k:',label='Intensity Uncertainties')
    plot(ps._mu,ps.profile_averages[0]-ps.profile_uncertainties[0],'k:')

    
    #plot(new_mu,oldps,'r',label='Shifted')
    #plot(new_mu,oldps,'ro',label='Shifted')
    legend(loc='lower right')
    
    if make_fig is True:
        savefig(star+'ldc.png', format='png')    
        close(fig0)
    
    fig1 = figure()
    title(star)
    xlabel('z')
    ylabel('I')
    plot(ps._z,ps.profile_averages[0],label='Model Intensity')
    plot(ps._z,ps.profile_averages[0],'bo',label='Sampled z')
    plot(ps._z,1-lnc[0,0]*(1-ps._mu),'y',label='Linear Fit')
    plot(ps._z,1-qc[0,0]*(1-ps._mu)-qc[0,1]*(1-ps._mu)**2,'g--',label = 'Quad Fit')
    #plot(new_mu,oldps,'r',label='Shifted')
    legend(loc='lower right')
    close(fig1)

    variables = {'sc':sc,'ps':ps,'qc':qc,'qe':qe,'qm':qm,'lnc':lnc,'lne':lne}
    
    if pickle_file is not None:
        pickle.dump(variables,open(pickle_file,'wb'),protocol=2)
        
    return lnc, lne, qc, qe
    
def k2_wrapper(start):
    import ldc_calc
    stellar_params = np.loadtxt('/Users/jlothrin/LDCs/catalog_stellar_photparams_1800pc.csv',delimiter=',',skiprows=1)
    star = stellar_params[:,0]
    teff = stellar_params[:,1]
    uteff = stellar_params[:,2]
    logg = stellar_params[:,3]
    ulogg = stellar_params[:,4]
    metal = stellar_params[:,5]
    umetal = stellar_params[:,6]
    
    out_list = []
    for i in range(start,len(star)):
        lnc, lne, qc, qe = ldc_calc.ldc_calc_k2(star=str(star[i]),teff=teff[i],tu=uteff[i],logg=logg[i],lu=ulogg[i],metal=metal[i],metu=umetal[i],out_file=None,pickle_file=None,make_txt=False,make_fig=True)
        #out_list.append([star,lnc,lne,qc,qe])
        b = open('K2_giant_run.csv','a')
        a = csv.writer(b)
        a.writerows([[star[i],lnc[0][0],lne[0],qc[0][0],qc[0][1],qe[0][0],qe[0][1]]])
        b.close()
    
    return
            
        
        
    
