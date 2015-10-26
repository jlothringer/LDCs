# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 12:50:49 2015

@author: jlothrin
"""

from astropy import constants as const
from astropy import units as u
import ldtk
from ldtk import (LDPSetCreator, BoxcarFilter,TabulatedFilter)

#LDC's for K2-9
mass = 0.30
umass = 0.14
radius = 0.31
uradius = 0.11

logg = np.log10((const.G.si * (mass) * const.M_sun.si / (((radius) * const.R_sun.si)**2)).cgs.value)

ug = np.sqrt(((const.G.si / const.R_sun.si**2)*umass*const.M_sun.si)**2+((-2*const.G.si*const.M_sun.si / const.R_sun.si**3)*uradius*const.R_sun.si)**2)

#g = [(const.G.si * (mass + j) * const.M_sun.si / (((radius +i) * const.R_sun.si)**2)).cgs.value for i in np.linspace(-uradius,uradius,25) for j in np.linspace(-umass,umass,25)]

g = [(const.G.si * (mass + j) * const.M_sun.si / (((radius +i) * const.R_sun.si)**2)).cgs.value for i in [-uradius,uradius] for j in [-umass,umass]]

#Now for LDCs
kep_trans = np.loadtxt('/Users/jlothrin/LDCs/kepler_response_lowres1.txt',delimiter='	',skiprows=9)

waves = kep_trans[:,0] * 1000.

throughs = kep_trans[:,1]

filters = [ldtk.TabulatedFilter('kepler',waves,throughs)]

sc = LDPSetCreator(teff=[3390,150],logg=[4.9,0.55],z=[-0.25,0.20],filters=filters)

ps = sc.create_profiles(nsamples=500)

plot(ps._mu,ps.profile_averages[0],'k',label='Model Intensity')
plot(ps._mu,ps.profile_averages[0],'ko',label='Sampled Mu')
oldmu = ps._mu
oldprof = ps.profile_averages[0]

ps.set_limb_z_try(sqrt(1-np.min(ps._mu)**2))

qc,qe = ps.coeffs_qd(do_mc=True)
qc,qm = ps.coeffs_qd(do_mc=True,return_cm=True)
nlc,nle = ps.coeffs_nl(do_mc=True)
lnc,lne = ps.coeffs_ln(do_mc=True)

#good_mu = ps_mu[2:]

figure()
title('K2-9')
xlabel('mu')
ylabel('I')
plot(oldmu,oldprof,'r',label='Model Intensity')
plot(oldmu,oldprof,'ro',label='Sampled Mu')
plot(ps._mu,ps.profile_averages[0],'k',label='Rescaled Model Intensity')
plot(ps._mu,ps.profile_averages[0],'ko',label='Rescaled Sampled Mu')
plot(ps._mu,1-lnc[0,0]*(1-ps._mu),'y',label='Linear Fit')
plot(ps._mu,1-qc[0,0]*(1-ps._mu)-qc[0,1]*(1-ps._mu)**2,'g--',label = 'Quad Fit')
#plot(new_mu,oldps,'r',label='Shifted')
#plot(new_mu,oldps,'ro',label='Shifted')
legend(loc='lower right')

figure()
title('K2-9')
xlabel('z')
ylabel('I')
plot(ps._z,ps.profile_averages[0],label='Model Intensity')
plot(ps._z,ps.profile_averages[0],'bo',label='Sampled z')
plot(ps._z,1-lnc[0,0]*(1-ps._mu),'y',label='Linear Fit')
plot(ps._z,1-qc[0,0]*(1-ps._mu)-qc[0,1]*(1-ps._mu)**2,'g--',label = 'Quad Fit')
#plot(new_mu,oldps,'r',label='Shifted')
legend(loc='lower right')


