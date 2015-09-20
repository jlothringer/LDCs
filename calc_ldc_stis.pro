pro calc_ldc_stis

;bins=[5500,6500,7500,8500,10000]

;bins = indgen(5,increment=(10100.-5315.)/4.,start=5315.)

bins = [5315.,10100.]

;This won't work anymore when we take out coeff from rdlimbg_55cnc
rdlimbg_55cnc,'lte5243_4.33_+0.25_55Cnc.irfout.27',lamb,flux,mu ;,coeff

ftab_ext,'p822207no_pht.fits','WAVELENGTH',wave
ftab_ext,'p822207no_pht.fits','THROUGHPUT',through
ftab_ext,'p822207no_pht.fits','OPT_ELEM',opt_elem

sorted = sort(lamb)
;stop
lamb = lamb(sorted)
flux=flux[*,sorted]
data=mrdfits('oco513010_raw_cleaned35x1d.fits',1)

good_through_waves = wave[where(wave[*,15] ne 0.0),15]
good_throughput = through[where(through[*,15] ne 0.0),15]

ind1 = where(lamb ge good_through_waves[0])
ind2 = where(lamb lt good_through_waves[-1])
ind = cgsetintersection(ind1,ind2)
good_waves = lamb(ind)
good_flux = flux(*,ind)

interpolated_throughput = interpol(good_throughput,good_through_waves,good_waves)

result = fltarr(n_elements(mu),n_elements(bins)-1)

for i = 0,n_elements(bins)-2 do begin
  bin_min = bins[i]
  bin_max = bins[i+1]
  w1 = where(good_waves gt bin_min)
  w2 = where(good_waves lt bin_max)
  w = cgsetintersection(w1,w2)
  ;for j = 0,n_elements(w)-1 do begin
    for k = 0,n_elements(mu)-1 do begin
      ;result(k,i) = total(reform(good_flux(k,w))*interpolated_throughput(w)*w)             ;Sum(I(mu,lamb)*S(lamb))
      result(k,i) = total(convol(reform(good_flux(k,w)),interpolated_throughput(w)))
    endfor
  ;endfor
endfor

Inorm = result

for k = 0,n_elements(mu)-1 do begin
  Inorm(k,*) = result(k,*) / result(-1,*)
endfor

;stop

;flux_through = dblarr(n_elements(mu),n_elements(good_waves))
;for i = 0,n_elements(mu)-1 do begin
;  flux_through(i,*) = good_flux(i,*) * interpolated_throughput * good_waves[*]   ;We want to use the 'observed' throughput instead
;endfor
;
;ind_mu = where(mu gt 0.1)
;good_mu = mu(ind_mu)
;good_flux_through = flux_through(ind_mu,*)

ldc = fltarr(2,n_elements(bins)-1)

for i = 0,n_elements(bins)-2 do begin
  a = [0.6,0.4]
  muzero = where(inorm(*,i) ge 1./exp(1))
  fit = curvefit(mu(muzero),inorm(muzero,i),weights,a,chi,/noderivative,function_name='limb_darkening')
  plot,mu(muzero),inorm(muzero,i),title='Range: ' + string(bins[i]) + ' to ' +string(bins[i+1])
  oplot,mu(muzero),fit,linestyle=5
  ldc[*,i]=a
  print, a
  ;stop
endfor
;bins = [5500,6500,7500,8500,10000]
;for i=0,n_elements(bins)-2 do begin
;  ;a=[0.6,0.4,0.1,0.1]
;  a=[0.6,0.4]
;  one2 = where(good_waves lt bins[i+1])
;  one1 = where(good_waves gt bins[i])
;  one = cgsetintersection(one1,one2)
;  one_sum = total(good_flux_through(*,one),2)
;  one_sum_norm = one_sum / one_sum(-1)
;  fit = curvefit(good_mu,one_sum_norm,weights,a,chi,/noderivative,function_name='limb_darkening')
;  plot,good_mu,one_sum_norm,title='Range: ' + string(bins[i]) + ' to ' +string(bins[i+1])
;  oplot,good_mu,fit,linestyle=5
;  ldc[*,i] = a
;  stop
;endfor
;Now do this for the other bin locations and fit coeff!
;Won't need ot use limb_params and fit in rdlimbg then?

stop
end