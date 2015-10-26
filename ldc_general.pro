pro ldc_general

  ;bins=[5500,6500,7500,8500,10000]
  bins = [6500,7500]

;Grab wavelength and flux information from the Phoenix model
  rdlimbg,'436_lowactivity_fort.27',lamb,flux,mu

;Grab throughput information from STIS reference files
  ftab_ext,'p822207no_pht.fits','WAVELENGTH',wave
  ftab_ext,'p822207no_pht.fits','THROUGHPUT',through
  ftab_ext,'p822207no_pht.fits','OPT_ELEM',opt_elem

  sorted = sort(lamb)
  ;stop
  lamb = lamb(sorted)
  flux=flux[*,sorted]
  ;Import spectral data from 55 Cancri STIS observations (could be used instead of throughput info?)
  data=mrdfits('oco513010_raw_cleaned35x1d.fits',1)

  good_through_waves = wave[where(wave[*,15] ne 0.0),15]
  good_throughput = through[where(through[*,15] ne 0.0),15]

  ind1 = where(lamb ge good_through_waves[0])
  ind2 = where(lamb lt good_through_waves[-1])
  ind = cgsetintersection(ind1,ind2)
  good_waves = lamb(ind)
  good_flux = flux(*,ind)

  interpolated_throughput = interpol(good_throughput,good_through_waves,good_waves)
  ;interpolated_throughput = fltarr(n_elements(good_waves))+1

  flux_through = dblarr(n_elements(mu),n_elements(good_waves))
  for i = 0,n_elements(mu)-1 do begin
    flux_through(i,*) = good_flux(i,*) * interpolated_throughput; * good_waves[*]   ;We may want to use the 'observed' throughput instead
  endfor

  ind_mu = where(mu gt 0.07979)
  good_mu = mu(ind_mu)
  good_flux_through = flux_through(ind_mu,*)

  ;stop

  ldc = fltarr(2,n_elements(bins)-1)
  ;bins = [5500,6500,7500,8500,10000]
  ;cgwindow,ysize=800,xsize=1200
  for i=0,n_elements(bins)-2 do begin
    ;a=[0.6,0.4,0.1,0.1]
    a=[0.6,0.4]
    one2 = where(good_waves lt bins[i+1])
    one1 = where(good_waves gt bins[i])
    one = cgsetintersection(one1,one2)
    one_sum = total(good_flux_through(*,one),2)
    one_sum_norm = one_sum / one_sum(-1)
    fit = curvefit(good_mu,one_sum_norm,weights,a,chi,/noderivative,function_name='limb_darkening')
    cgplot,good_mu,one_sum_norm,title = 'Josh Model - Teff = 5243 logg = 4.33 z=+0.25 - 650-750 nm',xtitle = 'Mu' , ytitle = 'Normalized Intensity',thick=3,/window,color='blue'
        ;title='Range: ' + string(bins[i]) + ' to ' +string(bins[i+1])
    cgoplot,good_mu,fit,linestyle=5,thick=3,/window
    cgoplot,good_mu,one_sum_norm,psym=4,thick=3,/window,color='blue'
    ldc[*,i] = a
    cgtext,0.4,0.4,'Q1 = '+string(a[0,i]),/window
    cgtext,0.4,0.36,'Q1 = '+string(a[1,i]),/window
    ;cgLegend, alignment=3,colors=['blue','black'],linestyles=[0,5],titles=['PHX Model','Quad Fit'],/window;,symcolors='blue',psyms=4,/window
  endfor
  ;Now do this for the other bin locations and fit coeff!
  ;Won't need ot use limb_params and fit in rdlimbg then?
  
  print, ldc
  stop
end