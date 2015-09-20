pro limb_params,mu,flux,coeff

;This procedure will loop through the fitting of limb darkening coeff for each wavelength
;We want to add up all of the intensity, weighted by the flux that we see in the data
;We will then normalize that and fit it

rdlimbg_55cnc,'lte5243_4.33_+0.25_55Cnc.irfout.27',lamb,flux,mu,coeff

n_mu = n_elements(mu)
size_flux = size(flux)
n_wavelengths = size_flux[2]
coeff = fltarr(2,n_wavelengths)
a=[0.657,0.115]
weights = fltarr(n_mu)

for i = 0,n_wavelengths-1 do begin
   ;weights = fltarr(n_mu)+1
   ;ratio = (flux[*,i])/(flux[n_mu,i])
   ;negatives = where(flux[*,i] le 0)
   ;weights[negatives] = 0
   ;weights(where(weights ne 0)) = 1
   tmp_flux = flux[*,i]
   ;ratio = flux[where(flux[*,i] gt 0)]/(flux[n_mu-1,i])
   ratio = tmp_flux(where(tmp_flux gt 0))/(tmp_flux(n_elements(tmp_flux)-1))
   mu_good = mu[where(tmp_flux gt 0)]
   weights = fltarr(n_elements(mu_good)) + 1   ;Should I weight for poisson error?
   fit = curvefit(mu_good,ratio,weights,a,chi,/noderivative,function_name='limb_darkening')
   coeff[0,i] = a[0]
   coeff[1,i] = a[1]
   ;stop
   ;values = 1.-a[0]*(1.-mu_good)-a[1]*((1.-mu_good)^2.)
   ;cgplot,mu_good,ratio
   ;cgoplot,mu_good,values
   stop
endfor

;stop
end