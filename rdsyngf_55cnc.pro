pro rdsyngf_55cnc,name,lamb,flux,bb,flux2
;---------------------------
; read synthetic spectra 
; With full resolution
;
lamb=dblarr(5000000l)
flux=lamb
flux2=lamb
bb=lamb
close,1 
;print,'reading file: ',name

;-- if in gzip format
if (strpos(name,'gz') ge 0) then begin
 spawn,string("gunzip -c "+name+' | egrep -v -i nan > ./rdsyn_7.tmp'),ierr
 openr, 1,"./rdsyn_7.tmp"
endif

if (strpos(name,'bz2') ge 0) then begin
 spawn,string("bunzip2 -c "+name+' | egrep -v -i nan > ./rdsyn_7.tmp'),ierr
 openr, 1,"./rdsyn_7.tmp"
endif

i=0l
dlamb = 0.d0
dflux = 0.d0
while not eof(1) do begin
on_ioerror, bad_rec
;readf,1,dlamb,dflux,dbb,df2,format='(D23.15,D22.15,D22.15,D22.15)'
readf,1,dlamb,dflux,dbb,df2
lamb(i) = dlamb
flux(i) = dflux
flux2(i) = df2
bb(i) = dbb
i=i+1l
endwhile
close,1
if (i eq 0) then begin
 bad_rec:
 lamb = 0
 flux = 0
 return
endif
ipoint=i-1l
lamb=lamb(0l:i-1l)
;flux=10.d0^(flux(0l:i-1l)-2.0*alog10(radius)+alog10(4.*!pi))
flux=10.d0^(flux(0l:i-1l))
flux2=10.d0^(flux2(0l:i-1l)-40.0)
bb=10.d0^(bb(0l:i-1l))
il = sort(lamb)
lamb=lamb(il)
flux=flux(il)
flux2=flux2(il)
bb=bb(il)
;spawn,string("rm ./rdsyn_7.tmp"),ierr
stop
end
