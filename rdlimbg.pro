pro rdlimbg,name,lamb,flux,mu
;----------------------------
; read angular radiation field
; With full resolution
;
close,1 
;print,'reading file: ',name
on_ioerror,l666
;spawn,string("gunzip -c ",name,' >rdyng_tmp || rm rdyng_tmp'),ierr
;openr, 1,"rdyng_tmp"
openr, 1,name
i=0l
dummy=" "
;readf,1,dummy
;readf,1,dummy
readf,1,nmu
print,nmu

mu = dblarr(nmu)
lamb=dblarr(500000l)
flux=dblarr(nmu,500000l)
dflux = dblarr(nmu)

readf,1,mu

while not eof(1) do begin
readf,1,dlamb
readf,1,dflux
;readf,1,dummy
lamb(i) = dlamb
flux(*,i) = dflux
i=i+1l
endwhile
l666:
close,1
ipoint=i-1l
lamb=lamb(0l:i-1l)
flux=flux(*,0l:i-1l)
;spawn,string("rm rdyng_tmp"),ierr
end
