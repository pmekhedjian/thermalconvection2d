PRO loaddat,file

COMMON movdat,nx,nz,aspect,tem,psi

IF N_ELEMENTS(file) EQ 0 THEN RETURN

OPENR,1,file

nx=0L
nz=0L
aspect=0.
ra=0.
pr=0.
READF,1,nx,nz,aspect,ra,pr

tem=FLTARR(nx,nz)
psi=FLTARR(nx,nz)
omg=FLTARR(nx,nz)

READF,1,tem,psi

CLOSE,1
FREE_LUN,1

end
