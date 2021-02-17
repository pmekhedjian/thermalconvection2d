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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 
PRO mov,count,file,T_MIN=t_min0,T_MAX=t_max0

;  count is the first frame number in the movie sequence 

COMMON movdat,nx,nz,aspect,tem,psi

nxwin=480

nfile=STRLEN(file)

if N_ELEMENTS(t_min0) ne 0 then begin
   t_min=t_min0
endif else begin
   t_min=0.
endelse
if N_ELEMENTS(t_max0) ne 0 then begin
   t_max=t_max0
endif else begin
   t_max=1.
endelse

LOADCT,33
cr=bindgen(256)
cg=bindgen(256)
cb=bindgen(256)
TVLCT,cr,cg,cb,/GET

DEVICE,true_color=32

MPEGFILE=file+'.mpg'
parfile=file+'.par'
OPENW,2,parfile
PRINTF,2,"PATTERN ibpb"
PRINTF,2,"OUTPUT "+mpegfile
PRINTF,2,"INPUT_DIR ."
PRINTF,2,"INPUT"

while FILE_TEST(file) do begin

   loaddat,file
   nywin=round(nxwin/aspect)

;  in next line, change tem to psi if want to plot psi
;  and enter the max (t_max) and min (t_min) values
;  of psi when running this procedure.
   t=congrid(tem,nxwin,nywin,cubic=-0.5,/minus_one)
   t=bytscl(t,min=t_min,max=t_max)
   t=reverse(t,2)
   y=bytarr(3,nxwin,nywin)
   y(0,*,*)=cr(t(*,*))
   y(1,*,*)=cg(t(*,*))
   y(2,*,*)=cb(t(*,*))
   WRITE_PPM,file+'.ppm',y

   PRINT,'count=',count
   PRINTF,2,file+'.ppm'
   count=count+1
   file=STRMID(file,0,nfile-1-(count ge 10)-(count ge 100)-  $
      (count ge 1000))+STRTRIM(count,2)

endwhile

PRINTF,2,"END_INPUT"
PRINTF,2,"BASE_FILE_FORMAT PPM"
PRINTF,2,"GOP_SIZE 30"
PRINTF,2,"INPUT_CONVERT *"
PRINTF,2,"FORCE_ENCODE_LAST_FRAME"
PRINTF,2,"REFERENCE_FRAME DECODED"
PRINTF,2,"PIXEL FULL"
PRINTF,2,"SLICES_PER_FRAME 1"
PRINTF,2,"IQSCALE  6"
PRINTF,2,"PQSCALE 12"
PRINTF,2,"BQSCALE 16"
PRINTF,2,"RANGE 8"
PRINTF,2,"PSEARCH_ALG EXHAUSTIVE"
PRINTF,2,"BSEARCH_ALG CROSS2"
CLOSE,2
FREE_LUN,2

SPAWN,"mpeg_encode "+parfile

end
