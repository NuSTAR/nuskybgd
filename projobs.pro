; pulls mast motions from caldb files and/or loads image of said motions

pro projobs,indir,obsid,ab,posim,refminx,refminy,bgddir,$
      clobber=clobber,silent=silent

dir=indir
if strmid(dir,strlen(dir)-1) ne '/' then dir=dir+'/'
cldir=dir+obsid+'/event_cl/'
if size(bgddir,/type) eq 0 then dir=cldir else begin
    dir=cldir+bgddir+'/'
    if not file_test(dir,/directory) then spawn,'mkdir '+dir
endelse
if not keyword_set(clobber) then clobber=0
lclobber=clobber

if lclobber eq 2 then lclobber=0
if not file_test(dir+'projobs'+ab+'.fits') or lclobber ne 0 then begin
    if not file_test(dir+'det1'+ab+'gtifilter.dat') or clobber $
           then begin
        det1=mrdfits(cldir+'nu'+obsid+ab+'_det1.fits',1,hdet,/silent)
        gti=mrdfits(cldir+'nu'+obsid+ab+'01_gti.fits',1,hgti,/silent)
        undefine,x
        undefine,y
        undefine,dt
        for i=0,n_elements(gti.start)-1 do begin
            ii=where(det1.time ge gti[i].start and det1.time lt gti[i].stop)
            if ii[0] ne -1 then begin
                push,x,det1[ii].x_det1
                push,y,det1[ii].y_det1
                push,dt,det1[ii+1].time-det1[ii].time
            endif
        endfor
        openw,lun,dir+'det1'+ab+'gtifilter.dat',/get_lun
        printf,lun,'# ',min(floor(x)),min(floor(y))
        for j=0L,n_elements(x)-1 do printf,lun,x[j],y[j],dt[j]
        free_lun,lun
        if not keyword_set(silent) then $
              print,'GTI filtered wobble file made for '+ab
    endif else begin
        readcol,dir+'det1'+ab+'gtifilter.dat',x,y,dt,/silent
        if not keyword_set(silent) then $
               print,'GTI filtered wobble file found for '+ab
    endelse

    posim=fltarr(1000,1000)
    for i=min(floor(x)),max(floor(x)) do for j=min(floor(y)),max(floor(y)) do begin
        ii=where(floor(x) eq i and floor(y) eq j)
        if ii[0] ne -1 then posim[i,j]=total(dt[ii])
    endfor
    fits_write,dir+'projobs'+ab+'.fits',posim
    if not keyword_set(silent) then print,'Reference projection made for '+ab
endif else if not keyword_set(silent) then $
        print,'Reference projection found for '+ab

openr,lun,dir+'det1'+ab+'gtifilter.dat',/get_lun
line=''
readf,lun,line
free_lun,lun
vals=strsplit(line,/extract)
refminx=float(vals[1])
refminy=float(vals[2])
fits_read,dir+'projobs'+ab+'.fits',posim


end
