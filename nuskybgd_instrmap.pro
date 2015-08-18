; NUSKYBGD_INSTRMAP
;
; Creates an instrument map excluding RAW pixels contained in bad pixel files.
; The 'files' argument should be ignored; if you notice additional pixels that
; look dead, the correct action is to add them to a user bad pixel list and
; rerun the pipeline, in which case this routine will automatically grab the
; additional bad pixels.  If you're impatient and unwise, files gives the name 
; (and full path) to additional bad pixel lists -- which can be a single string or 
; array of strings.  But this is like drinking from the ocean to end your thirst.


pro nuskybgd_instrmap,indir,obsid,ab,bgddir,files,clobber=clobber

if not keyword_set(clobber) then clobber=0



;grade weighting from NGC253 002 obs.
gw=[1.00000,    0.124902,    0.117130,    0.114720,    0.118038,   0.0114296,$
    0.0101738,   0.0113617,   0.0122017,   0.0157910,   0.0144079,   0.0145691,$
    0.0149934,  0.00165462,  0.00194312,  0.00156128,  0.00143400,  0.00210433,$
   0.00180735,  0.00140006,  0.00169704,  0.00189220,  0.00160371,  0.00150188,$
   0.00168007, 0.000296983, 0.000364864]

caldb=getenv('CALDB')+'/'
auxildir=getenv('NUSKYBGD_AUXIL')+'/'

dir=indir
if strmid(dir,strlen(dir)-1) ne '/' then dir=dir+'/'
cldir=dir+obsid+'/event_cl/'

if size(bgddir,/type) eq 0 then outdir=cldir else begin
    outdir=cldir+bgddir+'/'
    if not file_test(outdir,/directory) then spawn,'mkdir '+outdir
endelse

dirinstrmap=caldb+'data/nustar/fpm/bcf/instrmap/'
dirpixpos=caldb+'data/nustar/fpm/bcf/pixpos/'
dirbadpix=caldb+'data/nustar/fpm/bcf/badpix/'

fits_read,dirinstrmap+'nu'+ab+'instrmap20100101v003.fits',instrmap,header
instrmap=shift(instrmap,-1,-1)

; check if images already exist
if not file_test(auxildir+'pixmap'+ab+'.fits') or clobber then begin
    file=dirpixpos+'nu'+ab+'pixpos20100101v005.fits'
    pixmap=fltarr(360,360)
    pixmap[*,*]=-1
    detnum=fltarr(360,360)
    detnum[*,*]=-1
    allpdf=fltarr(360,360)
    for idet=0,3 do begin
      pixpos=mrdfits(file,idet+1,hh,/silent)
      for ix=0,31 do for iy=0,31 do begin
        ii=where(pixpos.ref_det1x ne -1 and pixpos.rawx eq ix and $
              pixpos.rawy eq iy and pixpos.grade le 26)
        thispdf=fltarr(360,360)
        if ii[0] ne -1 then for i=0,n_elements(ii)-1 do $
              if finite(total(pixpos[ii[i]].pdf)) then $
              thispdf[pixpos[ii[i]].ref_det1x:pixpos[ii[i]].ref_det1x+6,$
                    pixpos[ii[i]].ref_det1y:pixpos[ii[i]].ref_det1y+6]+= $
                    pixpos[ii[i]].pdf*gw[pixpos[ii[i]].grade]
        ii=where(thispdf gt allpdf)
        if ii[0] ne -1 then begin
            allpdf[ii]=thispdf[ii]
            pixmap[ii]=ix+iy*32
            detnum[ii]=idet
        endif

; Older code chunk, ignores grade
;    x=pixpos[ii].ref_det1x-1
;    y=pixpos[ii].ref_det1y-1
;    pdf=pixpos[ii].pdf
;    grade=pixpos[ii].grade
;    rawx=pixpos[ii].rawx
;    rawy=pixpos[ii].rawy
;    for i=0,n_elements(x)-1 do if finite(total(pdf[*,*,i])) then begin
;        jj=where(pdf[*,*,i] gt 0.0)
;;        jj2d=array_indices(intarr(7,7),jj)
;        thispdf=intarr(7,7)
;        thispdf[*,*]=-1
;        thisdet=thispdf
;        thispdf[jj]=rawx[i]+rawy[i]*32
;        thisdet[jj]=idet
;        pixmap[x[i]:x[i]+6,y[i]:y[i]+6]=thispdf
;        detnum[x[i]:x[i]+6,y[i]:y[i]+6]=thisdet
;    endif

      endfor
    endfor
    pixmap=shift(pixmap,-1,-1)
    detnum=shift(detnum,-1,-1)
    fits_write,auxildir+'pixmap'+ab+'.fits',pixmap,header
    fits_write,auxildir+'detnum'+ab+'.fits',detnum,header
endif else begin
    fits_read,auxildir+'pixmap'+ab+'.fits',pixmap,header
    fits_read,auxildir+'detnum'+ab+'.fits',detnum,header
endelse

push,files,cldir+'nu'+obsid+ab+'_bp.fits'
;push,files,cldir+'nu'+obsid+ab+'01_cl.evt'
;push,files,dirbadpix+'nu'+ab+'badpix20100101v002.fits'
;file=auxildir+'nu'+ab+'userbadpix20100101v001.fits'

for ifiles=0,n_elements(files)-1 do begin
; uncomment below if cleaned event file used instead of *_bp.fits file
;  if ifiles eq n_elements(files)-1 then icl=2 else $
        icl=0
  for idet=0,3 do begin
    badpix=mrdfits(files[ifiles],idet+1+icl,hh,/silent)
    if n_elements(size(badpix)) eq 4 then begin
        x=badpix.rawx
        y=badpix.rawy
        for i=0,n_elements(x)-1 do begin
            ii=where(pixmap eq x[i]+y[i]*32 and detnum eq idet)
            instrmap[ii]=0
        endfor
    endif
  endfor
  ii=where(instrmap gt 0)
  instrmap[ii]+=detnum[ii]
endfor

fits_write,outdir+'newinstrmap'+ab+'.fits',instrmap,header


end
