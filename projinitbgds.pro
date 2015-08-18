; writes detector and bgd aperture images for use by nuskybgd routines

pro projinitbgds,indir,obsid,header,ab,bgddir,pa=pa,oldpa=oldpa,clobber=clobber,$
      grxe=grxe
                
auxildir=getenv('NUSKYBGD_AUXIL')+'/'

dir=indir
if strmid(dir,strlen(dir)-1) ne '/' then dir=dir+'/' 
cldir=dir+obsid+'/event_cl/'
if size(bgddir,/type) eq 0 then dir=cldir else begin
    dir=cldir+bgddir+'/'
    if not file_test(dir,/directory) then spawn,'mkdir '+dir
endelse

if not keyword_set(grxe) then grxe=0
if not keyword_set(clobber) then clobber=0

imexist=0
for det=0,3 do begin
      imexist+=file_test(dir+'bgdap'+str(det)+ab+'.fits')
      imexist+=file_test(dir+'det'+str(det)+ab+'im.fits')
endfor
if clobber ne 0 or imexist ne 8 then begin

;JRM: If pa not set, then read the PA from the log file:
IF keyword_set(pa) NE 1 THEN BEGIN
   if file_test(dir+'pa_gtifilter.dat') then begin
      readcol,dir+'pa_gtifilter.dat',dt,pa,/silent
      pa=total(pa*dt)/total(dt)
   endif else begin
; If you want to use this (or another) log file for the PA, set stillneedpa
;  to 0 and uncomment the following lines
      stillneedpa=1
;      if file_test(auxildir+'NuSTAR_Observing_Log.txt') then begin
;         readcol, auxildir+'NuSTAR_Observing_Log.txt', $
;               log_obsid, log_pa, $
;               format='A,X,X,X,X,X,X,X,X,X,F', $
;               comment=';', /silent
;         o = WHERE(strtrim(log_obsid,2) EQ strtrim(obsid,2), n_o)
;         if o[0] ne -1 then begin
;            pa = mean(log_pa[o])-90.
;            print,'Using PA from log file: '+str(pa)
;         endif else stillneedpa=1
;      endif
      if stillneedpa then begin
         evt=mrdfits(cldir+'nu'+obsid+'A01_cl.evt',1,head,/silent)
; old PA definition
;         pa = sxpar(head,'PA_PNT')-90.
         pa = sxpar(head,'PA_PNT')
         print,'Using PA from event file header keyword PA_PNT: '+str(pa)
      endif
   endelse
ENDIF

readcol,auxildir+'nomapbgdparams.dat',paperbgd,/silent

if not keyword_set(oldpa) then arot=-pa*!pi/180. else arot=(pa-90.)*!pi/180.

projobs,indir,obsid,ab,posim,refminx,refminy,bgddir,clobber=clobber
if grxe then grxe2det,indir,obsid,header,ab,bgddir,pa,clobber=clobber
blah=max(posim,ii)
ii2d=array_indices(posim,ii)
refminx=ii2d[0]
refminy=ii2d[1]
ii=where(posim gt 0)
ii2d=array_indices(posim,ii)
refimap0=fltarr(1000,1000)
refimap1=fltarr(1000,1000)
refimap2=fltarr(1000,1000)
refimap3=fltarr(1000,1000)
refimg0=fltarr(1000,1000)
refimg1=fltarr(1000,1000)
refimg2=fltarr(1000,1000)
refimg3=fltarr(1000,1000)
refimi0=fltarr(1000,1000)
refimi1=fltarr(1000,1000)
refimi2=fltarr(1000,1000)
refimi3=fltarr(1000,1000)
fits_read,dir+'newinstrmap'+ab+'.fits',instr
intdet=fltarr(4,360,360)
blah=fltarr(360,360)
for idet=0,3 do begin
    kk=where(instr eq idet+1)
    blah[kk]=1.
    intdet[idet,*,*]=blah
    blah[kk]=0.
endfor
fits_read,auxildir+'det'+ab+'_det1.img',apim
if grxe then fits_read,dir+'grxe'+ab+'_det1.img',grxeim
if ab eq 'A' then ifpm=0 else ifpm=1
apim=rot(fshift(apim,paperbgd[4+ifpm*3]/0.12096,paperbgd[5+ifpm*3]/0.12096),$
      paperbgd[6+ifpm*3],/interp)*3.1865e-5*paperbgd[ifpm]*0.12096*0.12096
apim=apim[108:467,108:467]
apim+=0.00495*0.7
if grxe then begin
    grxeim=rot(fshift(grxeim,paperbgd[4+ifpm*3]/0.12096,paperbgd[5+ifpm*3]/0.12096),$
          paperbgd[6+ifpm*3],/interp)*3.1865e-5*0.12096*0.12096 ; need to scale?????
    grxeim=grxeim[108:467,108:467]
endif
ia=intarr(1000,1000)
ja=intarr(1000,1000)
for i=0,999 do begin
    ia[*,i]=indgen(1000)
    ja[i,*]=indgen(1000)
endfor
if not keyword_set(oldpa) then begin
    arot=-pa*!pi/180.
    detxa=350-(ia-refminx)*cos(arot)+(ja-refminy)*sin(arot)
    detya=350+(ia-refminx)*sin(arot)+(ja-refminy)*cos(arot)
    nudge=1
endif else begin
    arot=(pa-90.)*!pi/180.
    detxa=350+(ia-refminx)*cos(arot)-(ja-refminy)*sin(arot)
    detya=350-(ia-refminx)*sin(arot)-(ja-refminy)*cos(arot)
    nudge=-1
endelse
jj=where(detxa ge 0 and detxa lt 360 and $
      detya ge 0 and detya lt 360)
blah=min(ia[jj],iblah)
istart=ia[jj[iblah]]
blah=max(ia[jj],iblah)
iend=ia[jj[iblah]]
blah=min(ja[jj],iblah)
jstart=ja[jj[iblah]]
blah=max(ja[jj],iblah)
jend=ja[jj[iblah]]
totposim=total(posim)
;    print,'Making reference image for '+ab
;    for i=0,999 do for j=0,999 do begin
for i=istart,iend do for j=jstart,jend do begin
;        detx=350+(i-refminx)*cos(arot)-(j-refminy)*sin(arot)
;        dety=350-(i-refminx)*sin(arot)-(j-refminy)*cos(arot)
    detx=detxa[i,j]
    dety=detya[i,j]
    if detx ge 0 and detx lt 360 and dety ge 0 and $
          dety lt 360 then if instr[detx,dety] gt 0 then begin
;        if detx ge 0 and detx lt 360 and dety ge 0 and dety lt 360 then $
;            if instr[detx,dety] gt 0 then begin
;                apval=aperbgd(ifpm,detx,dety,paperbgd)/total(posim)
        apval=apim[detx,dety]/totposim
        refimap0[i,j]=apval*intdet[0,detx,dety]
        refimap1[i,j]=apval*intdet[1,detx,dety]
        refimap2[i,j]=apval*intdet[2,detx,dety]
        refimap3[i,j]=apval*intdet[3,detx,dety]
        if grxe then begin
            grxeval=grxeim[detx,dety]/totposim
            refimg0[i,j]=grxeval*intdet[0,detx,dety]
            refimg1[i,j]=grxeval*intdet[1,detx,dety]
            refimg2[i,j]=grxeval*intdet[2,detx,dety]
            refimg3[i,j]=grxeval*intdet[3,detx,dety]
        endif
        refimi0[i,j]=intdet[0,detx,dety]/totposim
        refimi1[i,j]=intdet[1,detx,dety]/totposim
        refimi2[i,j]=intdet[2,detx,dety]/totposim
        refimi3[i,j]=intdet[3,detx,dety]/totposim
    endif
;        if ((i-istart) mod ((iend-istart)/10)) eq 0 and j eq jstart then $
;              print,(i-istart)*100/(iend-istart),'%...',format='($,I3,A)'
endfor
;    print
bgdimap0=fltarr(1000,1000)
bgdimap1=fltarr(1000,1000)
bgdimap2=fltarr(1000,1000)
bgdimap3=fltarr(1000,1000)
if grxe then begin
    bgdimg0=fltarr(1000,1000)
    bgdimg1=fltarr(1000,1000)
    bgdimg2=fltarr(1000,1000)
    bgdimg3=fltarr(1000,1000)
endif
bgdimi0=fltarr(1000,1000)
bgdimi1=fltarr(1000,1000)
bgdimi2=fltarr(1000,1000)
bgdimi3=fltarr(1000,1000)
for i=0,n_elements(ii)-1 do begin
    bgdimap0+=shift(refimap0,ii2d[0,i]-refminx-nudge,ii2d[1,i]-refminy-nudge)* $
            posim[ii2d[0,i],ii2d[1,i]]
    bgdimap1+=shift(refimap1,ii2d[0,i]-refminx-nudge,ii2d[1,i]-refminy-nudge)* $
            posim[ii2d[0,i],ii2d[1,i]]
    bgdimap2+=shift(refimap2,ii2d[0,i]-refminx-nudge,ii2d[1,i]-refminy-nudge)* $
            posim[ii2d[0,i],ii2d[1,i]]
    bgdimap3+=shift(refimap3,ii2d[0,i]-refminx-nudge,ii2d[1,i]-refminy-nudge)* $
            posim[ii2d[0,i],ii2d[1,i]]
    if grxe then begin
        bgdimg0+=shift(refimg0,ii2d[0,i]-refminx-nudge,ii2d[1,i]-refminy-nudge)* $
                posim[ii2d[0,i],ii2d[1,i]]
        bgdimg1+=shift(refimg1,ii2d[0,i]-refminx-nudge,ii2d[1,i]-refminy-nudge)* $
                posim[ii2d[0,i],ii2d[1,i]]
        bgdimg2+=shift(refimg2,ii2d[0,i]-refminx-nudge,ii2d[1,i]-refminy-nudge)* $
                posim[ii2d[0,i],ii2d[1,i]]
        bgdimg3+=shift(refimg3,ii2d[0,i]-refminx-nudge,ii2d[1,i]-refminy-nudge)* $
                posim[ii2d[0,i],ii2d[1,i]]
    endif
    bgdimi0+=shift(refimi0,ii2d[0,i]-refminx-nudge,ii2d[1,i]-refminy-nudge)* $
            posim[ii2d[0,i],ii2d[1,i]]
    bgdimi1+=shift(refimi1,ii2d[0,i]-refminx-nudge,ii2d[1,i]-refminy-nudge)* $
            posim[ii2d[0,i],ii2d[1,i]]
    bgdimi2+=shift(refimi2,ii2d[0,i]-refminx-nudge,ii2d[1,i]-refminy-nudge)* $
            posim[ii2d[0,i],ii2d[1,i]]
    bgdimi3+=shift(refimi3,ii2d[0,i]-refminx-nudge,ii2d[1,i]-refminy-nudge)* $
            posim[ii2d[0,i],ii2d[1,i]]
;        if (i mod n_elements(ii)/10) eq 0 then $
;              print,i*100/n_elements(ii),'%...',format='($,I3,A)'
endfor
fits_write,dir+'bgdap0'+ab+'.fits',bgdimap0,header
fits_write,dir+'bgdap1'+ab+'.fits',bgdimap1,header
fits_write,dir+'bgdap2'+ab+'.fits',bgdimap2,header
fits_write,dir+'bgdap3'+ab+'.fits',bgdimap3,header
if grxe then begin
    fits_write,dir+'bgdgrxe0'+ab+'.fits',bgdimg0,header
    fits_write,dir+'bgdgrxe1'+ab+'.fits',bgdimg1,header
    fits_write,dir+'bgdgrxe2'+ab+'.fits',bgdimg2,header
    fits_write,dir+'bgdgrxe3'+ab+'.fits',bgdimg3,header
endif
fits_write,dir+'det0'+ab+'im.fits',bgdimi0,header
fits_write,dir+'det1'+ab+'im.fits',bgdimi1,header
fits_write,dir+'det2'+ab+'im.fits',bgdimi2,header
fits_write,dir+'det3'+ab+'im.fits',bgdimi3,header
print,'Finished projection for '+ab

endif else begin

print,'Projinitbgds not run: projected bgd files already exist.'
print,'  To remake the images, set the clobber keyword.'

endelse

end
