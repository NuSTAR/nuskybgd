; NUSKYBGD_FITAB
;
; Simultaneously fits an empirical model to multiple bgd spectra and writes
; an output parameter file that is used by nuskybgd_spec.pro and nuskybgd_image.pro.
;
; This code is adapted from nuskybgd_fit, allowing A and B spectra from similar
; regions to be fit simulataneously so that the fCXB normalizations can be tied
; together, reducing some degeneracy at lower energies.
;
; See nuskybgd.pro for typical inputs.
;
; OPTIONAL KEYWORDS
;
;   pa --> user set PA: Only used if reference background images not yet created,
;          or clobber keyword set.  If your PA follows the old convention, you
;          need to run projinitbgds.pro with /oldpa set instead.
;
;   paramfile --> name and full path for the output file:  If not set, the
;                 parameter file is written to 
;                     indir+'/'+obsid+'/event_cl/'+specdir+'/bgdfitparams[A/B].dat
;
;   nocheck --> Set this if you trust that the fit will be good.  Otherwise
;               the routine pauses inside XSPEC so you can evaluate and tweak
;               the fit.
;   paramstem --> stem of the input parameter file. Assumes
;               that you've run nuskybgd_fitab first.
;
;   fixfcxb --> Set to fix the focused CXB normalization to its nominal value.
;               WARNING: The CXB variance could be 100% depending on the size 
;                        of the regions, so using the nominal value is dangerous.
;               NOTE: Despite the warning above, the fCXB is linked between all
;                     the spectra, which is also likely to be incorrect.  When in 
;                     XSPEC, you can untie these if you like.



pro nuskybgd_setparams,indir,obsid,bgdreg,specdir,specname,ab,bgddir,header,$
                       clobber=clobber,pa=pa,paramfile=paramfile, $
                       grxe=grxe,iisrc=iisrc, paramstem=paramstem

absave=ab
auxildir=getenv('NUSKYBGD_AUXIL')+'/'
caldbdir=getenv('CALDB')+'/'
dir=indir
if strmid(dir,strlen(dir)-1) ne '/' then dir=dir+'/'
cldir=dir+obsid+'/event_cl/'
clspecdir=dir+obsid+'/event_cl/'+specdir+'/'
if size(bgddir,/type) eq 0 then dir=cldir else begin
    dir=cldir+bgddir+'/'
    if not file_test(dir,/directory) then spawn,'mkdir '+dir
endelse
if not keyword_set(grxe) then grxe=0
if not keyword_set(iisrc) then iisrc=replicate(0,n_elements(bgdreg))
if ~keyword_set(paramstem) then paramstem='bgdparams'


if ab eq 'AB' or ab eq 'ab' or ab eq 'Ab' or ab eq 'aB' then begin
    abarr=[replicate('A',n_elements(bgdreg)/2),replicate('B',n_elements(bgdreg)/2)]
    print,'WARNING: NUSKYBGD_FITAB is assuming you have supplied corresponding'
    print,'   !     spectra for A and B in corresponding order, e.g.,'
    print,'   !       bgdreg=["bgd1A.reg","bgd2A.reg","bgd1B.reg","bgd2B.reg"]'
    print,'   !     where bgd1A.reg and bgd1B.reg, etc., cover nearly the same'
    print,'   !     area on the sky and all the A spectra/regions are listed first.'
    if n_elements(bgdreg) mod 2 ne 0 or n_elements(specname) mod 2 ne 0 then $
          stop,'NUSKYBGD_FITAB: Silly person, you muddled your input -- try again!'
endif else abarr=replicate(ab,n_elements(bgdreg))



pt=loadnuabs(0)
czt=loadnuabs(1)

bgdbackscl=fltarr(n_elements(bgdreg))
bgddetfrac=fltarr(n_elements(bgdreg),4)
bgdapfrac=fltarr(n_elements(bgdreg),4)
aptot=fltarr(n_elements(bgdreg),4)
if grxe then bgdgrxefrac=fltarr(n_elements(bgdreg),4)
if grxe then grxetot=fltarr(n_elements(bgdreg),4)
dettot=fltarr(n_elements(bgdreg),4)
detfrac=fltarr(n_elements(bgdreg),4)
bgdpt=fltarr(n_elements(bgdreg))
bgdczt=fltarr(n_elements(bgdreg))
instrnorms=fltarr(n_elements(bgdreg),n_elements(aeline)+2,4)
for ibgd=0,n_elements(bgdreg)-1 do begin

   ab=abarr[ibgd]
   projinitbgds,indir,obsid,header,ab,bgddir,pa=pa,clobber=clobber,grxe=grxe

;print,'Pulling BACKSCAL (assumed to be region area as percentage of image) from:'
   print,'Checking that bgd spectrum exists:'
   print,'  '+clspecdir+specname[ibgd]
   if not file_test(clspecdir+specname[ibgd]) then begin
      stop,'Background spectrum does not exist.  Kinda hard to fit nothing...'
;    print,'Background spectrum does not exist.  Creating:'
;    namesplit=strsplit(bgdreg,'.',/extract)
;    spawn,'./getspecnoarf.py '+indir+' '+obsid+' '+namesplit[0]+' '+specdir+' '+ab
   endif

;bgdpha=mrdfits(clspecdir+specname[ibgd],1,hh,/silent)
;bgdbackscl[ibgd]=sxpar(hh,'BACKSCAL')

   bgdmask=reg2mask(dir+'bgdap0'+ab+'.fits',cldir+bgdreg[ibgd])

   for i=0,3 do begin
      fits_read,dir+'det'+str(i)+ab+'im.fits',detim
      fits_read,dir+'bgdap'+str(i)+ab+'.fits',apim
      if grxe then fits_read,dir+'bgdap'+str(i)+ab+'.fits',grxeim
      bgddetfrac[ibgd,i]=total(detim*bgdmask)
      dettot[ibgd,i]=total(detim)
      bgdapfrac[ibgd,i]=total(apim*bgdmask)
      aptot[ibgd,i]=total(apim)
      if grxe then bgdgrxefrac[ibgd,i]=total(grxeim*bgdmask)
      if grxe then grxetot[ibgd,i]=total(grxeim)
   endfor
   bgdbackscl[ibgd]=total(bgddetfrac[ibgd,*])/1000.^2
   detfrac[ibgd,*]=bgddetfrac[ibgd,*]/total(bgddetfrac[ibgd,*])

   if abarr[ibgd] eq 'A' then iab=0 else if abarr[ibgd] eq 'B' then iab=1 else stop,':('
   bgdpt[ibgd,*]=total(pt[*,iab]*detfrac[ibgd,*])
   bgdczt[ibgd,*]=total(czt[*,iab]*detfrac[ibgd,*])

endfor

readcol,auxildir+'ratiosA.dat',aeline,awidth,af0,af1,af2,af3,/silent
readcol,auxildir+'ratiosA.dat',aindex1,aindex2,ab0,ab1,ab2,ab3,aebreak,/silent
aneut=[aeline[n_elements(af0)-1],awidth[n_elements(af0)-1]]
aeline=aeline[0:n_elements(awidth)-3]
awidth=awidth[0:n_elements(awidth)-3]
aifactors_orig=[[af0],[af1],[af2],[af3]]
readcol,auxildir+'ratiosB.dat',beline,bwidth,bf0,bf1,bf2,bf3,/silent
readcol,auxildir+'ratiosB.dat',bindex1,bindex2,bb0,bb1,bb2,bb3,bebreak,/silent
bneut=[beline[n_elements(bf0)-1],bwidth[n_elements(bf0)-1]]
beline=beline[0:n_elements(bwidth)-3]
bwidth=bwidth[0:n_elements(bwidth)-3]
bifactors_orig=[[bf0],[bf1],[bf2],[bf3]]

b0=0
abref=abarr[0]
for ibgd =0, n_elements(abarr) -1 do begin 
  if abarr[ibgd] ne abref and b0 eq 0 then begin
    b0=ibgd
 endif
endfor




; write to file to be read for other programs to create images and spectra
if not keyword_set(paramfile) then paramfile=clspecdir+'bgdfitparams' else $
   paramfile=clspecdir+paramfile

if absave eq 'A' or absave eq 'B' then paramfile=paramfile+absave+'.dat' $
      else paramfile=paramfile+['A.dat','B.dat']

ab=['A', 'B']
for iab=0,1 do begin

if (absave eq 'A' and iab eq 0) or (absave eq 'B' and iab eq 1) or absave eq 'AB' $
      then begin

if iab eq 0 then begin
    eline=aeline 
    ifactors_orig=aifactors_orig
endif else begin
    eline=beline
    ifactors_orig=bifactors_orig
endelse
openw,lun,paramfile[iab],/get_lun


readcol,cldir+paramstem+ab[iab]+'.dat',params,/silent
line=''
openr,lun2,cldir+paramstem+ab[iab]+'.dat',/get_lun
readf,lun2,line
readf,lun2,line
free_lun,lun2
params1=float(strsplit(line,/extract))

; below is ratio of norm to CXB norm assumed in reference ap images times the
;  ratio of the total cts in 10^5s (used to create reference ap images) to
;  the cts of that image that fall in the extraction region

apnorm=(params[0]/0.002353)*32./total(bgdapfrac[iab*b0,*])
fcxbnorm=fltarr(n_elements(params1))
for i=0,n_elements(params1)-1 do $
      fcxbnorm[i]=params1[i]*total(dettot[i+iab*b0,*])/total(bgddetfrac[i+iab*b0,*])
;neutnorm=params[2]*total(dettot[iab*b0,*])/total(bgddetfrac[iab*b0,*])
if grxe then grxenorm=[(params[n_elements(params)-2])/total(bgdgrxefrac[iab*b0,*]),$
      (params[n_elements(params)-1])/total(bgdgrxefrac[iab*b0,*])]

ifactors=ifactors_orig
instrnorms=fltarr(n_elements(eline)+1,4)
for i=0,n_elements(instrnorms[*,0])-1 do $
      instrnorms[i,*]=params[i+3]/total(reform(ifactors[i,*],4)*reform(bgddetfrac[0+iab*b0,*],4)/reform(dettot[0+iab*b0,*],4))*reform(ifactors[i,*],4)
i=n_elements(eline)+1
neutnorm=fltarr(4)
neutnorm=params[2]/total(reform(ifactors[i,*],4)*reform(bgddetfrac[0+iab*b0,*],4)/reform(dettot[0+iab*b0,*],4))*reform(ifactors[i,*],4)

printf,lun,apnorm,$
  '   # Factor relative to reference norm of PL for aperture component',$
  format='(E12.4,A)'
printf,lun,fcxbnorm,$
  '   # Norm of PL in CXB xspec model for entire FOV',$
  format='('+str(n_elements(params1))+'E12.4,A)'
;printf,lun,neutnorm,$
; '   # Norm of PL in xspec model of soft componet for entire FOV',format='(E12.4,A)'
printf,lun,'# Norms of instrumental line components for 4 detectors: 0 1 2 3'
for i=0,n_elements(eline)-1 do $
      printf,lun,instrnorms[i,*],format='(4E12.4)'
printf,lun,'# Norms of instrumental continua for 4 detectors: 0 1 2 3'
printf,lun,instrnorms[n_elements(eline),*],format='(4E12.4)'
printf,lun,'# Norms of steep, low energy component for 4 detectors: 0 1 2 3'
printf,lun,neutnorm,format='(4E12.4)'
if grxe then begin
    printf,lun,'GRXE ',grxenorm[0],'   # GRXE Fe line normalization',$
          format='(A,E12.4,A)'
    printf,lun,'GRXE ',grxenorm[1],'   # GRXE WD normalization',$
          format='(A,E12.4,A)'
endif
free_lun,lun

endif

endfor

end
