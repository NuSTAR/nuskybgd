; NUSKYBGD_SPEC
;
; Read in a parameter file and generate a background spectrum (using fakeit in
; XSPEC) for any arbitrary region.
;
; OPTIONAL KEYWORDS
;
;   paramfile --> If paramdir is set and the parameter filename is the default
;                 (bgdfitparams[A/B].dat), this argument can be omitted.
;                 paramfile must include the full path, or at least the path
;                 relative to your IDL directory -- if set, paramdir is not used
;
;   specname --> Name and path (relative to the event_cl/ directory) of the
;                source spectrum if it does not follow the convention:
;                    event_cl/'+specdir+'/'+specdir+ab+'_sr_g30.pha'
;
;   expfactor --> Factor of the exposure time that the spectrum is simulated for.
;                 The default is 100x.  Smaller values can run into small number
;                 counts issues.  Setting expfactor=1 and /perfect can tell you
;                 the typical number of photons in a given energy range you expect
;                 in your source region.
;
;   perfect --> Suppresses the addition of shot noise to the bgd counts spectrum.
;
;   srcdir --> Relative to event_cl/, the output directory of the bgd spectra
;              if you want a location different from where your source spectrum is.
;
;   fakname --> If you don't want the bgd spectrum to have the same name as
;               the source spectrum (but with 'bgd' prefixed to it), you can
;               specify an alternative name to have 'bgd' prefixed on.
;
;   avgfcxb --> Instead of using the value for the fCXB normalization in the
;               parameter file, use the expected value.  In practice this doesn't
;               work very well, so I would only set this with care.


pro nuskybgd_spec_models,indir,obsid,srcreg,specdir,bgddir,ab,bgddirref, $
      paramfile=paramfile,specname=specname,expfactor=expfactor,grxe=grxe,$
      fakname=fakname,perfect=perfect,srcdir=srcdir,avgfcxb=avgfcxb,$
      forcermf=forcermf,savexcm=savexcm

auxildir=getenv('NUSKYBGD_AUXIL')+'/'
caldbdir=getenv('CALDB')+'/'


bearf = auxildir+'be.arf'
diagrmf=auxildir+'diag.rmf'

dir=indir
if strmid(dir,strlen(dir)-1) ne '/' then dir=dir+'/'
if ab eq 'A' then iab=0 else if ab eq 'B' then iab=1 else stop,'ab defined wrong'
if not keyword_set(expfactor) then expfactor=100.
if keyword_set(perfect) then ctstat='n' else ctstat='y'
if not keyword_set(grxe) then grxe=0

cldir=dir+obsid+'/event_cl/'+specdir+'/'
if size(bgddirref,/type) ne 0 then refdir=dir+obsid+'/event_cl/'+bgddirref+'/' $
      else begin
    if file_test(dir+obsid+'/event_cl/'+bgddir+'/bgdap0'+ab+'.fits') then $
          refdir=dir+obsid+'/event_cl/'+bgddir+'/' else $
        if file_test(dir+obsid+'/event_cl/bgdap0'+ab+'.fits') then $
          refdir=dir+obsid+'/event_cl/' else $
        stop,'Cannot find bgddirref: Please set correctly'
endelse
if not keyword_set(paramfile) then $
      paramfile=dir+obsid+'/event_cl/'+bgddir+'/bgdfitparams'+ab+'.dat'
if not keyword_set(specname) then specname=specdir+ab+'_sr_g30.pha'
;if not keyword_set(srcspec) then begin
;    src=strsplit(srcreg,',',/extract)
;    srcspec=dir+obsid+'/event_cl/'+src[0]+'/'+src[0]+ab+'_sr_g30.pha'
;endif

print,'Pulling BACKSCAL (assumed to be region area as percentage of image) from:'
print,'  '+cldir+specname
if not file_test(cldir+specname) then begin
    stop,'Spectrum does not exist.'
   ; write bit to automatically create, but a few extra steps to do that
endif

pha=mrdfits(cldir+specname,1,hh,/silent)
arf_file=fxpar(hh, 'ANCRFILE')
;backscl=sxpar(hh,'BACKSCAL')

livetime=sxpar(hh,'LIVETIME')

cldir=dir+obsid+'/event_cl/'
if size(srcreg,/type) eq 7 then $
      mask=reg2mask(refdir+'bgdap0'+ab+'.fits',cldir+srcreg) $
; addition of addabs2rmf below means we need the source region
;  else if size(srcreg,/type) eq 2 then mask=srcreg $
  else stop,'  NUSKYBGD_SPEC: Source region name/mask ill-defined.'
print, "Created mask"
backscl=total(mask)/1000.^2

detfrac=fltarr(4)
dettot=fltarr(4)
apreg=0.
grxereg=0.
for i=0,3 do begin
    fits_read,refdir+'det'+str(i)+ab+'im.fits',detim
    fits_read,refdir+'bgdap'+str(i)+ab+'.fits',apim
    if grxe then fits_read,refdir+'bgdgrxe'+str(i)+ab+'.fits',grxeim
    detfrac[i]=total(detim*mask)
    dettot[i]=total(detim)
    apreg+=total(apim*mask)
    if grxe then grxereg+=total(grxeim*mask)
endfor
dettotfrac=total(detfrac)/total(dettot)
detwt=detfrac/total(detfrac)
detfrac=detfrac/dettot

;readcol,auxildir+'ratios_lineE.dat',eline,width,/silent
;readcol,auxildir+'ratios_lineE.dat',blah,index1,ebreak,index2,$
;      format='(A,F,F,F)',/silent
;readcol,auxildir+'ratios'+ab+'.dat',f0,f1,f2,f3,/silent
readcol,auxildir+'ratios'+ab+'.dat',eline,width,f0,f1,f2,f3,/silent
readcol,auxildir+'ratios'+ab+'.dat',index1,index2,b0,b1,b2,b3,ebreak,/silent
;neut=[eline[n_elements(f0)-1],width[n_elements(f0)-1]]
eline=eline[0:n_elements(width)-3]
width=width[0:n_elements(width)-3]

readcol,paramfile,p,/silent
apnorm=p[0]*0.002353*apreg/32.
fcxbnorm=p[1]*dettotfrac
;neutnorm=p[2]*dettotfrac
readcol,paramfile,p0,p1,p2,p3,/silent,skipline=3
pinstr=fltarr(n_elements(p0))

for i=0,n_elements(pinstr)-1 do pinstr[i]=total([p0[i],p1[i],p2[i],p3[i]]*detfrac)
if grxe then begin
    readcol,paramfile,blah,gp,format='(A,F)',/silent
    ii=where(blah eq 'GRXE')
    if n_elements(ii) eq 2 and blah[ii[0]] eq 'GRXE' then $
          grxenorm=[gp[ii[0]],gp[ii[1]]]*grxereg $
    else stop,'NUSKYBGD_SPEC: Problem with GRXE values in parameter file'
endif

;namesplit=strsplit(specname,'.',/extract)
;if strmid(namesplit[0],strlen(namesplit[0])-4) eq '_g30' then $
;      rmfname=strmid(namesplit[0],0,strlen(namesplit[0])-4)+'.rmf' $
;      else rmfname=namesplit[0]+'.rmf'
h=headfits(cldir+specdir+'/'+specname,exten=1,/silent)
rmfname=sxpar(h,'RESPFILE')
if keyword_set(forcermf) then rmfname=forcermf

if file_test(cldir+specdir+'temp.rmf') then spawn,'rm -f '+cldir+specdir+'temp.rmf'
;addabs2rmf,cldir+specdir+'/'+rmfname,ab,refdir+'bgdap0'+ab+'.fits',cldir+srcreg,$
;      cldir+bgddir,cldir+specdir+'/temp.rmf',method=2

;; addabs2rmf,cldir+specdir+'/'+rmfname,ab,refdir+'bgdap0'+ab+'.fits',cldir+srcreg,$
;;       refdir,cldir+specdir+'/bgd_model.rmf',method=2

;; rmfname='bgd_model.rmf'

if not keyword_set(fakname) then fakname=specname
if not keyword_set(srcdir) then srcdir=cldir+specdir+'/' else srcdir=cldir+srcdir+'/'


openw,lun,srcdir+'/'+file_basename(specname, '.pha')+'_bgd_models.xcm',/get_lun
;openw,lun,'temp.xcm',/get_lun

dataset=0
data_ind = iab+1 + dataset * 2.
if ~grxe then begin
   model_ind = 6*(iab+dataset)
endif else begin
   model_ind = 7 *(iab+dataset)
endelse
if model_ind eq 0 then model_ind = 2

spec_ind = str(model_ind, format='(i0)')
data_ind = str(data_ind, format='(i0)')



printf, lun, 'response '+spec_ind+':'+data_ind+' '+rmfname
printf, lun, 'arf '+spec_ind+':'+data_ind+' '+bearf
;printf,lun,'model '+spec_ind+':apbgd'+ab+'_'+data_ind+' (po*highecut)'
;printf,lun,'1.29 -1'
;printf,lun,str(apnorm)+' -1'
;printf,lun,'1e-4 -1'
;printf,lun,'41.13 -1'
printf, lun, 'model '+spec_ind+':apbg'+ab+'_'+data_ind+' cutoffpl'
printf,lun,'1.29 -1'
printf,lun,'41.13 -1'
printf,lun,str(apnorm)+' -1'

model_ind++
spec_ind = str(model_ind, format='(i0)')

printf, lun, 'response '+spec_ind+':'+data_ind+' '+rmfname
printf, lun, 'arf '+spec_ind+':'+data_ind+' '+arf_file
;; printf,lun,'model '+spec_ind+':fxapbgd'+ab+'_'+data_ind+' (po*highecut)'
;; printf,lun,'1.29 -1'
;; if keyword_set(avgfcxb) then printf,lun,str(0.002353*(2.45810736/3600.*1000.)^2* $
;;       backscl)+' -1' else printf,lun,str(fcxbnorm)+' -1'
;; printf,lun,'1e-4 -1'
;; printf,lun,'41.13 -1'


printf,lun,'model '+spec_ind+':fxapbgd'+ab+'_'+data_ind+' cutoffpl'
printf,lun,'1.29 -1'
printf,lun,'41.13 -1'
if keyword_set(avgfcxb) then printf,lun,str(0.002353*(2.45810736/3600.*1000.)^2* $
      backscl)+' -1' else printf,lun,str(fcxbnorm)+' -1'
;printf,lun,'1e-4 -1'




model_ind++
spec_ind = str(model_ind, format='(i0)')

printf, lun, 'response '+spec_ind+':'+data_ind+' '+rmfname
printf, lun, 'arf '+spec_ind+':'+data_ind+' none'
printf,lun,'model '+spec_ind+':line_bgd'+ab+'_'+data_ind+' (',format='($,A)'

;printf,lun,'model nuabs*(',format='($,A)'
for i=0,n_elements(eline)-1 do printf,lun,'lorentz+',format='($,A)'
;printf,lun,'bknpo+phabs*po)'
printf,lun,'apec)'
for i=0,n_elements(eline)-1 do begin
   printf,lun,str(eline[i])+' -1'
   printf,lun,str(width[i])+' -1'
   printf,lun,str(pinstr[i])+ ' -1'
endfor
printf,lun,str(index1[0])+' -1'
printf,lun,str(index2[0])+' -1'
printf,lun,str(ebreak[0])+' -1'
printf,lun,str(pinstr[n_elements(eline)])+' -1'


model_ind++
spec_ind = str(model_ind, format='(i0)')
printf, lun, 'response '+spec_ind+':'+data_ind+' '+diagrmf
printf,lun,'model '+spec_ind+':particle_bgd'+ab+'_'+data_ind+' bknpo'
;printf,lun,'model bknpo'
printf,lun,str(index1[1])+' -1'
printf,lun,str(ebreak[1])+' -1'
printf,lun,str(index2[1])+' -1'
printf,lun,str(pinstr[n_elements(eline)+1])+ ' -1'









close, lun
free_lun, lun

print,'Background model created:'
print,'  '+srcdir+'bgd'+fakname




end
