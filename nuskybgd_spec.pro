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


pro nuskybgd_spec,indir,obsid,srcreg,specdir,bgddir,ab,bgddirref, $
      paramfile=paramfile,specname=specname,expfactor=expfactor,grxe=grxe,$
      fakname=fakname,perfect=perfect,srcdir=srcdir,avgfcxb=avgfcxb

auxildir=getenv('NUSKYBGD_AUXIL')+'/'
caldbdir=getenv('CALDB')+'/'
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
;backscl=sxpar(hh,'BACKSCAL')

livetime=sxpar(hh,'LIVETIME')

cldir=dir+obsid+'/event_cl/'
mask=reg2mask(refdir+'bgdap0'+ab+'.fits',cldir+srcreg)
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

pt=loadnuabs(0)
czt=loadnuabs(1)
spt=total(pt[*,iab]*detwt)
sczt=total(czt[*,iab]*detwt)

namesplit=strsplit(specname,'.',/extract)
if strmid(namesplit[0],strlen(namesplit[0])-4) eq '_g30' then $
      rmfname=strmid(namesplit[0],0,strlen(namesplit[0])-4)+'.rmf' $
      else rmfname=namesplit[0]+'.rmf'

if not keyword_set(fakname) then fakname=specname
if not keyword_set(srcdir) then srcdir=cldir+specdir+'/' else srcdir=cldir+srcdir+'/'

openw,lun,'temp.xcm',/get_lun
printf,lun,'lmod nuabs'
printf,lun,'model nuabs*(po*highecut)'
printf,lun,str(spt)+' -1'
printf,lun,str(sczt)+' -1'
printf,lun,'0. -1'
printf,lun,'0.9 -1'
printf,lun,'1.29 -1'
printf,lun,str(apnorm)
printf,lun,'1e-4 -1'
printf,lun,'41.13 -1'
;printf,lun,'4.8 -1'
;printf,lun,str(neutnorm)
spawn,'rm -f '+srcdir+'/bgdap'+fakname
printf,lun,'fakeit none & '+cldir+specdir+'/'+rmfname+' & '+auxildir+$
      'be.arf & '+ctstat+' &  & '+srcdir+'/bgdap'+fakname+' & '+$
      str(livetime*expfactor)
printf,lun,'data none'

printf,lun,'model nuabs*(po*highecut)'
printf,lun,str(spt)+' -1'
printf,lun,str(sczt)+' -1'
printf,lun,'0. -1'
printf,lun,'0.9 -1'
printf,lun,'1.29 -1'
if keyword_set(avgfcxb) then printf,lun,str(0.002353*(2.45810736/3600.*1000.)^2* $
      backscl) else printf,lun,str(fcxbnorm)
printf,lun,'1e-4 -1'
printf,lun,'41.13 -1'
spawn,'rm -f '+srcdir+'/bgdfcxb'+fakname
printf,lun,'fakeit none & '+cldir+specdir+'/'+rmfname+' & '+$
;      caldbdir+'data/nustar/fpm/bcf/arf/nu'+ab+'20100101v004.arf & '+$
      auxildir+'fcxb'+ab+'.arf & '+$
      ctstat+' &  & '+srcdir+'/bgdfcxb'+fakname+' & '+$
      str(livetime*expfactor)
printf,lun,'data none'

printf,lun,'model nuabs*(',format='($,A)'
for i=0,n_elements(eline)-1 do printf,lun,'lorentz+',format='($,A)'
;printf,lun,'bknpo+phabs*po)'
printf,lun,'apec)'
printf,lun,str(spt)+' -1'
printf,lun,str(sczt)+' -1'
printf,lun,'0. -1'
printf,lun,'0.9 -1'
for i=0,n_elements(eline)-1 do begin
    printf,lun,str(eline[i])+' -1'
    printf,lun,str(width[i])+' -1'
    printf,lun,str(pinstr[i])
endfor
printf,lun,str(index1[0])+' -1'
printf,lun,str(index2[0])+' -1'
printf,lun,str(ebreak[0])+' -1'
printf,lun,str(pinstr[n_elements(eline)])
;printf,lun,str(neut[0])
;printf,lun,str(neut[1])
;printf,lun,str(pinstr[n_elements(eline)+1])
spawn,'rm -f '+srcdir+'/bgdinstr'+fakname
printf,lun,'fakeit none & '+cldir+specdir+'/'+rmfname+' &  & '+ctstat+' &  & '+$
      srcdir+'/bgdinstr'+fakname+' & '+str(livetime*expfactor)
printf,lun,'data none'

printf,lun,'model bknpo'
printf,lun,str(index1[1])+' -1'
printf,lun,str(ebreak[1])+' -1'
printf,lun,str(index2[1])+' -1'
printf,lun,str(pinstr[n_elements(eline)+1])
printf,lun,'fakeit none & '+auxildir+'diag.rmf &  & '+ctstat+' &  & '+$
      srcdir+'/bgdintcont'+fakname+' & '+str(livetime*expfactor)
printf,lun,'data none'

if grxe then begin

printf,lun,'model nuabs*(gauss+atable{'+auxildir+'polarmodel.fits})'
printf,lun,str(spt)+' -1'
printf,lun,str(sczt)+' -1'
printf,lun,'0. -1'
printf,lun,'0.9 -1'
printf,lun,'6.7 -1'
printf,lun,'0.0 -1'
printf,lun,str(grxenorm[0])
printf,lun,'0.6 -1'
printf,lun,str(grxenorm[1])
spawn,'rm -f '+srcdir+'/bgdgrxe'+fakname
printf,lun,'fakeit none & '+cldir+specdir+'/'+rmfname+' & '+auxildir+$
      'be.arf & '+ctstat+' &  & '+srcdir+'/bgdgrxe'+fakname+' & '+$
      str(livetime*expfactor)

endif

printf,lun,'exit'
free_lun,lun

f = file_info(srcdir+'/bgdintcont'+fakname)
IF f.exists THEN spawn, 'rm '+srcdir+'/bgdintcont'+fakname

spawn,'xspec - temp.xcm'
spawn,'rm -f temp.xcm'



spawn,'fparkey '+str(backscl)+' '+srcdir+'/bgdap'+fakname+' BACKSCAL'
spawn,'fparkey '+str(backscl)+' '+srcdir+'/bgdfcxb'+fakname+' BACKSCAL'
spawn,'fparkey '+str(backscl)+' '+srcdir+'/bgdinstr'+fakname+' BACKSCAL'
spawn,'fparkey '+str(backscl)+' '+srcdir+'/bgdintcont'+fakname+' BACKSCAL'
if grxe then spawn,'fparkey '+str(backscl)+' '+srcdir+'/bgdgrxe'+fakname+' BACKSCAL'

spawn,'pwd > '+srcdir+'/temp'
cd,srcdir
thisdir='' ;cldir+specdir+'/'
if not grxe then begin
  
      print,'mathpha '+thisdir+'bgdap'+fakname+'+'+thisdir+'bgdfcxb'+fakname+$
            '+'+thisdir+'bgdinstr'+fakname+'+'+thisdir+'bgdintcont'+fakname+$
            ' COUNTS '+thisdir+'bgd'+fakname+' '+$
            str(livetime*expfactor)+' 1.0 0 clobber=yes'


      spawn,'mathpha '+thisdir+'bgdap'+fakname+'+'+thisdir+'bgdfcxb'+fakname+$
            '+'+thisdir+'bgdinstr'+fakname+'+'+thisdir+'bgdintcont'+fakname+$
            ' COUNTS '+thisdir+'bgd'+fakname+' '+$
            str(livetime*expfactor)+' 1.0 0 clobber=yes'

      print, expfactor*livetime

      spawn,'mathpha '+$
            thisdir+'bgdinstr'+fakname+'+'+thisdir+'bgdintcont'+fakname+$
            ' COUNTS '+thisdir+'bgd_allint_'+fakname+' '+$
            str(livetime*expfactor)+' 1.0 0 clobber=yes'


   ENDIF ELSE  spawn,'mathpha '+thisdir+'bgdap'+fakname+'+'+thisdir+'bgdfcxb'+fakname+$
                     '+'+thisdir+'bgdinstr'+fakname+'+'+thisdir+'bgdintcont'+fakname+$
                     '+'+thisdir+'bgdgrxe'+fakname+$
                     ' COUNTS '+thisdir+'bgd'+fakname+' '+$
                     str(livetime*expfactor)+' 1.0 0 clobber=yes'
spawn,'fparkey '+str(backscl)+' '+'bgd'+fakname+' BACKSCAL'
print,'fparkey '+str(backscl)+' '+'bgd'+fakname+' BACKSCAL'
spawn,'fparkey bgd'+fakname+' '+specname+' BACKFILE'
print,'fparkey bgd'+fakname+' '+specname+' BACKFILE'



spawn,'fparkey '+str(backscl)+' '+'bgd_allint_'+fakname+' BACKSCAL'
print,'fparkey '+str(backscl)+' '+'bgd_allint_'+fakname+' BACKSCAL'


line=''
openr,lun,'temp',/get_lun
readf,lun,line
free_lun,lun
spawn,'rm -f temp'
cd,line


print,str(expfactor)+'x faked background spectra created:'
print,'  '+srcdir+'bgd'+fakname


end
