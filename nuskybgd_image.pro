; NUSKYBGD_IMAGE
;
; Takes the nominal aperture bgd shape and flat, detector-dependent shapes for
;   the other 3 bgd components and creates bgd and bgd-subtracted images in
;   sky coordinates.
;
; See nuskybgd.pro for typical inputs.
;
; OPTIONAL KEYWORDS
;
;   paramfile --> If paramdir is set and the parameter filename is the default
;                 (bgdfitparams[A/B].dat), this argument can be omitted.
;                 paramfile must include the full path, or at least the path
;                 relative to your IDL directory -- if set, paramdir is not used
;
;   pa --> user set PA: Only used if reference background images not yet created,
;          or clobber keyword set.  If your PA follows the old convention, you
;          need to run projinitbgds.pro with /oldpa set instead.
;
;   noremakeinstr --> If you run this multiple times without changing the model 
;                     parameters for the instrumental component, you can suppress
;                     the recreation of the instrumental reference spectrum by
;                     setting this keyword.
;
;   instrnorms --> A 4-element array allowing the normalization of the instrumental 
;                  background in each chip to be varied -- to set this component
;                  10% higher than assumed by the parameter file values, set
;                  instrnorms=[1.1,1.1,1.1,1.1]
;
;   normap/normfcxb/normneut --> Scalar adjustments (relative to 1) to the
;                                normalizations of the 3 other components.
;
;   outim --> Alternative name for the output bgd image.
;
;   nobsub --> Set to suppress the creation of a bgd-subtracted image


pro nuskybgd_image,indir,obsid,imname,elow,ehigh,ab,bgddir,paramdir,paramfile,$
                pa=pa,noremakeinstr=noremakeinstr,instrnorms=instrnorms,$
                normap=normap,normfcxb=normfcxb,normneut=normneut,normgrxe=normgrxe,$
                outim=outim,nobsub=nobsub,clobber=clobber,grxe=grxe
                
auxildir=getenv('NUSKYBGD_AUXIL')+'/'

dir=indir
if strmid(dir,strlen(dir)-1) ne '/' then dir=dir+'/' 
cldir=dir+obsid+'/event_cl/'
if not size(bgddir,/type) then bgddir=''
dir=cldir+bgddir+'/'
if not keyword_set(grxe) then grxe=0

if not keyword_set(outim) then outim=cldir+'bgd'+imname $
      else outim=cldir+outim

if not size(paramfile,/type) then $
      paramfile=cldir+paramdir+'/bgdfitparams'+ab+'.dat' $
  else paramfile=paramdir+'/'+paramfile
instrspec=dir+'instrbgdfit'
if size(noremakeinstr,/type) eq 0 then noremakeinstr=0
if not noremakeinstr or not file_test(instrspec+ab+'_det0.pha') then $
      fakinstr,indir,obsid,ab,paramfile,bgddir
if not keyword_set(instrnorms) then instrnorms=[1.0,1.0,1.0,1.0]

usernorm=[1.0,1.0,1.0,1.0]
if size(normap,/type) ne 0 then usernorm[0]=normap
if size(normfcxb,/type) ne 0 then usernorm[1]=normfcxb
if size(normneut,/type) ne 0 then usernorm[2]=normneut
if size(normgrxe,/type) ne 0 then usernorm[3]=normgrxe

fits_read,cldir+imname,data,header
livetime=sxpar(header,'LIVETIME')
projinitbgds,indir,obsid,header,ab,bgddir,pa=pa,clobber=clobber

bgdim=fltarr(1000,1000)
for j=0,3 do begin
    fits_read,dir+'det'+str(j)+ab+'im.fits',detim
    if j eq 0 then detall=detim else detall+=detim
;    npix=sxpar(shead,'CORRSCAL')
    npix=total(detim)
    instrall=mrdfits(instrspec+ab+'_det'+str(j)+'.pha',1,shead,/silent)
    check=tag_names(instrall)
    if check[1] eq 'COUNTS' then tinstr=1e9 $
          else if check[1] eq 'RATE' then tinstr=1.0 $
          else stop,'What kind of infernal internal bgd spectrum did you make???'
    cmd='bgdim+=total(instrall[where(instrall.channel ge '+ $
          '(elow-1.6)/0.04+1 and instrall.channel le '+ $
          '(ehigh-1.6)/0.04)].'+check[1]+') * detim * '+ $
          'livetime/tinstr*instrnorms[j]/npix'
    check=execute(cmd)
    if not check then stop,'Stupid execute command being stupid.'
endfor

readcol,paramfile,norm0,norm1,norm2,norm3,/silent
intcont=[norm0[n_elements(norm1)-1],norm1[n_elements(norm1)-1],$
      norm2[n_elements(norm1)-1],norm3[n_elements(norm1)-1]]
readcol,auxildir+'ratios'+ab+'.dat',index1,index2,b0,b1,b2,b3,ebreak,/silent
if elow lt ebreak[1] and ehigh le ebreak[1] then begin
    term=(ehigh^(1.-index1[1])-elow^(1.-index1[1]))/(1.-index1[1])
endif else if elow ge ebreak[1] then begin
    ind=index2[1]-index1[1]
    term=(ehigh^(1.-index2[1])-elow^(1.-index2[1]))/(1.-index2[1])*ebreak[1]^ind
endif else begin
    ind=index2[1]-index1[1]
    term=(ebreak[1]^(1.-index1[1])-elow^(1.-index1[1]))/(1.-index1[1])
    term+=(ehigh^(1.-index2[1])-ebreak[1]^(1.-index2[1]))/(1.-index2[1])*$
          ebreak[1]^ind
endelse
for j=0,3 do begin
    fits_read,dir+'det'+str(j)+ab+'im.fits',detim
    npix=total(detim)
    bgdim+=intcont[j]*term*detim*livetime/npix
endfor

readcol,paramfile,norm,/silent
norm[1]=norm[1]/0.002353
prefixes=['aperbgd','fcxbbgd']  ;,'neutbgd']
if grxe then begin
    prefixes=[prefixes,'grxe1bgd','grxe2bgd']
    readcol,paramfile,blah,gp,format='(A,F)',/silent
    ii=where(blah eq 'GRXE')
    if n_elements(ii) eq 2 then grxenorm=[gp[ii[0]],gp[ii[1]]] $
      else stop,'NUSKYBGD_SPEC: Problem with GRXE values in parameter file'
endif
for k=0,n_elements(prefixes)-1 do begin
    for j=0,3 do begin
        spec=mrdfits(auxildir+prefixes[k]+ab+'_det'+str(j)+'.pha',$
              1,shead,/silent)
        check=tag_names(spec)
        if check[1] eq 'COUNTS' then tinstr=1.0 $
          else if check[1] eq 'RATE' then tinstr=1e9 $
          else stop,'What kind of infernal ref. bgd spectrum did you make???'
        cmd='refcts=total(spec[where(spec.channel ge (10.0-1.6)/0.04+1 '+ $
              'and spec.channel le (15.0-1.6)/0.04)].'+check[1]+')*tinstr'
        check2=execute(cmd)
        if not check2 then stop,'Stupid execute command being stupid.'
        cmd='cts=total(spec[where(spec.channel ge (elow-1.6)/0.04+1 '+ $
              'and spec.channel le (ehigh-1.6)/0.04)].'+check[1]+')*tinstr'
        check2=execute(cmd)
        if not check2 then stop,'Stupid execute command being stupid.'
        if k eq 0 then begin
            fits_read,dir+'bgdap'+str(j)+ab+'.fits',img
            bgdim+=cts/refcts*norm[k]*img*(livetime/1e5)*usernorm[k]
        endif else if k gt 1 then begin
            fits_read,dir+'bgdgrxe'+str(j)+ab+'.fits',img
            bgdim+=cts/refcts*norm[k]*img*(livetime/1e5)*usernorm[k]
        endif else begin
            fits_read,dir+'det'+str(j)+ab+'im.fits',img
            bgdim+=norm[k]*img*cts/1e9*livetime/total(detall)*usernorm[k]
        endelse
    endfor
endfor

fits_read,cldir+imname,data,header
fits_write,outim,bgdim,header
if not keyword_set(nobsub) then begin
    namesplit=strsplit(imname,'.',/extract)
                                ; Assume that the last . is one
                                ; assocaited with .img or .fits, but
                                ; allow other "."
    newname = cldir
    nfields = n_elements(namesplit)
    newname += namesplit[0]
    FOR i = 1, nfields -2 DO newname+='.'+namesplit[i]
    newname+='bsub.'+namesplit[nfields-1]


    fits_write, newname, data-bgdim, header

endif

end
