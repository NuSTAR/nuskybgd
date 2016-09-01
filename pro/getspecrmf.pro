pro getspecrmf,evtfile,ab,regions,bgddir,outdir,outbase,$
      cmprmf=cmprmf,grp=grp,$
      nospec=nospec,normf=normf,outspec=outspec,$
      imcheck=imcheck,e1=e1,e2=e2

if size(grp,/type) eq 0 then grp=3
if not file_test(outdir,/directory) then spawn,'mkdir '+outdir
if size(cmprmf,/type) eq 0 then cmprmf=1e-9 else if cmprmf eq 1 then $
      cmprmf=1e-6 else if cmprmf gt 1e-4 then $
      stop,'GETSPECRMF: CMPRMF value too large, set to 1 for value of 1e-6'
refspechdr=getenv('NUSKYBGD_AUXIL')+'/spechdr.txt'

; Make header file so regions can be read in correctly

header=headfits(evtfile,exten=0,/silent)
baseheader=headfits(evtfile,exten=1,/silent)
exptime=sxpar(baseheader,'EXPOSURE')
sxaddpar,header,'EXPOSURE',exptime

xtt='X' & ytt='Y' & i=1 & check=1
while i lt 100 and check ne 0 do begin
    if str(sxpar(baseheader,'TTYPE'+str(i))) eq xtt and $
          str(sxpar(baseheader,'TTYPE'+str(i+1))) eq ytt then check=0
    i++
endwhile
xx=str(i-1)
yy=str(i)

sxaddpar,header,'BITPIX',32
sxaddpar,header,'NAXIS',2
sxaddpar,header,'NAXIS1',$
      fix(sxpar(baseheader,'TLMAX'+xx))-fix(sxpar(baseheader,'TLMIN'+xx))+1
sxaddpar,header,'NAXIS2',$
      fix(sxpar(baseheader,'TLMAX'+yy))-fix(sxpar(baseheader,'TLMIN'+yy))+1
sxaddpar,header,'CTYPE1',sxpar(baseheader,'TCTYP'+xx)
sxaddpar,header,'CTYPE2',sxpar(baseheader,'TCTYP'+yy)
sxaddpar,header,'CRPIX1',sxpar(baseheader,'TCRPX'+xx)
sxaddpar,header,'CRPIX2',sxpar(baseheader,'TCRPX'+yy)
sxaddpar,header,'CRVAL1',sxpar(baseheader,'TCRVL'+xx)
sxaddpar,header,'CRVAL2',sxpar(baseheader,'TCRVL'+yy)
sxaddpar,header,'CDELT1',sxpar(baseheader,'TCDLT'+xx)
sxaddpar,header,'CDELT2',sxpar(baseheader,'TCDLT'+yy)

; Grab Reference Files

if not keyword_set(nospec) then begin
    evts=mrdfits(evtfile,1,/silent)
    undefine,hdr
    line=''
    openr,lun,refspechdr,/get_lun
    while ~eof(lun) do begin
        readf,lun,line
        push,hdr,line
    endwhile
    free_lun,lun
endif

if not keyword_set(normf) then begin
    bgddet=fltarr(1000,1000,4)
    for i=0,3 do begin
        file=bgddir+'/det'+str(i)+ab+'im.fits'
        if file_test(file) then fits_read,file,im $
              else stop,'GETSPECRMF: Det image file '+file+' not found.'
        bgddet[*,*,i]=im
    endfor
endif

; Convert region files to region image if necessary

if size(regions,/type) eq 7 then begin
    regim=fltarr(1000,1000)
    fits_write,'temp.fits',regim,header
    for i=0,n_elements(regions)-1 do begin
        mask=reg2mask('temp.fits',regions[i])
        ii=where(mask gt 0.5)
        regim[ii]=i+1
    endfor
    spawn,'rm -f temp.fits'
endif else regim=regions
nreg=max(regim)

; Grab indices of all regions

undefine,rr
undefine,nrr
for r=1,nreg do begin
    ii=where(regim eq r)
    if ii[0] eq -1 then stop,'regim needs to be contiguous 1 -> N_regions'
    push,rr,ii
    push,nrr,n_elements(ii)
endfor

; For each region, extract the spectrum and/or RMF

r0=0
for r=0,nreg-1 do begin

  if not keyword_set(nospec) then begin

    specstr=replicate({CHANNEL:0, COUNTS:0, QUALITY:0, GROUPING:1},4096)
    specstr.channel=indgen(4096)
    ii2d=array_indices(fltarr(1000,1000),rr[r0:r0+nrr[r]-1])
    for i=0,n_elements(ii2d[0,*])-1 do begin
        ii=where(evts.x-1 ge ii2d[0,i] and evts.x-1 lt ii2d[0,i]+1 and $
              evts.y-1 ge ii2d[1,i] and evts.y-1 lt ii2d[1,i]+1 and $
              evts.grade ge 0 and evts.grade le 26)
        if ii[0] ne -1 then for j=0,n_elements(ii)-1 do $
              specstr[evts[ii[j]].pi].counts++
    endfor
    bin=0
    first=1
    for i=0,n_elements(specstr)-1 do begin
        if not first then specstr[i].grouping=-1
        bin+=specstr[i].counts
        if bin lt grp then first=0 else begin
            first=1
            bin=0
        endelse
    endfor
    outspec=specstr.counts
    if keyword_set(imcheck) then begin
        ii=where(regim0 eq r+1)
        totim=total(imcheck[ii])
        totspec=total(outspec[round((e1-1.6)/0.04+1):round((e2-1.6)/0.04)-1])
        if totim ne totspec then begin
            print,totim,totspec
            stop,'GETSPECRMF: Spectrum =/= imcheck'
        endif
    endif
    sxaddpar,hdr,'EXPOSURE',exptime
    sxaddpar,hdr,'LIVETIME',exptime
    sxaddpar,hdr,'INSTRUME','FPM'+ab
    sxaddpar,hdr,'BACKSCAL',nrr[r]/1000./1000.
    if not keyword_set(normf) then sxaddpar,hdr,'RESPFILE',$
          outbase+ab+str(r+1)+'.rmf'
    mwrfits,specstr,outdir+'/'+outbase+ab+str(r+1)+'.pha',hdr,/create,/silent
    print,'Spectrum '+outbase+ab+str(r+1)+'.pha created'

  endif

  if not keyword_set(normf) then begin

    mask=intarr(1000,1000)
    ii=rr[r0:r0+nrr[r]-1]
    mask[ii]=1
    bgdval=fltarr(4)
    for i=0,3 do bgdval[i]=total(mask*bgddet[*,*,i])
    bgdval/=total(bgdval)
    ii=where(bgdval lt 0.01)
    bgdval[ii]=0.
    ii=where(bgdval ge 0.01)
    if ii[0] eq -1 then $
          stop,'GETSPECRMF: Region does not overlap with any detectors'
    bgdval/=total(bgdval)
    matrix=fltarr(4096,4096)
;    abs=fltarr(4096)
;    abs[*]=1.
    rmf1str=mrdfits(getcaldbfile('rmf',ab,0),2,r2h,/silent)
    for i=0,n_elements(ii)-1 do begin
        absstr=mrdfits(getcaldbfile('detabs',ab,ii[i]),ii[i]+1,dh,/silent)
        abs=absstr.detabs
;        abs+=absstr.detabs*bgdval[ii[i]]
        rmfstr=mrdfits(getcaldbfile('rmf',ab,ii[i]),1,rh,/silent)
        thismatrix=rmfstr.matrix
        for e=0,4095 do thismatrix[*,e]=thismatrix[*,e]*abs[e]
        matrix+=thismatrix*bgdval[ii[i]]
    endfor
    mwrfits,rmf1str,outdir+'/'+outbase+ab+str(r+1)+'_temp.rmf',r2h,/silent,/create
    rmfstr.matrix=matrix
    sxaddpar,rh,'LO_THRES',0.0,'modified by cmprmf'
    mwrfits,rmfstr,outdir+'/'+outbase+ab+str(r+1)+'_temp.rmf',rh,/silent

    if cmprmf then spawn,'cmprmf '+outdir+'/'+outbase+ab+str(r+1)+'_temp.rmf '+ $
          outdir+'/'+outbase+ab+str(r+1)+'.rmf '+str(cmprmf)+' clobber=yes'
    print,'RMF '+outbase+ab+str(r+1)+'.rmf created'
    spawn,'rm -f '+outdir+'/'+outbase+ab+str(r+1)+'_temp.rmf'

  endif

  r0=r0+nrr[r]

endfor


end
