function reg2mask,fitsfile,regfile,outfile=outfile,notype=notype


fits_read,fitsfile,im,header
xpix0 = sxpar(header,'CRPIX1')
ypix0 = sxpar(header,'CRPIX2')
xdeg0 = sxpar(header,'CRVAL1')
ydeg0 = sxpar(header,'CRVAL2')
dx = sxpar(header,'CDELT1')
dy = sxpar(header,'CDELT2')
nx=n_elements(im[*,0])
ny=n_elements(im[0,*])
undefine,params
undefine,nparams
undefine,inclexcl
undefine,type

;print,"Reading ",regfile

openr,lun,regfile,/get_lun
line=''
if not keyword_set(notype) then begin
    while ~ eof(lun) do begin
        readf,lun,line
        if strcmp(str(line),'fk5') or strcmp(str(line),'image') or $
              strcmp(str(line),'physical') then break
    endwhile
    if (not strcmp(str(line),'fk5')) and (not strcmp(str(line),'image')) and $
              (not strcmp(str(line),'physical')) then begin
        free_lun,lun
    stop,'REG2MASK: Region files must be in FK5, IMAGE, or PHYSICAL units'
    endif
    if line eq 'fk5' then wcs=1 else wcs=0
endif else begin
    for blah=1,notype do readf,lun,line
    wcs=1
endelse
while ~ eof(lun) do begin
    readf,lun,line
    t0=strsplit(line,')',/extract)
    t1=strsplit(t0[0],'(',/extract)
    if strmid(t1[0],0,1) ne '-' then begin
        push,inclexcl,0 
        push,type,t1[0]
    endif else begin
        push,inclexcl,1
        push,type,strmid(t1[0],1)
    endelse
    t2=t1[1]
    t3=strsplit(t2,',',/extract)
    push,nparams,n_elements(t3)
    for i=0,n_elements(t3)-1 do begin
;        if i eq n_elements(t3)-1 then t3[i]=strmid(t3[i],0,strlen(t3[i])-1)
        push,params,t3[i]
    endfor
endwhile
free_lun,lun

inclmask=intarr(nx,ny)
exclmask=intarr(nx,ny)
for ir=0,n_elements(type)-1 do begin
  if ir ne 0 then np+=nparams[ir-1] else np=0
  case type[ir] of 
    'polygon': begin
      xy=fltarr(nparams[ir]/2.,2)
      for i=0,nparams[ir]/2.-1 do begin
        if wcs then begin
          sexi=strsplit(params[i*2+np],':',/extract)
          if n_elements(sexi) gt 1 then $
              ra=15*ten(float(sexi[0]),float(sexi[1]),float(sexi[2])) $
                 else ra=float(params[i*2+np])
          sexi=strsplit(params[i*2+1+np],':',/extract)
          if n_elements(sexi) gt 1 then $
              dec=ten(float(sexi[0]),float(sexi[1]),float(sexi[2])) $
                 else dec=float(params[i*2+1+np])
          adxy,header,ra,dec,xx,yy
          xx=xx+1
          yy=yy+1
        endif else begin
          xx=float(params[i*2+np])
          yy=float(params[i*2+1+np])
        endelse
        xy[i,*]=[xx,yy]
      endfor
      end
    'circle': begin
      if wcs then begin
        sexi=strsplit(params[0+np],':',/extract)
        if n_elements(sexi) gt 1 then $
            ra=15*ten(float(sexi[0]),float(sexi[1]),float(sexi[2])) $
               else ra=float(params[0+np])
        sexi=strsplit(params[1+np],':',/extract)
        if n_elements(sexi) gt 1 then $
            dec=ten(float(sexi[0]),float(sexi[1]),float(sexi[2])) $
               else dec=float(params[1+np])
        adxy,header,ra,dec,x,y
        unit=strmid(params[2+np],strlen(params[2+np])-1,1)
        r=float(strmid(params[2+np],0,strlen(params[2+np])-1))
        if unit eq '"' then r=r/3600. else if unit eq "'" then r=r/60. $
              else stop,'REG2MASK: arcmin or arcsec unit expected.'
        r=r/dy
        x=x+1
        y=y+1
      endif else begin
        x=float(params[0+np])
        y=float(params[1+np])
        r=float(params[2+np])
      endelse
      ang=findgen(72)*5.*!pi/180.
      xy=fltarr(n_elements(ang),2)
      for k=0,n_elements(ang)-1 do xy[k,*]=r*[cos(ang[k]),sin(ang[k])]+[x,y]
      end
    'ellipse': begin
      if wcs then begin
        sexi=strsplit(params[0+np],':',/extract)
        if n_elements(sexi) gt 1 then $
            ra=15*ten(float(sexi[0]),float(sexi[1]),float(sexi[2])) $
               else ra=float(params[0+np])
        sexi=strsplit(params[1+np],':',/extract)
        if n_elements(sexi) gt 1 then $
            dec=ten(float(sexi[0]),float(sexi[1]),float(sexi[2])) $
               else dec=float(params[1+np])
        adxy,header,ra,dec,x,y
        unit=strmid(params[2+np],strlen(params[2+np])-1,1)
        r1=float(strmid(params[2+np],0,strlen(params[2+np])-1))
        if unit eq '"' then r1=r1/3600. else if unit eq "'" then r1=r1/60. $
              else stop,'REG2MASK: arcmin or arcsec unit expected.'
        r1=r1/dy
        unit=strmid(params[3+np],strlen(params[3+np])-1,1)
        r2=float(strmid(params[3+np],0,strlen(params[3+np])-1))
        if unit eq '"' then r2=r2/3600. else if unit eq "'" then r2=r2/60. $
              else stop,'REG2MASK: arcmin or arcsec unit expected.'
        r2=r2/dy
        x=x+1
        y=y+1
      endif else begin
        x=float(params[0+np])
        y=float(params[1+np])
        r1=float(params[2+np])
        r2=float(params[3+np])
      endelse
      rang=float(params[4+np])*!pi/180.
      ang=findgen(72)*5.*!pi/180.
      xy=fltarr(n_elements(ang),2)
      for k=0,n_elements(ang)-1 do $
            xy[k,*]=[[cos(rang),sin(rang)],[-sin(rang),cos(rang)]]#$
                [r1*cos(ang[k]),r2*sin(ang[k])]+[x,y]
      end
    'box': begin
      if wcs then begin
        sexi=strsplit(params[0+np],':',/extract)
        if n_elements(sexi) gt 1 then $
            ra=15*ten(float(sexi[0]),float(sexi[1]),float(sexi[2])) $
               else ra=float(params[0+np])
        sexi=strsplit(params[1+np],':',/extract)
        if n_elements(sexi) gt 1 then $
            dec=ten(float(sexi[0]),float(sexi[1]),float(sexi[2])) $
               else dec=float(params[1+np])
        adxy,header,ra,dec,x,y
        unit=strmid(params[2+np],strlen(params[2+np])-1,1)
        r1=float(strmid(params[2+np],0,strlen(params[2+np])-1))
        if unit eq '"' then r1=r1/3600. else if unit eq "'" then r1=r1/60. $
              else stop,'REG2MASK: arcmin or arcsec unit expected.'
        r1=r1/dy
        unit=strmid(params[3+np],strlen(params[3+np])-1,1)
        r2=float(strmid(params[3+np],0,strlen(params[3+np])-1))
        if unit eq '"' then r2=r2/3600. else if unit eq "'" then r2=r2/60. $
              else stop,'REG2MASK: arcmin or arcsec unit expected.'
        r2=r2/dy
        x=x+1
        y=y+1
      endif else begin
        x=float(params[0+np])
        y=float(params[1+np])
        r1=float(params[2+np])
        r2=float(params[3+np])
      endelse
      rang=float(params[4+np])*!pi/180.
      xy=fltarr(4,2)
      xy[0,*]=[[cos(rang),sin(rang)],[-sin(rang),cos(rang)]]#[-r1/2.,r2/2.]+[x,y]
      xy[1,*]=[[cos(rang),sin(rang)],[-sin(rang),cos(rang)]]#[r1/2.,r2/2.]+[x,y]
      xy[2,*]=[[cos(rang),sin(rang)],[-sin(rang),cos(rang)]]#[r1/2.,-r2/2.]+[x,y]
      xy[3,*]=[[cos(rang),sin(rang)],[-sin(rang),cos(rang)]]#[-r1/2.,-r2/2.]+[x,y]
      end
    else: stop,'REG2MASK: Region type '+type[ir]+' not supported.'
  endcase
  ii=polyfillv(xy[*,0],xy[*,1],nx,ny)
  if inclexcl[ir] then exclmask[ii]=1 else inclmask[ii]=1
endfor

ii=where(exclmask gt 0.5)
if ii[0] ne -1 then inclmask[ii]=0

if keyword_set(outfile) then fits_write,outfile,inclmask,header
return,inclmask

end
