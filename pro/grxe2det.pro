pro grxe2det,indir,obsid,header,ab,bgddir,pa,clobber=clobber

common caldb, caldbmodule, caldbinstrument, caldbmetrology, caldbOA

auxildir=getenv('NUSKYBGD_AUXIL')+'/'

dir=indir
if strmid(dir,strlen(dir)-1) ne '/' then dir=dir+'/'
cldir=dir+obsid+'/event_cl/'
if size(bgddir,/type) eq 0 then dir=cldir else begin
    dir=cldir+bgddir+'/'
    if not file_test(dir,/directory) then spawn,'mkdir '+dir
endelse

if ab eq 'A' then imod=0 else if ab eq 'B' then imod=1 else $
      stop,'GRXE2DET: ab defined stupidly, must be either A or B'

if not keyword_set(clobber) then clobber=0
;imexist=0
;for det=0,3 do imexist+=file_test(dir+'bgdgrxe'+str(det)+ab+'.fits')

;if clobber ne 0 or imexist ne 4 then begin
if clobber ne 0 or not file_test(dir+'grxe'+ab+'_det1.img') then begin


; amount of blurring of optics bench edge to account for mass motions; not
;  a significant effect
blursig=10.

; how many det1 pixels to bin up to speed up run time of code
nbinpix=4

; dimensions (x & y) of sky plane; the optics bench is placed in its center
dim=2000

; align GXRE image given orientation of telescopes

fits_read,auxildir+'modelrxte_04deg.fits',grxeim,hdr

; numbers here taken from Roman, may want to update later
dr = !PI/180.
rot_mat = [ [ cos(pa*dr), -sin(pa*dr)], $
              [ sin(pa*dr),  cos(pa*dr)] ]
crpix=[38.,38.]
cdelt=[12.3/3600.0, 12.3/3600]
crval=[sxpar(header,'RA_OBJ'),sxpar(header,'DEC_OBJ')]
make_astr, astr, CD=rot_mat, DELTA = cdelt, CRPIX = crpix, CRVAL = crval
xy2ad,31.25,31.25,astr,center_ra,center_dec

delt=(1.0/9907.6)/dr
crpix=[1000.0,1000.0]
crval=[center_ra,center_dec]
make_astr,astr_optim, CD=rot_mat, DELT = [delt, delt], CRPIX = crpix, CRVAL = crval
mkhdr, hdr_optim, 4, [dim,dim]
putast, hdr_optim, astr_optim, CD_TYPE=2

fits_read,auxildir+'modelrxte_04deg.fits',grxeim,hdr
EXTAST, hdr, astr

if clobber ne 0 or not file_test(dir+'grxe_image.fits') then begin
    ; map GRXE fits image to our visible sky
    print
    print,'Rotating GRXE image for your obs -- this will take awhile'
    perc=dim/100
    grxe_image=fltarr(dim,dim)
    for ii=0,dim-1 do begin
       if ii mod perc eq 0 then print,str(ii/perc)+'%..',format='($,A)'
       for jj=0,dim-1 do begin
          xy2ad,ii,jj,astr_optim,ra,dec
          euler, ra, dec, lon, lat, 1
          ad2xy,lon,lat,astr,xx,yy
          grxe_image[ii,jj]=grxeim[xx,yy]
       endfor
    endfor
    print
    writefits, dir+'grxe_image.fits',grxe_image,hdr_optim
endif else fits_read,dir+'grxe_image.fits',grxe_image,hdr_optim

; outline of optics module
omoutline = [[-64,-176,-237,-234,-118,173,334,626,747,748,689,570,-64], $
             [-240,-140,-44,92,205,588,588,205,92,-44,-140,-240,-240]]
; aperture stop radius
rap = 58./2.
module=['A','B']
caldbdir=getenv('CALDB')+'/'
alignfile=caldbdir+'data/nustar/fpm/bcf/align/'
if not file_test(alignfile+'nuCalign20100101v001.fits') then $
      stop,'GRXE2DET: Fatal error, caldb changed since code last updated.'+$
            '  Grab latest code from svn repository; if code up to date, '+$
            'then email daniel.r.wik@nasa.gov this message or edit code to '+$
            'point to the correct alignment file.'
version=2
while file_test(alignfile+'nuCalign20100101v'+string(version,format='(I03)')+$
      '.fits') do version++
alignfile=alignfile+'nuCalign20100101v'+string(version-1,format='(I03)')+'.fits'
; get relative z-axis positions from alignment file
readCALDB,alignfile,module=module[imod]
if imod eq 0 then begin
    apcent=caldbinstrument.v_fb_asa[0:1]
    apz=caldbinstrument.v_fb_asa[2]-caldbmodule.v_fb_fpm[2]
    omz=caldbinstrument.v_ob_oma[2]
endif else begin
    apcent=caldbinstrument.v_fb_asb[0:1]
    apz=caldbinstrument.v_fb_asb[2]-caldbmodule.v_fb_fpm[2]
    omz=caldbinstrument.v_ob_omb[2]
    omoutline[*,0]-=caldbinstrument.v_ob_omb[0]-caldbinstrument.v_ob_oma[0]
endelse

; mask out part of sky blocked by optics bench
obscure=polyfillv(omoutline[*,0]+dim/2,omoutline[*,1]+dim/2,dim,dim)
optimblur=fltarr(dim,dim)
optimblur[*,*]=1
optimblur[obscure]=0
optimblur=gblur(optimblur,blursig,3)

; polygon outline of aperture stop
ang=[findgen(72)*5.*!pi/180.,0.]
apxy=fltarr(n_elements(ang),2)
for k=0,n_elements(ang)-1 do apxy[k,*]=rap*[cos(ang[k]),sin(ang[k])]+apcent
vfb=fltarr(n_elements(ang),2)

; set up grid on detector plane
submm=0.12096*nbinpix
dx=findgen(70/submm)*submm-35.
dy=findgen(70/submm)*submm-35.
imdet=fltarr(n_elements(dx),n_elements(dy))
optimall=fltarr(dim,dim)
; loop over detector plane sample points (if nbinpix=4, then the straylight
;  contribution is found for every 4th pixel)
;  NOTE: this does NOT directly correspond to det1 pixels, but instead positions
;  centered around a 7x7 cm box where the detectors are supposed to be -- there
;  is a slight shift between the center of the detectors and this position
for i=0,n_elements(dx)-1 do for j=0,n_elements(dy)-1 do begin
    vfb2=qtvrotation([dx[i],dy[j],0],caldbmodule.q_fb_fpm,/inv)+caldbmodule.v_fb_fpm
    vfb[*,0]=vfb2[0]
    vfb[*,1]=vfb2[1]
    vis=(apxy-vfb)*(caldbinstrument.v_fb_ob[2]+2*omz-caldbmodule.v_fb_fpm[2])/apz
; ORIGINAL Z-DISTANCE, EQUALS FOCAL LENGTH BUT THIS IS AT THE CENTER OF THE
;      OPTICS, NOT THE PART CLOSEST TO THE FOCAL PLANE WHICH IS WHAT WE WANT
;    vis=(apxy-vfb)*(caldbinstrument.v_fb_ob[2]+omz-caldbmodule.v_fb_fpm[2])/apz

; find part of sky visible through the aperture stop
    optim=fltarr(dim,dim)
    visible=polyfillv(vis[*,0]+dim/2,vis[*,1]+dim/2,dim,dim)
    optim[visible]=1
;    optim[obscure]=0

; remove any part of aperture that is blocked by the optics modules
    optim=optim*optimblur*grxe_image
    optimall+=optim
    tvscl,congrid(optimall,500,500)
    imdet[i,j]=total(optim)
    if j eq 0 and i mod 2 eq 0 then $
          print,float(i)/n_elements(dx)*100.,'%..',format='($,F5.1,A)'
endfor
print
print

; interpolates so the output image has pixels equivalent in size to det1 pixels
imdet=congrid(imdet,n_elements(imdet[*,0])*nbinpix,$
      n_elements(imdet[0,*])*nbinpix,/interp)

; sky coverage output image
fits_write,dir+'skygrxe'+module[imod]+'_det1.img',optimall,hdr_optim
; detector plane image, pixel size same as det1 coords, but over larger area
fits_write,dir+'grxe'+module[imod]+'_det1.img',imdet


endif else begin

print,'Pre-existing GRXE image found, not recreated.'
print,'  To remake the images, set the clobber keyword.'

endelse

end
