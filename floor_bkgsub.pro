pro floor_bkgsub

cd, '/Users/bwgref/science/casA/image_analysis/split_files/bgd_sub'
f = file_search('nu*bsub.fits')
expdir = '/Users/bwgref/science/casA/image_analysis/expo_maps/'

for i = 0, n_elements(f) -1 do begin
   print, f[i]
   image_file = f[i]
   
   prefix = (strsplit(f[i], '_', /extract))[0]
   
   exp_file = expdir+prefix+'_ex_5.img'



   cts_image = mrdfits(image_file, 0, header)
   exp_image = mrdfits(exp_file, 0)

   not_finite = where(~finite(cts_image / exp_image))

   cts_image[NOT_FINITE] = 0

;   plot, cts_image / exp_image
   
;   cgimage, cts_image / exp_image, /scale, /keep_aspect
;   stop

   
   low_exposure = where(exp_image lt 0.2 * mean(exp_image[where(exp_image ne 0)]))   
;   plot, exp_image
;   oplot, low_exposure * [0, n_elements(exp_image)]

   cts_image[low_exposure] = 0.

   filename = file_basename(image_file, '.fits')+'_floor.fits'
   fi = file_info(filename)
;   if fi.exists then continue

   print, filename
   
   mwrfits, cts_image, filename, header, /create
endfor

end


