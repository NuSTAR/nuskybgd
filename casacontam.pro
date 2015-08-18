pro casacontam,dir,obsid,reg,pl,kt
	
if n_elements(reg) ne n_elements(pl) and n_elements(pl) ne 1 then $
      stop,'CASACONTAM: input arrays screwy silly person'

caldbdir=getenv('CALDB')+'/'
ab=['A','B']


openw,lun,dir+'/'+obsid+'/event_cl/contam.xcm',/get_lun
for iab=0,0 do for ir=0,n_elements(reg)-1 do begin
     printf,lun,'resp 1:'+str(n_elements(reg)/2*iab+ir+1)+' '+$
           dir+'/'+obsid+'/event_cl/bgdspec/'+reg[ir]+'_sr.rmf'
     printf,lun,'arf 1:'+str(n_elements(reg)/2*iab+ir+1)+' '+caldbdir+$
           'data/nustar/fpm/bcf/arf/nu'+ab[iab]+'20100101v004.arf'
endfor
printf,lun,'model nuabs*(gauss+po+mekal)'
;for iab=0,0 do for ir=0,n_elements(reg)-1 do begin
for iab=0,1 do for ir=0,n_elements(reg)/2-1 do begin
    nuabs=getnuabs(dir,obsid,reg[ir]+'.reg','bgd',ab[iab])
    printf,lun,str(nuabs[0])+' -1'
    printf,lun,str(nuabs[1])+' -1'
    printf,lun,'0. -1'
    printf,lun,'0.9 -1'
    printf,lun,'6.64 -1'
    printf,lun,'0.1 -0.1'
    printf,lun,'1e-5 1e-6'
    printf,lun,str(pl)+' -1'
    printf,lun,'3e-4 1e-5'
    printf,lun,str(kt)+' -1'
    printf,lun,'1e-6 -1'
    printf,lun,'0. -1'
    printf,lun,'0. -1'
    printf,lun,'0'
    printf,lun,'1e-2 1e-3'
; code from thermal model, good reference if you fit A and B together
;    if iab eq 0 then printf,lun,str(kt[ir])+' -1' else printf,lun,'='+str(5+ir*10)
;    printf,lun,'1e-6 -1'
;    if iab eq 0 then printf,lun,'0.2 0.01' else printf,lun,'='+str(7+ir*10)
;    printf,lun,str(z)+' -0.001'
;    printf,lun,'0'
;    if iab eq 0 then printf,lun,'1e-3 0.01' else printf,lun,'='+str(10+ir*10)
endfor
free_lun,lun

end
