; NUSKYBGD_IMREFSPEC
;
; Creates reference spectra of the 3 soft bgd components, used by nuskybgd_image, 
; to produce counts images for arbitrary energy bands.  Only needs to be run once
; unless the nuabs parameters change in the auxil directory.


pro nuskybgd_imrefspec

caldbdir=getenv('CALDB')+'/'
auxildir=getenv('NUSKYBGD_AUXIL')+'/'
pt=loadnuabs(0)
czt=loadnuabs(1)

; CXB params, from Boldt '87
pl=str(1.29)
norm=str(0.002353)
ecut=str(41.13)

; Soft neutron component photon index
neutpl=str(4.8)

; GRXE params
glenergy=str(6.7)
wdmass=str(0.6)

; Simulated exposure time (assumed elsewhere, do not change!)
texp=str(1e9)

openw,lun,'temp.xcm',/get_lun
printf,lun,'lmod nuabs'
printf,lun,'set pt "'+str(pt[0]),format='($,A)'
for i=1,7 do printf,lun,' '+str(pt[i]),format='($,A)'
printf,lun,'"'
printf,lun,'set czt "'+str(czt[0]),format='($,A)'
for i=1,7 do printf,lun,' '+str(czt[i]),format='($,A)'
printf,lun,'"'
printf,lun,'foreach {i} {0 1 2 3} {'
printf,lun,'  foreach {j} {A B} {'
printf,lun,'    if {$j == "B"} {set off 4} else {set off 0}'
printf,lun,'    model nuabs*po*highecut & [lindex $pt [expr $off+$i]] & '+$
      '[lindex $czt [expr $off+$i]] & 0. & 0.9 & '+pl+' & '+norm+' & 1e-4 & '+ecut
printf,lun,'    fakeit none & '+caldbdir+'data/nustar/fpm/cpf/rmf/'+$
      'nu${j}cutdet${i}_20100101v001.rmf & '+auxildir+'be.arf & no & & '+$
      auxildir+'aperbgd${j}_det${i}.pha & '+texp+', 0.'
printf,lun,'    fakeit none & '+caldbdir+'data/nustar/fpm/cpf/rmf/'+$
      'nu${j}cutdet${i}_20100101v001.rmf & '+auxildir+'fcxb${j}.arf & no & & '+$
      auxildir+'fcxbbgd${j}_det${i}.pha & '+texp+', 0.'
;printf,lun,'    fakeit none & '+caldbdir+'data/nustar/fpm/cpf/rmf/'+$
;      'nu${j}cutdet${i}_20100101v001.rmf & '+caldbdir+'data/nustar/fpm/bcf/arf/'+$
;      'nu${j}20100101v004.arf & no & & '+auxildir+$
;      'fcxbbgd${j}_det${i}.pha & '+texp+', 0.'
;;;printf,lun,'    model nuabs*po & [lindex $pt [expr $off+$i]] & [lindex $czt '+$
;;;      '[expr $off+$i]] & 0. & 0.9 & '+neutpl+' & 1.0'
;;;printf,lun,'    fakeit none & '+caldbdir+'data/nustar/fpm/cpf/rmf/'+$
;;;      'nu${j}cutdet${i}_20100101v001.rmf & '+auxildir+'be.arf & no & & '+$
;;;      auxildir+'neutbgd${j}_det${i}.pha & '+texp+', 0.'
; GRXE stuff -- causing issues for some reason, unclear why -- comment until needed
;printf,lun,'    model nuabs*gauss & [lindex $pt [expr $off+$i]] & '+$
;      '[lindex $czt [expr $off+$i]] & 0. & 0.9 & '+glenergy+' & 0.0 & 1.0'
;printf,lun,'    fakeit none & '+caldbdir+'data/nustar/fpm/cpf/rmf/'+$
;      'nu${j}cutdet${i}_20100101v001.rmf & '+auxildir+'be.arf & no & & '+$
;      auxildir+'grxe1bgd${j}_det${i}.pha & '+texp+', 0.'
;printf,lun,'    model nuabs*atable{'+auxildir+'polarmodel.fits} & '+$
;      '[lindex $pt [expr $off+$i]] & '+$
;      '[lindex $czt [expr $off+$i]] & 0. & 0.9 & '+wdmass+' & 1.0'
;printf,lun,'    fakeit none & '+caldbdir+'data/nustar/fpm/cpf/rmf/'+$
;      'nu${j}cutdet${i}_20100101v001.rmf & '+auxildir+'be.arf & no & & '+$
;      auxildir+'grxe2bgd${j}_det${i}.pha & '+texp+', 0.'
printf,lun,'  }'
printf,lun,'}'
printf,lun,'exit'
free_lun,lun

spawn,'xspec - temp.xcm'
spawn,'rm -f temp.xcm'

end
