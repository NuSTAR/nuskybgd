pro fakinstr,indir,obsid,ab,paramfile,bgddir
      

auxildir=getenv('NUSKYBGD_AUXIL')+'/'
dir=indir
if strmid(dir,strlen(dir)-1) ne '/' then dir=dir+'/'
if ab eq 'A' then iab=0 else if ab eq 'B' then iab=1 else stop,'ab defined wrong'
ctstat='n'

if not size(bgddir,/type) then bgddir=''
cldir=dir+obsid+'/event_cl/'+bgddir+'/'
if not size(paramfile,/type) then $
      paramfile=cldir+'bgdfitparams'+ab+'.dat'
print,paramfile
if not file_test(paramfile) then stop,'FAKINSTR: Cannot find parameter file.'

;readcol,auxildir+'ratios_lineE.dat',eline,width,/silent
;readcol,auxildir+'ratios_lineE.dat',blah,index1,ebreak,index2,$
;      format='(A,F,F,F)',/silent
;readcol,auxildir+'ratios'+ab+'.dat',f0,f1,f2,f3,/silent
readcol,auxildir+'ratios'+ab+'.dat',eline,width,f0,f1,f2,f3,/silent
readcol,auxildir+'ratios'+ab+'.dat',index1,index2,b0,b1,b2,b3,ebreak,/silent
neut=[eline[n_elements(f0)-1],width[n_elements(f0)-1]]
eline=eline[0:n_elements(width)-3]
width=width[0:n_elements(width)-3]

caldb=getenv('CALDB')+'/'
pt=loadnuabs(0)
czt=loadnuabs(1)
readcol,paramfile,p0,p1,p2,p3,/silent,skipline=3

for idet=0,3 do begin

if idet eq 0 then p=p0 else if idet eq 1 then p=p1 else if idet eq 2 then p=p2 $
      else p=p3

spt=pt[idet,iab]
sczt=czt[idet,iab]

rmfname=caldb+'data/nustar/fpm/cpf/rmf/nu'+ab+'cutdet'+str(idet)+'_20100101v001.rmf'
fakname='instrbgdfit'+ab+'_det'+str(idet)+'.pha'

openw,lun,'temp.xcm',/get_lun
printf,lun,'lmod nuabs'
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
    printf,lun,str(p[i])
endfor
printf,lun,str(index1[0])+' -1'
printf,lun,str(index2[0])+' -1'
printf,lun,str(ebreak[0])+' -1'
printf,lun,str(p[n_elements(eline)])
;printf,lun,str(neut[0])
;printf,lun,str(neut[1])
;printf,lun,str(p[n_elements(eline)+1])
spawn,'rm -f '+cldir+fakname
printf,lun,'fakeit none & '+rmfname+' &  & '+ctstat+' &  & '+$
      cldir+fakname+' & 1e9'
printf,lun,'exit'
free_lun,lun

spawn,'xspec - temp.xcm'
spawn,'rm -f temp.xcm'

endfor

end
