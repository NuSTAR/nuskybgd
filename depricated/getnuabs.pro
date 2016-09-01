function getnuabs,dir,obsid,srcreg,bgddir,ab

pt=loadnuabs(0)
czt=loadnuabs(1)

cldir=dir+'/'+obsid+'/event_cl/'
mask=reg2mask(cldir+bgddir+'/det0'+ab+'im.fits',cldir+srcreg)
if ab eq 'A' then iab=0 else iab=1

detfrac=fltarr(4)
dettot=fltarr(4)
for i=0,3 do begin
    fits_read,cldir+bgddir+'/det'+str(i)+ab+'im.fits',detim
    detfrac[i]=total(detim*mask)
    dettot[i]=total(detim)
endfor
detwt=detfrac/total(detfrac)

spt=total(pt[*,iab]*detwt)
sczt=total(czt[*,iab]*detwt)

print,'Pt =  '+str(spt)
print,'CZT = '+str(sczt)

return,[spt,sczt]

end
