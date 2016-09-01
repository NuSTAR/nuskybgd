function gblur,inputim,sig,extent

gdim=ceil(sig*extent)
xim=fltarr(2*gdim+1,2*gdim+1)
yim=fltarr(2*gdim+1,2*gdim+1)
for i=0,2*gdim do begin
    xim[*,i]=findgen(2*gdim+1)-gdim
    yim[i,*]=findgen(2*gdim+1)-gdim
endfor
rim2=xim*xim+yim*yim
gim=exp(-rim2/2./sig^2)
gim=gim/total(gim)

xdim=n_elements(inputim[*,0])
ydim=n_elements(inputim[0,*])
newim=fltarr(xdim,ydim)
for i=gdim,xdim-1-gdim do for j=gdim,ydim-1-gdim do $
      newim[i-gdim:i+gdim,j-gdim:j+gdim]+=inputim[i,j]*gim
newim[0:2*gdim,*]=inputim[0:2*gdim,*]
newim[*,0:2*gdim]=inputim[*,0:2*gdim]
newim[xdim-2*gdim:xdim-1,*]=inputim[xdim-2*gdim:xdim-1,*]
newim[*,ydim-2*gdim:ydim-1]=inputim[*,ydim-2*gdim:ydim-1]

return,newim

end
