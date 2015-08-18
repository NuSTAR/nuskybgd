function loadnuabs,ptorczt

auxildir=getenv('NUSKYBGD_AUXIL')+'/'
readcol,auxildir+'nuabs.dat',pta,czta,ptb,cztb,/silent

if ptorczt eq 0 then return,[[pta],[ptb]] else if ptorczt eq 1 then $
      return,[[czta],[cztb]] else stop,'LOADNUABS: incorrect argument'

end
