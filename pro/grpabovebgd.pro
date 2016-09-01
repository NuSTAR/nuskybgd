pro grpabovebgd,inspec,bgdspec,outspec,ncts,chan=chan,backcor=backcor

if not keyword_set(chan) then chan=0

spec=mrdfits(inspec,1,head,/silent)
bgd=mrdfits(bgdspec,1,bhead,/silent)

src_exp = sxpar(head,'EXPOSURE')
src_backscal = sxpar(head, 'BACKSCAL')

bgd_exp = sxpar(bhead, 'EXPOSURE')
bgd_backscal = sxpar(bhead, 'BACKSCAL')

bspec = float(bgd.counts) * ( float(src_backscal) * float(src_exp)) /  ( float(bgd_backscal) * float(bgd_exp)) 

;expfactor=sxpar(head,'EXPOSURE')/sxpar(bhead,'EXPOSURE')
;bspec=float(bgd.counts)*expfactor
;bspec = float(bgd.counts) 

if keyword_set(backcor) then bspec=(float(bgd.counts)+float(bgd.counts)*backcor)*expfactor

cts=0
bcts=0
firstval=1
for c=0,n_elements(bspec)-1 do begin
    cts+=spec[c].counts
    bcts+=bspec[c]
    if cts-bcts ge ncts or c eq chan-1 then begin
        cts=0
        bcts=0
        if firstval then spec[c].grouping = 1 else spec[c].grouping = -1
        firstval=1
    endif else begin
        if firstval eq 1 then spec[c].grouping = 1 else spec[c].grouping = -1
        firstval=0
    endelse
endfor

mwrfits,spec,outspec,head,/create

end
