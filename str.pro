function str,s,l=l,_extra=extra

snew = strcompress(string(s,_extra=extra),/remove_all)

if keyword_set(l) then begin
    len = strlen(snew)
    while len lt l do begin
	snew = ' '+snew
	len = strlen(snew)
    endwhile
    return,snew
endif else return,snew

end
