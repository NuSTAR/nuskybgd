function getcaldbfile,type,ab,par1,par2,par3,par4,par5

indx=mrdfits(getenv('CALDB')+'/data/nustar/fpm/caldb.indx',1,/silent)
if ab ne 'A' and ab ne 'B' then inst='FPM' else inst='FPM'+ab
if type eq 'psf' then begin
    cnam='GRPPSF'
    ii=where(str(indx.cal_cnam) eq cnam and indx.cal_qual eq 0 and $
          str(indx.instrume) eq inst)
    psf=mrdfits(getenv('CALDB')+'/'+str(indx[ii].cal_dir)+'/'+$
          str(indx[ii].cal_file),1,/silent)
    jj=where(par1 ge psf.energ_lo and par1 lt psf.energ_hi)
    if jj[0] eq -1 then stop,'GETCALDBFILE: error retrieving psf file'
    returnstr=getenv('CALDB')+'/'+str(indx[ii].cal_dir)+'/'+str(psf[jj]._2dpsffile)
endif else if type eq 'arf' then begin
    cnam='SPECRESP'
    ii=where(str(indx.cal_cnam) eq cnam and indx.cal_qual eq 0 and $
          str(indx.instrume) eq inst)
    returnstr=getenv('CALDB')+'/'+str(indx[ii].cal_dir)+'/'+str(indx[ii].cal_file)
endif else if type eq 'vign' then begin
    cnam='TVIGNET'
    ii=where(str(indx.cal_cnam) eq cnam and indx.cal_qual eq 0 and $
          str(indx.instrume) eq inst)
    returnstr=getenv('CALDB')+'/'+str(indx[ii].cal_dir)+'/'+str(indx[ii].cal_file)
endif else if type eq 'nuabs' or type eq 'detabs' then begin
    cnam='DETABS'
    if size(par1,/type) ne 0 then $
          ii=where(str(indx.cal_cnam) eq cnam and indx.cal_qual eq 0 and $
            str(indx.instrume) eq inst and str(indx.detnam) eq 'DET'+str(par1)) $
      else ii=where(str(indx.cal_cnam) eq cnam and indx.cal_qual eq 0 and $
            str(indx.instrume) eq inst)
    returnstr=getenv('CALDB')+'/'+str(indx[ii].cal_dir)+'/'+str(indx[ii].cal_file)
endif else if type eq 'instrmap' then begin
    cnam='INSTRMAP'
    ii=where(str(indx.cal_cnam) eq cnam and indx.cal_qual eq 0 and $
          str(indx.instrume) eq inst)
    returnstr=getenv('CALDB')+'/'+str(indx[ii].cal_dir)+'/'+str(indx[ii].cal_file)
endif else if type eq 'pixpos' then begin
    cnam='PIXPOS'
    ii=where(str(indx.cal_cnam) eq cnam and indx.cal_qual eq 0 and $
          str(indx.instrume) eq inst)
    returnstr=getenv('CALDB')+'/'+str(indx[ii].cal_dir)+'/'+str(indx[ii].cal_file)
endif else if type eq 'rmf' then begin
;    cnam='MATRIX'
    cnam='GRPRMF'
    ii=where(str(indx.cal_cnam) eq cnam and indx.cal_qual eq 0 and $
          str(indx.instrume) eq inst and str(indx.detnam) eq 'DET'+str(par1) and $
          strmid(indx.cal_cbd,9,3) eq 'NOM')
    rmf=mrdfits(getenv('CALDB')+'/'+str(indx[ii].cal_dir)+'/'+$
          str(indx[ii].cal_file),par1+1,/silent)
    returnstr=getenv('CALDB')+'/'+str(indx[ii].cal_dir)+'/'+str(rmf[0].rmffile)
endif else if type eq 'apstop' then begin
    cnam='APERTURE'
    ii=where(str(indx.cal_cnam) eq cnam and indx.cal_qual eq 0 and $
          str(indx.instrume) eq inst)
    returnstr=getenv('CALDB')+'/'+str(indx[ii].cal_dir)+'/'+str(indx[ii].cal_file)
endif else stop,'GETCALDBFILE: requested item unavailable'


return,returnstr

end
