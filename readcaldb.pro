PRO readCALDB, file, module=module

  
  if file eq 'none' then $
    file = '/Users/kristin/data/IOCtoolsKKM/IOCtoolsKKM/align_caldb/nuCalign20100101v003_4.fits'

  common caldb, caldbmodule, caldbinstrument, caldbmetrology, caldbOA

  caldbinstrument = mrdfits(file,1,hh,/silent)
  caldbmetrology = mrdfits(file,2,hh,/silent)
  caldbOA = mrdfits(file,3,hh,/silent)

  if keyword_set(module) then begin 
 
    if module eq 'A' then begin
      v_fb_fpm = caldbinstrument.V_FB_FPMA
      q_fb_fpm = caldbinstrument.Q_FB_FPMA
      v_ob_om = caldbinstrument.V_OB_OMA
      q_ob_om = caldbinstrument.Q_OB_OMA
      v_fpm_det1 = caldbinstrument.V_FPMA_DET1
      q_fpm_det1 = caldbinstrument.Q_FPMA_DET1
      v_det2_ob = caldbinstrument.V_DET2A_OB 
      q_det2_ob = caldbinstrument.Q_DET2A_OB 
      om_point = caldbOA.OMA_POINT
    end
    if module eq 'B' then begin
      v_fb_fpm = caldbinstrument.V_FB_FPMB
      q_fb_fpm = caldbinstrument.Q_FB_FPMB
      v_ob_om = caldbinstrument.V_OB_OMB
      q_ob_om = caldbinstrument.Q_OB_OMB
      v_fpm_det1 = caldbinstrument.V_FPMB_DET1
      q_fpm_det1 = caldbinstrument.Q_FPMB_DET1
      v_det2_ob = caldbinstrument.V_DET2B_OB 
      q_det2_ob = caldbinstrument.Q_DET2B_OB 
      om_point = caldbOA.OMB_POINT
    end


    caldbmodule = {$
      V_IN_SC:caldbinstrument.V_IN_SC,$
      Q_IN_SC:caldbinstrument.Q_IN_SC,$
      V_SC_FB:caldbinstrument.V_SC_FB,$
      Q_SC_FB:caldbinstrument.Q_SC_FB,$
      V_FB_FPM:V_FB_FPM,$
      Q_FB_FPM:Q_FB_FPM,$
      V_FB_OB:caldbinstrument.V_FB_OB,$ 
      Q_FB_OB:caldbinstrument.Q_FB_OB,$ 
      V_OB_OM:V_OB_OM,$
      Q_OB_OM:Q_OB_OM,$
      V_OB_ST:caldbinstrument.V_OB_ST,$ 
      Q_OB_ST:caldbinstrument.Q_OB_ST,$ 
      V_FPM_DET1:V_FPM_DET1,$
      Q_FPM_DET1:Q_FPM_DET1,$
      V_DET2_OB:V_DET2_OB,$ 
      Q_DET2_OB:Q_DET2_OB,$ 
      OM_POINT:om_point $
    } 

  endif else caldbmodule = 0 



END

