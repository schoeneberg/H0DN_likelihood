Pro print_residuals, details, prefix=prefix
;
; Structured print of residual files to separate files.  Files are
; named for the kind of residual in each, with an optional prefix to
; distinguish the run.
;

if (keyword_set(prefix) eq 0) then prefix = './'

; Identify hosts

mas                = details.mas_details.mas
ddh                = details.host_data_details
hosts              = ddh.hosts
n2                 = ddh.n2
host2              = ddh.host2
mu_host2_error     = ddh.mu_host2_error
mas_index          = ddh.mas_index
mas_host_error     = ddh.mas_host_error
anchor_dist_error  = ddh.anchor_dist_error
res                = details.residuals
c = !const.c/1.e3


if (details.do_cull eq 0) then begin
   ; Match host2 with the hosts file, and remove if 
   ; the host is not used for any calibrators
   accept_host = bytarr (n2)
   for k = 0, n_elements(hosts) - 1 do begin
      wh = where (ddh.host2 eq hosts[k], nwh)
      if (details.host_used[k] gt 0 and nwh gt 0) then accept_host[wh] = 1B
   endfor
endif else accept_host = replicate (1B, n2)
   
; For each type of residual (each element in mas) output the residuals
; to an appropriately named file

mas_type = !null
mas_offset = !null
mas_error = !null
mas_number = !null

for k = 0, n_elements(mas)-1 do begin
   mas_name = mas[k]
   mas_pieces = strsplit (mas_name, '&', /extract)
   mas_pieces = strtrim(mas_pieces, 2)
   npieces = n_elements(mas_pieces)
   filename = prefix + mas_pieces[0]
   for j = 1, npieces-1 do filename = filename + '__' + mas_pieces[j]
   filename = filename + '.dat'
   ; Determine which residuals belong to this MAS
   whsub = where (mas_index eq k and accept_host, nsub)
   if (nsub le 0) then goto, printed
   subres = res [whsub]   ; This assumes that the host equations start at 0
   subhost = host2[whsub]
   suberror = mu_host2_error[whsub]  ; The error for this host measurement
   referror = mas_host_error[whsub[0]] ; The common mas reference error
   disterror = anchor_dist_error[whsub[0]] ; The common anchor distance error
   ; All information available.  Print.
   w = 1./suberror^2
   mean_offset = total (w*subres) / total(w)
   mean_error = sqrt (1/total(w))
   mean_error = sqrt (mean_error^2 + referror^2 + disterror^2) 
   mas_type = [mas_type, mas_name]
   mas_number = [mas_number, nsub]
   mas_offset = [mas_offset, mean_offset]
   mas_error = [mas_error, mean_error]
   openw, lun, filename, /get_lun
   printf, lun, '# Residuals for '+mas_name
   printf, lun, '# Mean, uncertainty: ', mean_offset, mean_error, $
           format = '(a, 2f10.3)'
   for j = 0, nsub-1 do begin
      printf, lun, subhost[j], subres[j], suberror[j], referror, disterror, $
              format = '(a20, 4x, 4f10.4)'
   endfor
   close, lun
   free_lun, lun
   printed:
endfor

if (n_elements(mas_type) gt 0) then begin
   openw, lun, prefix+'host_summary.dat', /get_lun
   for k = 0, n_elements(mas_type)-1 do printf, lun, mas_offset[k], $
      mas_error[k], mas_number[k], mas_type[k], format='(2f10.4, i4, 2x, a)'
   close, lun
   free_lun, lun
endif



; Now go through the other categories of residuals.

; Coma

n_coma              = details.coma_details.n_coma             
if (n_coma gt 0) then begin
   coma_sn_name       = details.coma_details.coma_sn_name      
   ieq_coma_start     = details.coma_details.ieq_coma_start    
   n_coma              = details.coma_details.n_coma             
   coma_sn_err        = details.coma_details.coma_sn_err
   mu_coma_value      = details.coma_details.mu_coma_value
   mu_coma_error      = details.coma_details.mu_coma_error
                        ; uncertainty in Coma distance
   filename = prefix + 'Coma.dat'
   wh = ieq_coma_start + indgen(n_coma) ; Coma residuals
   ; these are supernova distance moduli minus the Coma
   ; distance modulus
   subres = res [wh]
   source_names = coma_sn_name
   suberr = coma_sn_err         ; Source-by-source uncertainty for Coma SNe
   openw, lun, filename, /get_lun
   for j = 0, n_coma-1 do printf, lun, source_names[j], subres[j], suberr[j], $
                                  mu_coma_error, format = '(a20, 4x, 3f10.4)'
   close, lun
   free_lun, lun
endif

; Residuals in the Hubble flow.

; These can be shown as the residual between the individual distance
; estimates based on the calibration of each method and the distance
; based on the corrected redshift.  That information for Hubble Flow
; objects is saved by compute_alpha starting in version 3.15.
; The information for MM is as used in make_residual_plots_v3, since
; we do not yet use compute_alpha for MM (would need to be run for
; angular diameter distance)

; mm
n_mm              = details.mm_details.n_mm
if (n_mm gt 0) then begin
   mm_name        = details.mm_details.mm_name
   ieq_mm_start   = details.mm_details.ieq_mm_start
   mm_logh_error  = details.mm_details.mm_logh_error
   mm_v_corr      = details.mm_details.mm_v_corr
   vcorr_mm       = details.mm_details.vcorr_mm ; correction code
   wh = ieq_mm_start + indgen(n_mm) ; SBF residuals
   subres = res[wh]
   suberr = mm_logh_error
   filename = prefix + 'mm_hf.dat'
   openw, lun, filename, /get_lun
   printf, lun, '# Megamaser residuals in log(H0) vs redshift'
   printf, lun, ' # Mean, dispersion = ', avg(subres), sigma(subres), format='(a, 2f12.5)'
   zhd = sqrt((1+mm_v_corr/c)/(1-mm_v_corr/c))-1
   for j = 0, n_mm-1 do printf, lun, mm_name[j], zhd[j], subres[j], suberr[j], format='(a20, 2x, 3f12.4)'
   close, lun
   free_lun, lun
endif

; SN1A in Hubble flow
; Information from alphadet_sn1a
det = details.sn1a_details
if (det.do_sn1a) then if (details.sn1a_details.read_sn1a_hf) then begin
   alphadet = det.alphadet_sn1a
   n_sn1a_hf = alphadet.number
   zhd = alphadet.zhd
   alphavec = alphadet.alphavec ; individual values of alpha
   alphavec_err = sqrt(alphadet.alphavar) ; diagonal element of covariance
   alpha_value = details.sn1a_details.alpha_sn1a_value
   ieq = details.sn1a_details.ieq_sn1a_end ; equation for alpha_sn1a
   ; The equation involving alpha_sn1a has the form:
   ; log(H0) - 0.2*M0 = alpha_sn1a + 5 == yval[ieq]
   ; The residual is defined as residual = yval[ieq] - log(H0) + 0.2*M0
   ; So the inferred value for alphavec is 
   res_value = res[ieq]
   alphaoff = alphavec - alpha_value + res_value
   filename = prefix + 'sn1a_hf.dat'
   openw, lun, filename, /get_lun
   printf, lun, '# SN Ia  residuals in log(H0) vs redshift'
   printf, lun, '# Mean, dispersion = ', avg(alphaoff), $
           sigma(alphaoff), format='(a, 2f12.5)'
   name = det.name_sn1a_hf
   for j = 0, n_sn1a_hf-1 do printf, lun, name[j], zhd[j], alphaoff[j], $
                                     alphavec_err[j], format='(a20, 2x, 3f12.4)'
   close, lun
   free_lun, lun
endif

; SN2 in Hubble flow
; Information from alphadet_sn2
det = details.sn2_details
if (det.do_sn2) then begin
   alphadet = det.alphadet_sn2
   n_sn2_hf = alphadet.number
   zhd = alphadet.zhd
   alphavec = alphadet.alphavec ; individual values of alpha
   alphavec_err = sqrt(alphadet.alphavar) ; diagonal element of covariance
   alpha_value = det.alpha_sn2_value
   ieq = det.ieq_sn2_end ; equation for alpha_sn1a
   ; The equation involving alpha_sn1a has the form:
   ; log(H0) - 0.2*M0 = alpha_sn1a + 5 == yval[ieq]
   ; The residual is defined as residual = yval[ieq] - log(H0) + 0.2*M0
   ; So the inferred value for alphavec is 
   res_value = res[ieq]
   alphaoff = alphavec - alpha_value + res_value
   filename = prefix + 'sn2_hf.dat'
   openw, lun, filename, /get_lun
   printf, lun, '# SN II  residuals in log(H0) vs redshift'
   printf, lun, '# Mean, dispersion = ', avg(alphaoff), sigma(alphaoff), format='(a, 2f12.5)'
   name = det.name_sn2_hf
   for j = 0, n_sn2_hf-1 do printf, lun, name[j], zhd[j], alphaoff[j], alphavec_err[j], format='(a20, 2x, 3f12.4)'
   close, lun
   free_lun, lun
endif

; sbf.  Redone to account for change of the SBF solution.  Since
; Version 3.4, Hubble flow SBF objects are in details.sbf_details.alphadet_sbf
; Information from alphadet_sbf.

det = details.sbf_details
if (det.do_sbf) then begin
   alphadet = det.alphadet_sbf
   n_sbf_hf = alphadet.number
   zhd = alphadet.zhd
   alphavec = alphadet.alphavec ; individual values of alpha
   alphavec_err = sqrt(alphadet.alphavar) ; diagonal element of covariance
   alpha_value = det.alpha_sbf_value
   ieq = det.ieq_sbf_end ; equation for alpha_sbf
   ; The equation involving alpha_sbf has the form:
   ; log(H0) - 0.2*M0 = alpha_sbf + 5 == yval[ieq]
   ; The residual is defined as residual = yval[ieq] - log(H0) + 0.2*M0
   ; So the inferred value for alphavec is 
   res_value = res[ieq]
   alphaoff = alphavec - alpha_value + res_value
   filename = prefix + 'sbf_hf.dat'
   openw, lun, filename, /get_lun
   printf, lun, '# SBF  residuals in log(H0) vs redshift'
   printf, lun, '# Mean, dispersion = ', avg(alphaoff), sigma(alphaoff), $
           format='(a, 2f12.5)'
   name = det.name_sbf_hf
   for j = 0, n_sbf_hf-1 do printf, lun, name[j], zhd[j], alphaoff[j], $
                                    alphavec_err[j], format='(a20, 2x, 3f12.4)'
   close, lun
   free_lun, lun
endif

; TF in Hubble flow
det = details.tf_details
if (det.do_tf) then begin
   alphadet = det.alphadet_tf
   n_tf_hf = alphadet.number
   zhd = alphadet.zhd
   alphavec = alphadet.alphavec ; individual values of alpha
   alphavec_err = sqrt(alphadet.alphavar) ; diagonal element of covariance
   alpha_value = det.alpha_tf_value
   ieq = det.ieq_tf_end ; equation for alpha_tf
   ; The equation involving alpha_sn1a has the form:
   ; log(H0) - 0.2*M0 = alpha_sn1a + 5 == yval[ieq]
   ; The residual is defined as residual = yval[ieq] - log(H0) + 0.2*M0
   ; So the inferred value for alphavec is 
   res_value = res[ieq]
   alphaoff = alphavec - alpha_value + res_value
   filename = prefix + 'tf_hf.dat'
   openw, lun, filename, /get_lun
   printf, lun, '# TF  residuals in log(H0) vs redshift'
   printf, lun, '# Mean, dispersion = ', avg(alphaoff), sigma(alphaoff), format='(a, 2f12.5)'
   name = det.name_tf_hf
   for j = 0, n_tf_hf-1 do printf, lun, name[j], zhd[j], alphaoff[j], alphavec_err[j], format='(a20, 2x, 3f12.4)'
   close, lun
   free_lun, lun
endif

; epm
; ; n_epm            = details.epm_details.n_epm
; ; if (n_epm gt 0) then begin
; ;    epm_name        = details.epm_details.epm_name
; ;    ieq_epm_start   = details.epm_details.ieq_epm_start
; ;    epm_logh_error  = details.epm_details.epm_logh_error
; ;    epm_v_corr      = details.epm_details.epm_v_corr
; ;    vcorr_epm       = details.epm_details.vcorr_epm ; correction code
; ;          wh = ieq_epm_start + indgen(n_epm) ; EPM residuals
; ;          subres = residuals[wh]
; ;          suberr = epm_logh_error
; ;          ; includes single-galaxy distance and velocity correction errors
; ;          ; removes the covariant term in the distance calibration
; ;          ; cal_sbf_error (the distance calibration error) has also 
; ;          ; been added to the distance error from Jensen+ 2021
; ;          vrange = max(epm_v_corr)-min(epm_v_corr)
; ;          xrange = [min(epm_v_corr)-0.05*vrange, $
; ;                    max(epm_v_corr)+0.05*vrange]
; ;          yrange = [min(subres-suberr) < 0, $
; ;                    max(subres+suberr) > 0]
; ;          plot, epm_v_corr, subres, xrange=xrange, /xstyle, $
; ;                yrange=yrange, $
; ;                xtitle='Velocity corrected by '+vcorr_epm, $
; ;                ytitle='Residual in log10(H0)', $
; ;                title='SNe II caibrated via EPM in Hubble flow', $
; ;                psym=4, symsize=symsz, thick=thick, charsize=chsize
; ;          ; if (printerr) then begin
; ;          ;    printf, lun, k, 'epm ',total(finite(epm_v_corr) eq 0),$
; ;          ;            total(finite(subres) eq 0), xrange, yrange, symsz, $
; ;          ;            thick, chsize, format='(i4,a10,2i4,4f12.4,3f8.3)'
; ;          ; endif
; ;          errplot, epm_v_corr, subres-suberr, subres+suberr
; ;          for j = 0, n_epm-1 do begin
; ;             dplus = [epm_v_corr, xrange[1]] - epm_v_corr[j]
; ;             dminus = [epm_v_corr, xrange[0]] - epm_v_corr[j]
; ;             whp = where (dplus gt 0, np)
; ;             whm = where (dminus lt 0, nm)
; ;             dplus_closest = min(dplus[whp])
; ;             dminus_closest = max(dminus[whm])
; ;             hoff = +1.
; ;             if (abs(dplus_closest) lt abs(dminus_closest)) then hoff =-1.
; ;             xpos = epm_v_corr[j] + ch_height * hoff * (3+hoff)/2.
; ;             ; results in an offset of 2.0 ch_height to the plus side
; ;             ; or -1.0 ch_height to the minus side
; ;             align = 0.5
; ;             ; stop
; ;             xyouts, xpos, subres[j], epm_name[j], $
; ;                     align=align, charsize=0.8*subchsize, orient=90.
; ; endif

return
end








      
