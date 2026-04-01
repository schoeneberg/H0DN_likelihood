pro print_extrapars, det_list, descr_list, latexlabel=latexlabel, $
                     inlatex=inlatex, outdir=outdir
;                    hvaluecmb=hvaluecmb, $
;                    herrorcmb=herrorcmb, ndig_var=ndig_var, $
;                    ndig_table=ndig_table, ndig_sigma=ndig_sigma

nvars = n_elements (det_list)
if (keyword_set(outdir) eq 0) then outdir = './'
  
openw, lun4, outdir+'extrapars.tex', /get_lun
for k = 0, nvars-1 do begin
   dd = det_list[k]
   label = descr_list[k]
   if (inlatex[k]) then begin  ; only print line if requested by inlatex
      ; Extract short label (V00 etc)
      label=descr_list[k]
      icolon = strpos(label, ':')
      label1 = strmid(label,0,icolon)   ; variant number
      label2 = strmid(label, icolon+1, strlen(label)-icolon-1) ; variant descr
      ;      
      ; Create the strings for absolute calibration and alpha values
      sabs = ''
      salpha=''
      sncalib = ''
      ; Add for SN1a, SN2, SBF, TF
      ; sn1a
      if (dd.sn1a_details.do_sn1a) then begin
         sabs += string (' & $ ', dd.params_value[dd.iabs_sn1a], ' \pm ', $
                         sqrt(dd.params_var[dd.iabs_sn1a,dd.iabs_sn1a]), $
                         ' $ ', format='(2(a, f7.3), a)')
         salpha += string (' & $ ', dd.sn1a_details.alpha_sn1a_value, $
                           ' \pm ', dd.sn1a_details.alpha_sn1a_error, $
                           ' $ ', format = '(2(a, f7.3), a)')
         sncalib += string (' & ', dd.sn1a_details.n_sn1a, format='(a,i4)')
      endif else begin
         sabs += string (' & ', '---', format='(a, 9x, a, 9x)') 
         salpha += string (' & ', '---', format='(a, 9x, a, 9x)') 
         sncalib += string (' & ', 0, format='(a, i4)') 
      endelse
      ; SBF
      if (dd.sbf_details.do_sbf) then begin
         sabs += string (' & $ ', dd.params_value[dd.iabs_sbf], ' \pm ', $
                         sqrt(dd.params_var[dd.iabs_sbf,dd.iabs_sbf]), $
                         ' $ ', format='(2(a, f7.3), a)')
         salpha += string (' & $ ', dd.sbf_details.alpha_sbf_value, $
                           ' \pm ', dd.sbf_details.alpha_sbf_error, $
                           ' $ ', format = '(2(a, f7.3), a)')
         sncalib += string (' & ', dd.sbf_details.n_sbf, format='(a,i4)')
      endif else begin
         sabs += string (' & ', '---', format='(a, 9x, a, 9x)') 
         salpha += string (' & ', '---', format='(a, 9x, a, 9x)') 
         sncalib += string (' & ', 0, format='(a, i4)') 
      endelse
      ; SN2
      if (dd.sn2_details.do_sn2) then begin
         sabs += string (' & $ ', dd.params_value[dd.iabs_sn2], ' \pm ', $
                         sqrt(dd.params_var[dd.iabs_sn2,dd.iabs_sn2]), $
                         ' $ ', format='(2(a, f7.3), a)')
         salpha += string (' & $ ', dd.sn2_details.alpha_sn2_value, $
                           ' \pm ', dd.sn2_details.alpha_sn2_error, $
                           ' $ ', format = '(2(a, f7.3), a)')
         sncalib += string (' & ', dd.sn2_details.n_sn2, format='(a,i4)')
      endif else begin
         sabs += string (' & ', '---', format='(a, 9x, a, 9x)') 
         salpha += string (' & ', '---', format='(a, 9x, a, 9x)') 
         sncalib += string (' & ', 0, format='(a, i4)') 
      endelse
      ; TF
      if (dd.tf_details.do_tf) then begin
         sabs += string (' & $ ', dd.params_value[dd.iabs_tf], ' \pm ', $
                         sqrt(dd.params_var[dd.iabs_tf,dd.iabs_tf]), $
                         ' $ ', format='(2(a, f7.3), a)')
         salpha += string (' & $ ', dd.tf_details.alpha_tf_value, $
                           ' \pm ', dd.tf_details.alpha_tf_error, $
                           ' $ ', format = '(2(a, f7.3), a)')
         sncalib += string (' & ', dd.tf_details.n_tf, format='(a,i4)')
      endif else begin
         sabs += string (' & ', '---', format='(a, 9x, a, 9x)') 
         salpha += string (' & ', '---', format='(a, 9x, a, 9x)') 
         sncalib += string (' & ', 0, format='(a, i4)') 
      endelse
      printf, lun4, label1, sabs, salpha, sncalib, ' \\', format='(5a)'
   endif
endfor

close, lun4
free_lun, lun4

return
end
