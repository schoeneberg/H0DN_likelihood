Pro print_tables, det_list, descr_list, latexlabel=latexlabel, $
                  inlatex=inlatex, hvaluecmb=hvaluecmb, $
                  herrorcmb=herrorcmb, ndig_var=ndig_var, $
                  ndig_table=ndig_table, ndig_sigma=ndig_sigma

; Utility to print the tables used for diagnostics and for the paper.
; 1) Basic text table (default name vars_out.txt)
; 2) Main variants table (default name vartable.tex)
; 3) Details table (default name var_details.tex)
; 4) Definition of latex variables (default name vardef.tex)
;

; Number of decimal places for H0 in the table, H0 in the text, and
; sigma estimates

if (keyword_set(ndig_table) eq 0) then ndig_table = 2
if (keyword_set(ndig_var) eq 0)   then ndig_var = 3
if (keyword_set(ndig_sigma) eq 0) then ndig_sigma = 1

alpha_in_main = 0B
absmag_in_main = 0B
; If 1, then alpha and/or M0B for SNe Ia are included in the main
; variant table.  Otherwise they are skipped.

; Output file names

openw, lun, outdir+'vars_out.txt', /get_lun
openw, lun2, outdir+'vartable.tex', /get_lun
openw, lun3, outdir+'var_details.tex', /get_lun
openw, lun4, outdir+'vardef.tex', /get_lun


; Header for vartable.tex

if (absmag_in_main) then begin
   hh1_1 = '$ M_0 $ & '
   hh1_2 = ' mag & '
endif else begin
   hh1_1 = ''
   hh1_2 = ''
endelse

if (alpha_in_main) then begin
   hh2_1 = '$ \alpha_{SN Ia} $ & '
   hh2_2 = ' mag & '
endif else begin
   hh2_1 = ''
   hh2_2 = ''
endelse

printf, lun2, '\newcommand{\ptt}[1]{\parbox[t]{4cm}{\raggedright #1}}'
printf, lun2, '\begin{table}'
printf, lun2, '\caption{Main variants for $ H_0 $ calculation ' + $
              '\label{tab:variants}} % Updated 2025/10/01.'  
printf, lun2, '    \begin{center}'
printf, lun2, '        \small \begin{tabular}{lccrrrrrrl}'
printf, lun2, '\hline\hline'
printf, lun2, '\# & $ H_0 $ & 1--$\sigma$ &  $ \chi^2 $ & $ N_{\rm dof} $ '+$
              '& Reduced & ' + hh1_1 + hh2_1 + ' $ N_{\rm calib} $ '+$
              '& Description \\'
printf, lun2, '           & \multicolumn{2}{c}{$ \rm km/s/Mpc $} & & & '+$
        '$ \chi^2 $ & ' + hh1_2 + hh2_2 + ' (SNe Ia)   \\'
printf, lun2, '\hline'

; General definitions
printf, lun3, '% Adopted CMB values'
printf, lun3, '\newcommand{\Hvaluecmb}{'+$
           string(hvaluecmb,format='(f0.2)')+'}'
printf, lun3, '\newcommand{\Herrorcmb}{'+$
           string(herrorcmb,format='(f0.2)')+'}'
inot = 0
for k = 0, nvar-1 do begin
   dd = det_list[k]
   label = descr_list[k]
   if (dd.sn1a_details.do_sn1a) then begin
      alpha_sn1a = dd.sn1a_details.alpha_sn1a_value
      ealpha_sn1a = dd.sn1a_details.alpha_sn1a_error
      absmag_sn1a = dd.abs_sn1a_value
      eabsmag_sn1a = dd.abs_sn1a_error
      n_sn1a = dd.sn1a_details.n_sn1a
   endif else begin
      alpha_sn1a = 0.d0
      ealpha_sn1a = 0.d0
      absmag_sn1a = 0.d0
      eabsmag_sn1a = 0.d0
      n_sn1a = 0
   endelse

   if (dd.sn2_details.do_sn2) then begin
      alpha_sn2 = dd.sn2_details.alpha_sn2_value
      ealpha_sn2 = dd.sn2_details.alpha_sn2_error
      absmag_sn2 = dd.abs_sn2_value
      eabsmag_sn2 = dd.abs_sn2_error
      n_sn2 = dd.sn2_details.n_sn2
   endif else begin
      alpha_sn2 = 0.d0
      absmag_sn2 = 0.d0
      ealpha_sn2 = 0.d0
      eabsmag_sn2 = 0.d0
      n_sn2 = 0
   endelse
   ; stop
   printf, lun, k+1, dd.h0_value, dd.h0_error, dd.chi2, dd.ndof, $
           dd.chi2/dd.ndof, absmag_sn1a, absmag_sn1a, alpha_sn1a, $
           ealpha_sn1a, n_sn1a, $
           absmag_sn2, eabsmag_sn2, alpha_sn2, ealpha_sn2, n_sn2, label, $
           format='(i2, ") ", 2f10.4, 2x, f10.4, i7, f10.4, 3x, ' + $
           '2(2f10.3, 2f10.4, i4), 2x, a)'
   print, k+1, dd.h0_value, dd.h0_error, dd.chi2, dd.ndof, $
          dd.chi2/dd.ndof, absmag_sn1a, eabsmag_sn1a, alpha_sn1a, $
          ealpha_sn1a, n_sn1a, absmag_sn2, eabsmag_sn2, $
          alpha_sn2, ealpha_sn2, n_sn2, label, $
          format='(i2, ") ", 2f10.4, 2x, f10.4, i7, f10.4, 3x, ' + $
          '2(2f10.3, 2f10.4, i4), 2x, a)'
   ; For LaTeX
   form_table = '(f'+string(4+ndig_table,format='(i0)')+'.'+$
                 string(ndig_table, format='(i0)') + ')'
   form_var   = '(f'+string(4+ndig_var,format='(i0)')+'.'+$
                 string(ndig_var, format='(i0)') + ')'
   form_sigma = '(f'+string(4+ndig_sigma,format='(i0)')+'.'+$
                 string(ndig_sigma, format='(i0)') + ')'
   s_value = string(dd.h0_value, format=form_table)
   s_absmag = string(absmag_sn1a, format='(f7.3)')
   s_eabsmag = string(eabsmag_sn1a, format='(f7.3)')
   s_alpha = string(alpha_sn1a, format='(f7.4)')
   s_ealpha = string(ealpha_sn1a, format='(f7.4)')
   if (n_sn1a eq 0) then begin
      s_absmag = '  ---  '
      s_alpha = '  ---  '
      s_eabsmag = '  ---   '
      s_ealpha = '  ---   '
   endif
   icolon = strpos(label, ':')
   label1 = strmid(label,0,icolon)
   label2 = strmid(label, icolon+1, strlen(label)-icolon-1)
   ; Table uses 3 digits for variants: f6.3 for H0 value, f5.2 for error
   if (inlatex[k]) then printf, lun2, label1, ' & ', dd.h0_value, ' & ', $
           dd.h0_error, ' & ', dd.chi2, ' & ', dd.ndof, ' & ', $
           dd.chi2/dd.ndof, ' & ', s_absmag, ' & ', s_eabsmag, ' & ', $
           s_alpha, ' & ', s_ealpha, ' & ', $
           n_sn1a, ' & ', label2, ' \\', $
           format='(a5, a3, f6.3, a3, f5.3, a3, f8.4, a3, i3, a3, f7.5, ' + $
           ' 2(a3, a7), 2(a3, a8), a3, i3, a3, a, a3)'
   ; LaTeX variable definitions.  Use 2 decimal places for H0 and error
   llab = latexlabel[k]
   hsigma = abs(dd.h0_value-hvaluecmb)/sqrt(herrorcmb^2+dd.h0_error^2)
   printf, lun3, '% Quantities for '+label
   printf, lun3, '\newcommand{\Hvalue'+llab+'}{'+$
           string(dd.h0_value,format='(f0.2)')+'}'
   printf, lun3, '\newcommand{\Herror'+llab+'}{'+$
           string(dd.h0_error,format='(f0.2)')+'}'
   printf, lun3, '\newcommand{\Hsigma'+llab+'}{'+$
           string(hsigma,format='(f0.1)')+'}'
   printf, lun3, '\newcommand{\chisq'+llab+'}{'+$
           string(dd.chi2,format='(f0.3)')+'}'
   printf, lun3, '\newcommand{\ndof'+llab+'}{'+$
           string(dd.ndof,format='(i0)')+'}'
   printf, lun3, '\newcommand{\redchisq'+llab+'}{'+$
           string(dd.chi2/dd.ndof,format='(f0.4)')+'}'
   printf, lun3, '\newcommand{\absmagsn'+llab+'}{'+$
           string(absmag_sn1a,format='(f0.3)')+'}'
   printf, lun3, '\newcommand{\alphasn'+llab+'}{'+$
           string(alpha_sn1a,format='(f0.3)')+'}'
   printf, lun3, '\newcommand{\ncalibsn'+llab+'}{'+$
           string(n_sn1a,format='(i0)')+'}'
endfor

printf, lun2, '\hline\hline'
printf, lun2, '        \end{tabular}'
printf, lun2, '    \end{center}'
printf, lun2, '\end{table}'

close, lun
free_lun, lun
close, lun2
free_lun, lun2
close, lun3
free_lun, lun3

openw, lun4, 'extrapars.tex', /get_lun
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
                         sqrt(dd.params_var[dd.iabs_sn1a,dd_iabs_sn1a]), $
                         ' $ ', format='(2(a, f7.3), a)')
         salpha += string (' & $ ', dd.sn1a_details.alpha_sn1a_value, $
                           ' \pm ', dd.sn1a_details.alpha_sn1a_error, $
                           ' $ ', format = '(2(a, f7.3), a)')
         sncalib += string (' & ', dd_sn1a_details.n_sn1a, format='(a,i4)')
      endif else begin
         sabs += string (' & ', '---', format='(a, 9x, a, 9x)') 
         salpha += string (' & ', '---', format='(a, 9x, a, 9x)') 
         sncalib += string (' & ', 0, format='(a, i4)') 
      endelse
      ; SBF
      if (dd.sbf_details.do_sbf) then begin
         sabs += string (' & $ ', dd.params_value[dd.iabs_sbf], ' \pm ', $
                         sqrt(dd.params_var[dd.iabs_sbf,dd_iabs_sbf]), $
                         ' $ ', format='(2(a, f7.3), a)')
         salpha += string (' & $ ', dd.sbf_details.alpha_sbf_value, $
                           ' \pm ', dd.sbf_details.alpha_sbf_error, $
                           ' $ ', format = '(2(a, f7.3), a)')
         sncalib += string (' & ', dd_sbf_details.n_sbf, format='(a,i4)')
      endif else begin
         sabs += string (' & ', '---', format='(a, 9x, a, 9x)') 
         salpha += string (' & ', '---', format='(a, 9x, a, 9x)') 
         sncalib += string (' & ', 0, format='(a, i4)') 
      endelse
      ; SN2
      if (dd.sn2_details.do_sn2) then begin
         sabs += string (' & $ ', dd.params_value[dd.iabs_sn2], ' \pm ', $
                         sqrt(dd.params_var[dd.iabs_sn2,dd_iabs_sn2]), $
                         ' $ ', format='(2(a, f7.3), a)')
         salpha += string (' & $ ', dd.sn2_details.alpha_sn2_value, $
                           ' \pm ', dd.sn2_details.alpha_sn2_error, $
                           ' $ ', format = '(2(a, f7.3), a)')
         sncalib += string (' & ', dd_sn2_details.n_sn2, format='(a,i4)')
      endif else begin
         sabs += string (' & ', '---', format='(a, 9x, a, 9x)') 
         salpha += string (' & ', '---', format='(a, 9x, a, 9x)') 
         sncalib += string (' & ', 0, format='(a, i4)') 
      endelse
      ; TF
      if (dd.tf_details.do_tf) then begin
         sabs += string (' & $ ', dd.params_value[dd.iabs_tf], ' \pm ', $
                         sqrt(dd.params_var[dd.iabs_tf,dd_iabs_tf]), $
                         ' $ ', format='(2(a, f7.3), a)')
         salpha += string (' & $ ', dd.tf_details.alpha_tf_value, $
                           ' \pm ', dd.tf_details.alpha_tf_error, $
                           ' $ ', format = '(2(a, f7.3), a)')
         sncalib += string (' & ', dd_tf_details.n_tf, format='(a,i4)')
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


      
   
if (keyword_set(savefile)) then save, file=savefile, det_list, descr_list

return
end
