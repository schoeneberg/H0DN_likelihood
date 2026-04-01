Pro make_residual_plot_v3, plotfile, details, rota=rota, maxnames=maxnames, $
                           plabel=plabel

psave = !p
dsave = !d

; Extract variables out of DETAILS structure and substructures

; General
residuals          = details.residuals
h0_value           = details.h0_value
h0_error           = details.h0_error
abs_sn1a_value     = details.abs_sn1a_value
abs_sn1a_error     = details.abs_sn1a_error
mu_coma_value      = details.mu_coma_value
mu_coma_error      = details.mu_coma_error

; Hosts (was Rung2)
mas                = details.mas_details.mas
n2                 = details.host_data_details.n2
host2              = details.host_data_details.host2
mu_host2_error     = details.host_data_details.mu_host2_error
mas_index          = details.host_data_details.mas_index
mas_host_error     = details.host_data_details.mas_host_error
anchor_dist_error  = details.host_data_details.anchor_dist_error
hosts              = details.host_data_details.hosts

; sbf
n_sbf                 = details.sbf_details.n_sbf
if (n_sbf gt 0) then begin
   ieq_sbf_start      = details.sbf_details.ieq_sbf_start
   sbf_v_corr         = details.sbf_details.sbf_v_corr ; corrected velocities
   sbf_v_cmb          = details.sbf_details.sbf_v_cmb  ; corrected velocities
   vcorr_sbf          = details.sbf_details.vcorr_sbf  ; correction code
endif

; mm (was direct)

n_mm            = details.mm_details.n_mm
if (n_mm gt 0) then begin
   mm_name        = details.mm_details.mm_name
   ieq_mm_start   = details.mm_details.ieq_mm_start
   mm_logh_error  = details.mm_details.mm_logh_error
   mm_v_corr      = details.mm_details.mm_v_corr
   vcorr_mm       = details.mm_details.vcorr_mm ; correction code
endif

n_epm            = details.epm_details.n_epm
if (n_epm gt 0) then begin
   epm_name        = details.epm_details.epm_name
   ieq_epm_start   = details.epm_details.ieq_epm_start
   epm_logh_error  = details.epm_details.epm_logh_error
   epm_v_corr      = details.epm_details.epm_v_corr
   vcorr_epm       = details.epm_details.vcorr_epm ; correction code
endif

; coma
n_coma              = details.coma_details.n_coma             
if (n_coma gt 0) then begin
   coma_sn_name       = details.coma_details.coma_sn_name      
   ieq_coma_start     = details.coma_details.ieq_coma_start    
   n_coma              = details.coma_details.n_coma             
   coma_sn_err        = details.coma_details.coma_sn_err
   mu_coma_value      = details.coma_details.mu_coma_value
   mu_coma_error      = details.coma_details.mu_coma_error
endif

if (tag_exist (details, 'remove_unused')) then $
   remove_unused = details.remove_unused else remove_unused = 0B

if (remove_unused) then begin
   ; Match host2 with the hosts file, and remove if 
   ; the host is not used for any calibrators
   accept_host = bytarr (n2)
   for k = 0, n_elements(hosts) - 1 do begin
      wh = where (host2 eq hosts[k], nwh)
      if (details.host_used[k] gt 0 and nwh gt 0) then accept_host[wh] = 1B
   endfor
endif else accept_host = replicate (1B, n2)


; Plot setup
if (keyword_set(plabel) eq 0) then plabel=''

if (strlowcase(plotfile) eq 'screen') then begin
   set_plot, 'x'
   !p.font = -1
   chsize = 2.0
   subchsize = 0.9
   symsz = 1.4
   thick = 2
   screenplot = 1B
endif else begin
   set_plot, 'ps'
   !p.font = 0
   chsize = 0.9
   subchsize = 0.6
   symsz = 1.1
   thick = 1
   screenplot = 0B
   device, file=plotfile, /color, /landscape, /times
endelse

plot_subj = mas
; if (n_sbf gt 0) then plot_subj = [plot_subj , 'SBF']
; SBF plotting needs to be redone...
if (n_mm gt 0) then plot_subj = [plot_subj , 'MM']
if (n_epm gt 0) then plot_subj = [plot_subj , 'EPM']
if (n_coma gt 0) then plot_subj = [plot_subj , 'Coma']
nplots = n_elements(plot_subj)

if (keyword_set(maxnames) eq 0) then maxnames = 20
   ; if there are more than MAXNAMES systems in a plot, do not label them
if (nplots eq 1) then !p.multi=0
if (nplots gt 1 and nplots le 2) then !p.multi=[0, 2, 1]
if (nplots gt 2 and nplots le 4) then !p.multi=[0, 2, 2]
if (nplots gt 4 and nplots le 6) then !p.multi=[0, 3, 2]
if (nplots gt 6 and nplots le 9) then !p.multi=[0, 3, 3]
if (nplots gt 9 and nplots le 12) then !p.multi=[0, 4, 3]
if (nplots gt 12 and nplots le 16) then !p.multi=[0, 4, 4]
if (nplots gt 16) then !p.multi=[0,5,4]
yrange = [min(1.1*residuals) < (-0.1), max(1.1*residuals) > 0.1]
ydelta = yrange[1]-yrange[0]
ymid = 0.5*(yrange[1]+yrange[0])
for k = 0, nplots-1 do begin
   case plot_subj[k] of
      'SBF': begin
         wh = ieq_sbf_start + indgen(n_sbf) ; SBF residuals
         subres = residuals[wh]
         suberr = sqrt (sbf_logh_error^2-(0.2*cal_sbf_error)^2) 
         ; includes single-galaxy distance and velocity correction errors
         ; removes the covariant term in the distance calibration
         ; cal_sbf_error (the distance calibration error)
         ; has also been added to the distance error in Jensen+ 2021
         ; 
         xrange = [min(sbf_v_corr), 1.05*max(sbf_v_corr)]
         yrange = [min(subres-suberr) < 0, max(subres+suberr) > 0]
         plot, sbf_v_corr, subres, $
               xtitle='Velocity corrected by '+vcorr_sbf, $
               ytitle='Residual in log10(H0)', $
               title='SBF distance matched to Hubble flow', $
               psym=4, symsize=symsz, thick=thick, charsize=chsize
         oplot, [0,1.2*xrange[1]], [0,0]
         errplot, sbf_v_corr, subres-suberr, subres+suberr
         errplot, [0.99*yrange[1]],[-cal_sbf_error], $
                   [cal_sbf_error], psym=4, color=cgcolor('red'), $
                   thick = thick
         ; no individual galaxy names for now
      end
      
      'EPM': begin
         wh = ieq_epm_start + indgen(n_epm) ; EPM residuals
         subres = residuals[wh]
         suberr = epm_logh_error
         ; includes single-galaxy distance and velocity correction errors
         ; removes the covariant term in the distance calibration
         ; cal_sbf_error (the distance calibration error) has also 
         ; been added to the distance error from Jensen+ 2021
         vrange = max(epm_v_corr)-min(epm_v_corr)
         xrange = [min(epm_v_corr)-0.05*vrange, $
                   max(epm_v_corr)+0.05*vrange]
         yrange = [min(subres-suberr) < 0, $
                   max(subres+suberr) > 0]
         plot, epm_v_corr, subres, xrange=xrange, /xstyle, $
               yrange=yrange, $
               xtitle='Velocity corrected by '+vcorr_epm, $
               ytitle='Residual in log10(H0)', $
               title='SNe II caibrated via EPM in Hubble flow', $
               psym=4, symsize=symsz, thick=thick, charsize=chsize
         ; if (printerr) then begin
         ;    printf, lun, k, 'epm ',total(finite(epm_v_corr) eq 0),$
         ;            total(finite(subres) eq 0), xrange, yrange, symsz, $
         ;            thick, chsize, format='(i4,a10,2i4,4f12.4,3f8.3)'
         ; endif
         errplot, epm_v_corr, subres-suberr, subres+suberr
         oplot, xrange, [0,0]
         ; NEEDS FIXING!!!
         ; errplot, [n_epm-0.3], [-abs_sn1a_error], [abs_sn1a_error], psym=4, $
         ;          color=cgcolor('red'), thick=thick
         ; Vertical names; choose plus or minus
         ; depending on what is closest
         ch_height = float(!d.y_ch_size) / !d.x_vsize * $
                     0.8*subchsize * (xrange[1]-xrange[0]) / $
                     (!x.window[1]-!x.window[0])
         ; This is the height of a character rotated 90 deg 
         ; in units of the horiziontal coordinate
         for j = 0, n_epm-1 do begin
            dplus = [epm_v_corr, xrange[1]] - epm_v_corr[j]
            dminus = [epm_v_corr, xrange[0]] - epm_v_corr[j]
            whp = where (dplus gt 0, np)
            whm = where (dminus lt 0, nm)
            dplus_closest = min(dplus[whp])
            dminus_closest = max(dminus[whm])
            hoff = +1.
            if (abs(dplus_closest) lt abs(dminus_closest)) then hoff =-1.
            xpos = epm_v_corr[j] + ch_height * hoff * (3+hoff)/2.
            ; results in an offset of 2.0 ch_height to the plus side
            ; or -1.0 ch_height to the minus side
            align = 0.5
            ; stop
            xyouts, xpos, subres[j], epm_name[j], $
                    align=align, charsize=0.8*subchsize, orient=90.
         endfor
      end
      ;      
      'MM': begin
         wh = ieq_mm_start + indgen(n_mm) ; MM residuals
         subres = residuals[wh]
         suberr = mm_logh_error
         ; includes single-galaxy distance and velocity correction errors
         ; removes the covariant term in the distance calibration
         ; cal_sbf_error (the distance calibration error) has also 
         ; been added to the distance error from Jensen+ 2021
         vrange = max(mm_v_corr)-min(mm_v_corr)
         xrange = [min(mm_v_corr)-0.05*vrange, $
                   max(mm_v_corr)+0.05*vrange]
         yrange = [min(subres-suberr) < 0, $
                   max(subres+suberr) > 0]
         plot, mm_v_corr, subres, xrange=xrange, /xstyle, $
               yrange=yrange, $
               xtitle='Velocity corrected by '+vcorr_mm, $
               ytitle='Residual in log10(H0)', $
               title='Megamaser systems in Hubble flow', $
               psym=4, symsize=symsz, thick=thick, charsize=chsize
         ; if (printerr) then begin
         ;    printf, lun, k, 'mm ',total(finite(mm_v_corr) eq 0),$
         ;            total(finite(subres) eq 0), xrange, yrange, symsz, $
         ;            thick, chsize, format='(i4,a10,2i4,4f12.4,3f8.3)'
         ; endif
         errplot, mm_v_corr, subres-suberr, subres+suberr
         oplot, xrange, [0,0]
         ; NEEDS FIXING!!!
         ; errplot, [n_coma-0.3], [-abs_sn1a_error], [abs_sn1a_error], psym=4, $
         ;          color=cgcolor('red'), thick=thick
         ; Vertical names; choose plus or minus
         ; depending on what is closest
         ch_height = float(!d.y_ch_size) / !d.x_vsize * $
                     0.8*subchsize * (xrange[1]-xrange[0]) / $
                     (!x.window[1]-!x.window[0])
         ; This is the height of a character rotated 90 deg 
         ; in units of the horiziontal coordinate
         for j = 0, n_mm-1 do begin
            dplus = [mm_v_corr, xrange[1]] - mm_v_corr[j]
            dminus = [mm_v_corr, xrange[0]] - mm_v_corr[j]
            whp = where (dplus gt 0, np)
            whm = where (dminus lt 0, nm)
            dplus_closest = min(dplus[whp])
            dminus_closest = max(dminus[whm])
            hoff = +1.
            if (abs(dplus_closest) lt abs(dminus_closest)) then hoff =-1.
            xpos = mm_v_corr[j] + ch_height * hoff * (3+hoff)/2.
            ; results in an offset of 2.0 ch_height to the plus side
            ; or -1.0 ch_height to the minus side
            align = 0.5
            ; stop
            xyouts, xpos, subres[j], mm_name[j], $
                    align=align, charsize=0.8*subchsize, orient=90.
         endfor
      end
      'Coma': begin
         wh = ieq_coma_start + indgen(n_coma) ; Coma residuals
         ; these are supernova distance moduli minus the Coma
         ; distance modulus
         subres = residuals[wh]
         source_names = coma_sn_name
         suberr = coma_sn_err
         xrange = [-1, n_coma]
         yrange = [min(subres-coma_sn_err) < 0, $
                   max(subres+coma_sn_err) > 0]
         ; stop
         plot, indgen(n_coma), subres, xrange=xrange, yrange=yrange, $
               ytitle='Distance modulus residuals', $
               title='Type Ia SNe in Coma', $
               psym=4, symsize=symsz, thick=thick, charsize=chsize, $
               xticks=n_coma-1, xtickv=indgen(n_coma), $
               xtickname=replicate('  ', n_coma) ; blank x label
         ; if (printerr) then begin
         ;    printf, lun, k, '   Coma ', total(finite(subres) eq 0), $
         ;            xrange, yrange, symsz, $
         ;            thick, chsize, format='(i4,a10,4x,i4,4f12.4,3f8.3)'
         ; endif
         errplot, indgen(n_coma), subres-suberr, subres+suberr
         oplot, [-1,n_coma], [0,0]
         errplot, [n_coma-0.3], [-mu_coma_error], [mu_coma_error], $
                  color=cgcolor('red'), thick=thick
         ; add SN names as x labels if there are not too many
         ; code from mas section
         if (n_coma le maxnames) then begin
            lchsize = subchsize
            label_levels = 1
            deltastep = 0.06
            deltaoffset = 0.08
            if (n_coma gt 5) then begin
               label_levels = 2
               deltaoffset *= 0.8
            endif
            if (n_coma gt 10) then begin
               label_levels = 3
               lchsize *= 0.9
               deltaoffset *= 0.8
            endif
            if (n_coma gt 20) then begin
               label_levels = 4
               lchsize *= 0.8
               deltaoffset *= 0.9
               deltastep = 0.04
            endif
            angle = 0
            align = 0.5
            if (keyword_set(rota)) then begin
               angle=90
               label_levels = 1
               align=1.05
               deltaoffset = 0.0
               lchsize = 0.7*subchsize
            endif 
            for j=0,n_coma-1 do begin
               xpos = j
               ypos = !y.crange[0] - deltaoffset*ydelta - $
                      (j mod label_levels)*deltastep*ydelta
               ; m1 = check_math()
               xyouts, xpos, ypos, source_names[j], align=align, $
                       charsize=lchsize, orient=angle
               ; m2 = check_math()
               ; if (printerr) then printf, lun, k, j, plot_subj[k], $
               ;    xpos, ypos, subhost[j], lchsize, align, angle, m1, m2, $
               ;    format='(2i4,2x,a20,2f10.2,2x,a10,3f10.2, 2i10)'
            endfor
         endif            
      end
      else: begin
         ; it's one of the MAS; k is the index in the MAS array
         ; (this overloads and assumes that al MAS are at the start;
         ; using a separate variable might be preferable
         ddh = details.host_data_details
         wh = where (ddh.mas_index eq k and accept_host, nsub)
         if (nsub le 0) then goto, skip_plot
         subres = details.residuals[wh]
         subhost = ddh.host2[wh]
         suberrors = ddh.mu_host2_error[wh]
            ; This is just the error bar for the individual host measurement
         common_error = sqrt(ddh.mas_host_error[wh[0]]^2 + $
                             ddh.anchor_dist_error[wh[0]]^2)
            ; This is the error bar in common to all measurements
            ;  for this group: anchor + anchor_method_source
            ; It is obtained from the first element as it should be 
            ; the same for all elements
         label = details.mas_details.mas[k]
         ; Warning: if nsub > 60, it cannot be used for xticks etc.
                                ; Use nsubsm obtained by dividing nsub
                                ; by 2 until it is less than 60
         nsubsm = nsub
         while (nsubsm gt 60) do nsubsm = nsubsm/2
         plot, indgen(nsub), subres, xrange=[-0.4,nsub-0.6], $
               /xstyle, yrange=yrange, /ystyle, title=label, $
               charsize=chsize, psym=4, symsize=symsz, thick=thick, $
               xticks=nsubsm-1, xtickv=indgen(nsubsm), $
               xtickname=replicate('  ', nsubsm) ;  subhost
         errplot, findgen(nsub), subres-suberrors, subres+suberrors
         oplot,[-1,nsub+1],[0,0]
         oplot, [nsub-0.8], [0], psym=4, symsize=0.8*symsz, $
                thick=thick, color=cgcolor('red')
         ; indicative common error bar
         errplot, [nsub-0.8], [-common_error], [common_error], $
                  thick=thick, color=cgcolor('red')
         if (nsub le maxnames) then begin
            lchsize = subchsize
            label_levels = 1
            deltastep = 0.06
            deltaoffset = 0.08
            if (nsub gt 5) then begin
               label_levels = 2
               deltaoffset *= 0.8
            endif
            if (nsub gt 10) then begin
               label_levels = 3
               lchsize *= 0.9
               deltaoffset *= 0.8
            endif
            if (nsub gt 25) then begin
               label_levels = 4
               lchsize *= 0.8
               if (screenplot eq 0) then lchsize *= 0.9
               deltaoffset *= 0.9
               deltastep = 0.04
            endif
            angle = 0
            align = 0.5
            if (keyword_set(rota)) then begin
               angle=90
               label_levels = 1
               align=1.05
               deltaoffset = 0.0
               lchsize = subchsize*0.7
            endif 
            for j=0,nsub-1 do begin
               xpos = j
               ypos = yrange[0] - deltaoffset*ydelta - $
                      (j mod label_levels)*deltastep*ydelta
               ; m1 = check_math()
               xyouts, xpos, ypos, subhost[j], align=align, $
                       charsize=lchsize, orient=angle
               ; m2 = check_math()
               ; if (printerr) then printf, lun, k, j, plot_subj[k], $
               ;    xpos, ypos, subhost[j], lchsize, align, angle, m1, m2, $
               ;    format='(2i4,2x,a20,2f10.2,2x,a10,3f10.2, 2i10)'
            endfor
         endif
      end
      ; Old commands:
      ; if (subres[j] gt ymid) then ypos =
      ; subres[j] + 0.07*ydelta else $
      ; ypos = subres[j] - 0.12*ydelta
      ; xpos = j
      ; if (j eq 0) then xpos = 0.2
      ; if (j eq nsub-1) then xpos = nsub-1.2
      ; xyouts, xpos, ypos, subhost[j], align=0.5, charsize=subchsize
   endcase
   skip_plot:
endfor

if (plabel ne '') then begin
   !p.multi=0
   xyouts, 0.5, -0.04, /norm, plabel, align=0.5, charsize=1.2
endif

if (screenplot eq 0) then device, /close

; tn = tag_names(!p)
; for k = 0, n_elements(tn)-1 do !p.k = psave.k
; !p = psave
; !d = dsave

return
end
