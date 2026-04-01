Function compute_alpha, mag, mag_err, vcorr=vcorr, vpec=vpec, vcmb=vcmb, $
                        vhel=vhel, q0=q0, j0=j0, vdisp=vdisp, covar=covar, $
                        alpha_err=alpha_err, chisq=chisq, ndof=ndof, $
                        ignore_offdiag=ignore_offdiag, $
                        redshift_range=redshift_range, $
                        simplified=simplified, debug=debug, fail=fail, $
                        optical=optical, ignore_cosmology=ignore_cosmology, $
                        verbose=verbose, details=details
;
; This function computes the alpha value, related to the intercept of the
; apparent magnitude-redshift relation in the Hubble flow, for a set
; of objects with measured apparent magnitudes and corrected velocities.
;
; The exact use of this relation depdns on how the luminosity relation
; is calibrated.  For SNe, the calibration yields a reference
; absolute magnitude in a specific filter (B for SNe Ia, I fo SNe II)
; and the expected value of the apparent magnitude is simply the
; reference absolute magnitude plus the luminosity distance modulus.
; For objects in the Tully-Fisher flow, the apparent magnitude is
; obtained as an absolute magnitude M, which relates to the velocity
; width - and therefore is not the same value for all objects.  As a
; consequence, the value of alpha for TF is based on m-M, not on m
; alone.

; Specific equations are found in DESCRIBE_V3 and in the main
; procedure.

; Upgrades:
; 1) Use the Davis+ 2019 formulation (same as Nils)
; 2) Include covariance if provided
; 3) Return chi^2 for the Hubble flow only (not including any
;    discrepancies in intercept vs Hubble value; that is in the network chi^2)
; The value f alpha thus defined is essentially the mean apparent magnitude
; at a fiducial distance.   It can be redefined as the mean magniude
; that minimizes the residualsin the Hubble Flow representation, and
; it can include covariance matrix.
  
; 6/24: Return values for individual objects in the Hubble flow.  This
; is in preparation for residuals plots and to include deviations in
; the full solution.  Returned in structure DETAILS.  Include:
; z_hd 

; 20250807: additions to handle specific cases
; optical: if set, adopt the optical convention that z = v/c.  Default
;    is to use z = sqrt((1+v/c)/(1-v/c))-1.  Megamasers and SBF adopt the
;    optical convention.  (however, megamasers and EPM do not use
;    compute_alpha at the moment)
; ignore_cosmology: if set, the quantity kz (usually the second-order
;    expansion of the luminosity distance vs redshift) is set to unity.
;    The reason is that SBF already include the cosmological
;    correction in their reported velocities, therefore we do not want
;    to apply the correction twice.
; In some future version, we can also add a switch to compute angular
;    diameter distance instead of luminosity distance.  Thi s is not
;    needed here because only megamasers use angular diameter
;    distances, and those are handled in the main code.

fail = 1   ; set to 0 before successful return
if (keyword_set(debug) eq 0) then do_debug = 0 else do_debug = debug
if (keyword_set(optical) eq 0) then optical = 0B
if (keyword_set(ignore_cosmology) eq 0) then ignore_cosmology=0B
if (keyword_set(simplified) eq 0) then simplified = 0B
speed_of_light = !const.c / 1.e3 ; !const.c in m/s; convert to km/s

; Set redshift_selection if redshift_range is set and at least one of
; its components is non-zero.

redshift_selection = 0B
if (keyword_set(redshift_range)) then begin
   if (redshift_range[0] ne 0 or redshift_range[1] ne 0) then $
      redshift_selection = 1B
endif
; Note that redshift_range applies to zcorr (or zhd)

; `1) Ugrade: exact formula for redshift-luminosity relation
  
if (keyword_set(q0) eq 0) then q0= -0.55d0
if (keyword_set(j0) eq 0) then j0= 1.d0

if (keyword_set(vdisp) eq 0) then vdisp = 150.d0
if (keyword_set(verbose) eq 0) then verbose=0
  
if (keyword_set(vcmb)) then goto, advanced
; The "advanced" version uses the more advanced calculation outlined by Nils
; Schoneberg and based on Davis+ 2019.  Also allows a covariance matrix
; to be given.  It is triggered by providing vcmb separately from
; vcorr.  If vcmb is not given, vcorr is required,  If vcmb is given,
; either vcorr or vpec can be provided. 


; Dispersion in velocity corrections.  Correlations are ignored for now.

if (optical) then zhd = vcorr/speed_of_light else $
   zhd = sqrt ((1+vcorr/speed_of_light)/(1-vcorr/speed_of_light))-1

useflag = replicate (1, n_elements(zhd))
if (redshift_selection) then useflag = (zhd ge redshift_range[0] and $
                                        zhd le redshift_range[1])
if (total(useflag) eq 0) then begin
   print, ' ERROR: no objects left after redshift selection'
   print, ' ALPHA set to -99'
   alpha=-99.d0
   alpha_err = 0
   fail=1
   if (do_debug) then stop
   details = {fail:1B}
   return, alpha
endif

if (ignore_cosmology) then kz=1.d0 else $
   kz = 1 + 0.5d0 * (1-q0) * zhd - 1.d0/6.d0 * (1-q0-3*q0^2+j0)*zhd^2

alphavec = alog10(vcorr) + alog10(kz) - 0.2*mag
errsqvec = 0.2^2*mag_err^2 + (alog10(vcorr+vdisp)-alog10(vcorr))^2

; The equation for logh0 in the non-advanced case is
; logh0 = log(c z_hd kz / D_L) = log10(c z_hd kz) - 0.2 * (m-M-25)
;       = log(c z_hd kz) - 0.2*m + 0.2*(M+25)
; compute_alpha returns the value of the terms on the rhs EXCEPT for
; M, which it does not know.  Also, for numerical convenience, the
; value does not include 0.2*25 = 5.
;

w = useflag/errsqvec
result = total (w*alphavec) / total(w)
alpha_err = 1.d0 / sqrt (total(w))
chisq = total(w*(alphavec-result)^2)
fail=0
ndof = total(useflag gt 0)-1
if (verbose gt 0) then print, result, alpha_err, chisq, ndof, format='(4f12.6)'
if (do_debug) then stop
wh = where(useflag, nwh)
details = {number:nwh, zhd:zhd[wh], kz:kz[wh], alphavec:alphavec[wh], $
           alphavar:errsqvec[wh], advanced:0B, simplified:simplified, $
           vdisp:vdisp, ignore_cosmology:ignore_cosmology, $
           optical:optical, has_covar:0B, covar:0, fail:0B}
return, result

advanced:

; Here if vcmb is set.

if (optical) then zcmb = vcmb/speed_of_light else $
   zcmb = sqrt ((1+vcmb/speed_of_light)/(1-vcmb/speed_of_light))-1

; compute zpec, zcorr from either vpec or vcorr
; Assume that vcorr (== vhd) = vcmb - vpec
;          or vpec = vcmb - vcorr
; Redshift cmbination formula: z3 = (1+z1)*(1+z2)-1

if (keyword_set(vpec) eq 0) then begin
   if (keyword_set(vcorr) eq 0) then begin
      print, ' Need one and only one of VPEC, VCORR; bailing out'
      if (do_debug) then stop
      return, 0
   endif
   ; Use the correct transformations based ion redshift, not velocity
   ; Here if vcorr is given
   if (optical) then zcorr = vcorr/speed_of_light else $
      zcorr = sqrt ((1+vcorr/speed_of_light)/(1-vcorr/speed_of_light)) - 1  
   zpec = (1+zcmb)/(1+zcorr)-1                  
   if (optical) then vpec = speed_of_light * zpec else $
      vpec = speed_of_light * ((1+zpec)^2-1)/((1+zpec)^2+1)
endif else begin
   ; Here if vpec is given
   if (optical) then zpec = vpec / speed_of_light else $
      zpec = sqrt ((1+vpec/speed_of_light)/(1-vpec/speed_of_light))-1
   zcorr = (1+zcmb)/(1+zpec)-1
   if (optical) then vcorr = speed_of_light * zcorr else $
      vcorr = speed_of_light * ((1+zcorr)^2-1)/((1+zcorr)^2+1)
endelse

if (keyword_set(vhel) eq 0) then begin
   ; Approximate vhel as vcmb (not quite right)
   vhel = vcmb
endif

if (optical) then begin
   zhel = vhel/speed_of_light
   zpec = vpec/speed_of_light
   zcmb = vcmb/speed_of_light
endif else begin
   zhel = sqrt ((1+vhel/speed_of_light)/(1-vhel/speed_of_light)) - 1
   zpec = sqrt ((1+vpec/speed_of_light)/(1-vpec/speed_of_light)) - 1
   zcmb = sqrt ((1+vcmb/speed_of_light)/(1-vcmb/speed_of_light)) - 1
endelse

zhd = (1+zcmb)/(1+zpec) - 1  ; equivalent to vcorr = vobs-vpec
t1 = (1+zhel) / (1+zhd)
t2 = speed_of_light * zhd
if (ignore_cosmology) then t3 = 1.d0 else $
   t3 = 1 + 0.5*(1-q0)*zhd - (1./6.d0)*(1-q0-3*q0^2+j0)*zhd^2
kz = t3
; print, t1, t2, t3
; stop

m_model = 5 * alog10 (t1*t2*t3)

useflag = replicate (1, n_elements(zhd))
if (redshift_selection) then useflag = (zhd ge redshift_range[0] and $
                                        zhd le redshift_range[1])
if (total(useflag) eq 0) then begin
   print, ' ERROR: no objects left after redshift selection'
   print, ' ALPHA set to -99'
   alpha=-99.d0
   alpha_err = 0
   fail=1
   if (do_debug) then stop
   return, alpha
endif

if (simplified) then begin
   zhd = vcorr / speed_of_light
   if (ignore_cosmology) then kz=1.d0 else $
      kz = 1 + 0.5d0 * (1-q0) * zhd - 1.d0/6.d0 * (1-q0-3*q0^2+j0)*zhd^2
   m_model = 5 * (alog10(vcorr) + alog10(kz))
endif

vel_var = (alog10(vcorr+vdisp)-alog10(vcorr))^2

; Now alpha is the result of optimizing the residuals between
; m_obs and m_model.  More specifically, alpha is the value that minimizes
; chi^2 = sum (m_obs-m_model+5*alpha)^T # invcov # (m_obs-m_model+5*alpha)
;
; This is a function of alpha only, and can be easily optimized by any
; of a number of simple methods.  In the current description, the
; problem is strictly linear in alpha, and with the covariance matrix
; approximation, is also assumed to be Gaussian.  Therefore a simple
; generalized least squares in one parameter suffices.
;
; If the covariance matrix is given, assume it applies to magnitudes
; only, and add uncorrelated, linearized velocity errors on the
; diagonal.  If the covariance matrix already includes velocity
; errors, simply give a very small number for vdisp.
;
; As a reminder, for the generalized least quares problem, the
; parameter array P that minimizes (D-A*P)^T # INVERT(C) # (D-A*P)
; (where C is the covariance matrix, A is the array of
; coefficients. and D is the array of data) is given by:
;    P = INVERT (A^T # S # A) # A^T # S # D
; where S = INVERT(C), and
;    VAR(P) = INVERT (A^T # S # A)
; Must use care because INVERT fails on a 1x1 matrix, which applies
; to A^T # S # A in a one-parameter case like this.
;

if (keyword_set(ignore_offdiag) eq 0) then ignore_offdiag = 0B
alphavec = 0.2d0*(m_model-mag)

if (keyword_set(covar) eq 0) then begin
   has_covar = 0B
   if (n_params() lt 2) then begin
      print, ' Need magnitude errors or covariance matrix; bailing out.'
      if (do_debug) then stop
      return, 0
   endif
   errsqvec = (0.2d0*mag_err)^2 + vel_var
   weight = useflag/errsqvec
   alpha = total (weight*alphavec) / total(weight)
   alpha_err = 1.d0 / sqrt (total(weight))
   chisq = total (weight*(alphavec-alpha)^2)
   ndof = nint(total(useflag gt 0)) - 1
   lcovv = 0
   wh = where (useflag, nwh)
   details = {number:nwh, zhd:zhd[wh], kz:kz[wh], alphavec:alphavec[wh], $
              alphavar:errsqvec[wh], has_covar:0B, covar:0, $
              vdisp:vdisp, ignore_cosmology:ignore_cosmology, $
              optical:optical, advanced:1B, simplified:simplified, fail:0B}
endif else begin
   ;
   ; In this description, the quantity to
   ; minimize is (m_model-mag-5*alpha).
   ; The COEFFS array is the [5, 5, 5,...] NEQ times
   ; LCOVV is the covariance matrix plus the velocity term
   ; The TEMP matrix is defined as TEMP = COEFFS^T # INV(LCOVV) # COEFFS
   ; 
   ; The input covariance matrix is in magnitudes.  The units of
   ; alphavec are 0.2*magnitude (or log10(distance), log10(velocity).
   ; Thus there is a factor 1/5, or 1/25 for the variance expressed
   ; in mag or mag^2.
   ;
   ; Note: before 20250624, the equations in this section were in
   ; magnitudes, unlike the log10(distance, velocity) used in different
   ; sections of the routine. Changed for uniformity.   

   has_covar = 1B
   wh = where (useflag gt 0, neq)
   lcovv = dblarr (neq, neq)
   if (ignore_offdiag eq 0) then begin
      for j = 0,neq-1 do for k=0,neq-1 do lcovv[j,k] = covar[wh[j],wh[k]]/25.d0
   endif else begin
      for k = 0,neq-1 do lcovv[k,k]=covar[wh[k],wh[k]]/25.d0
   endelse
   data = alphavec[wh]
   for k = 0, neq-1 do lcovv[k,k] += vel_var[wh[k]]
   alphavar = dblarr(neq)
   for k = 0, neq-1 do alphavar[k] = lcovv[k,k]
   coeffs = replicate (1.d0, neq)
   invcovv = la_invert(lcovv)
   temp = transpose(coeffs) # invcovv # coeffs
   if (n_elements(temp) eq 1) then invtemp = 1.d0/temp else $
      invtemp=la_invert(temp)
   alpha = invtemp # transpose(coeffs) # invcovv # data
   alpha=alpha[0]
   alpha_err = (sqrt(invtemp))[0]
   chisq = transpose(data-alpha*coeffs) # invcovv # (data-alpha*coeffs)
   ndof = neq-1
   details = {number:neq, zhd:zhd[wh], kz:kz[wh], alphavec:alphavec[wh], $
              alphavar:alphavar, has_covar:1b, covar:lcovv, $
              vdisp:vdisp, ignore_cosmology:ignore_cosmology, $
              optical:optical, advanced:1B, simplified:simplified, fail:0B}
endelse

fail = 0
if (do_debug) then stop
if (verbose gt 0) then print, alpha, alpha_err, chisq, ndof, format='(4f12.6)'

return, alpha

end
