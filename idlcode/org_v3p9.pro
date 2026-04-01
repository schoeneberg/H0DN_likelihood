Pro org_v3p9, config=config, anchor_data=anchor_data, host_data=host_data, $
              sn1a_calib=sn1a_calib, $
              sn1a_hf=sn1a_hf, sn1a_covar=sn1a_covar, $
              sn1a_ignore_offdiag=sn1a_ignore_offdiag, $
              sn1a_intrinsic=sn1a_intrinsic, $ 
              sn2_calib=sn2_calib, sn2_hf=sn2_hf, $
              tf_calib=tf_calib, tf_hf=tf_hf, groups_file=groups_file, $
              alpha_sn1a_value=alpha_sn1a_value, $
              alpha_sn1a_error=alpha_sn1a_error, $
              alpha_sn2_value=alpha_sn2_value, $
              alpha_sn2_error=alpha_sn2_error, $
              alpha_tf_value=alpha_tf_value, $
              alpha_tf_error=alpha_tf_error, $
              exclude_list=exclude_list, sn_exclude=sn_exclude, $
              include_list=include_list, sn_include=sn_include, $
              details=details, coma=coma, mm=mm, epm=epm, $
              sbf_calib=sbf_calib, sbf_hf=sbf_hf, q0=q0, j0=j0, $
              prior_coma=prior_coma, sigma_mag_fp=sigma_mag_fp, $
              h0_coma_nom=h0_coma_nom, d_coma_nom=d_coma_nom, $
              vcorr_mm=vcorr_mm, vcorr_epm=vcorr_epm, $
              vcorr_sbf=vcorr_sbf, vcorr_sn1a=vcorr_sn1a, $
              vpec_error=vpec_error, sbf_vpec_error=sbf_vpec_error, $
              logfile=logfile, fitsfile=fitsfile, $
              plotfile=plotfile, savefile=savefile, $
              min_redshift=min_redshift, max_redshift=max_redshift, $
              basedir=basedir, plabel=plabel, force_diag=force_diag, $
              debug=debug, nocull=nocull

; Version 3.9: Initial release version.
;
; GENERAL COMMENTS
; ================
;
; Approach based on published data.  Results in some duplications in
; data, identifiable as singlarities (linearly dependent rows/columns)
; in the data covariance matrix.  This is handled by using the
; Moore-Penrose inverse, which zeroes out the part of the problem
; corresponding to very small eigenvalues (numerically consistent with zero).
;
; In this approach, all the information regarding hosts is read at the
; outset, and all host distances are optimizable parameters.
;
; Then there are two main categories of data:
; - Host-based data.  These are objects for which a host distance is
;   determined and used to calibrate the reference absolute luminosity
;   of the class.  In some cases these are truly hosts (e.g., SNe Ia
;   and SNe II).  In other cases the host can be the same as the
;   galaxy itself (TF and many SBF calibrators).  In all cases the
;   distance of the host is used to help calibrate the distance
;   indicator.  The calibration is an optimizable parameter.
;  
; - Hubble flow data.  These are distance indicator at redshift large
;   enough that the peculiar velocity is small, and the object can be
;   considered part of the Hubble flow.  Some peculiar velocity
;   correction is usually still needed.  For objects in the Hubble
;   flow, only the value of alpha (an intercept of the apparent
;   magnitude/luminosity distance modulus, equivalent to the quantity
;   called a_X in SN analysis) is computed, together with its error
;   and the related chi^2.  The value of alpha is then used n the H0
;   optimization.
;
; In addition to these main categories, which apply to SNe Ia, SNe II,
; and TF, there are some special cases:
;
; - Direct distance indicators.  These are objects for which the
;   distance is obtained directly from a physics-based analysis, and
;   does not depend on anchors or host distances.  Examples include
;   megamasers (an angular diameter distance obtained from the
;   keplerian motion of sources near the center of the system) and SNe
;   II calibrated with the expanding photosphere method.  These have
;   no optimizable parameter associated with them, and directly
;   provide a constraint on H0.  We treat them as separate equations,
;   although it would also be acceptable to conflate them into an
;   equivalent intercept (essentially a value of alpha).
; - Fundamental Plane.  This is equivalent to other systems in that it
;   has its own set of Hubble Flow measurement, but the only host used
;   to date is the Coma Cluster, and the Hubble Flow information is
;   not currently handled directly.  The distance to the Coma Cluster is
;   obtained using the equivalent of a host file in which all objects
;   (Type Ia supernovae) are in the same host, and an additional SBF
;   constraint on Coma distance is used.  Consequently the Coma
;   Cluster adds one optimizable parameter (the Coma distance itself)
;   and two equations, one for the SBF constraint and one for the map
;   to the Hubble flow.  The latter is in the form of a delta distance
;   compared to a nominal Coma calibration.
; - SBF.  These are treted as precalibrated distance estimates, based
;   on TRGB measurements; the estimated distance moduli are covariant
;   with all host distance equations using the same measurements.
;   Ideally one would include calibrators in the distance network, but
;   that has not been done yet, in part because some calibrators share
;   a distance (to the Virgo or Fornax clusters) which require
;   aditional parameters and equations.  Ad hoc code for this case is
;   deferred to V3.3. 
;
; From a coding perspective, we first read all possible hosts; then
; read additional distance indicator information, including both a
; host/calibrator file and a Hubble flow file (these go in pairs); and
; finally add the specialized information for megamasers and SNe II
; using EPM.  For each block we add the relevant parameters, typically
; the absolute reference magnitude of the calibrator, and the
; corresponding equations, one per calibrator plus one for the Hubble
; Flow mapping to the Hubble constant.  We keep track of each equation
; with three string arrays: eq_descr, which names the equation;
; eq_shape, which identifies the terms in the equation and
; specifically which parameters are involved; and eq_covar, which
; describes which terms contribute to the covariance between data.
; The covariance needs to be tracked carefully.
;
; An alternate approach would make each model quantity, including
; observed magntiudes etc, a separate optimizable parameter,
; constrained by available measurements.  As long as the measurements
; are independent, the covariance matrix would be diagonal.  We did
; not use this approach because the relevant observational information
; is not always directly available, and would often require
; interpretation of the published results.  It would also increase
; significantly the number of optimizable parameters, and create a
; substantial burden in interpreting the corrections that often have
; to be applied to measurements.  We may attempt this approach in the
; future; some collaborators have discussed potential implementations,
; but difficulties remain.
;
; The distance measurements for each host, based on anchor and method,
; are all collected into a single file (default name host_data.dat),
; which also includes codes for the source of the measurement.  We
; generally treat measurements coming from different sources as
; uncorrelated, although some exceptions exist, e.g., for the SBF
; calibration.   We define equations and covariances based on anchor,
; method, and source for each measurement.  In this case, the diagonal
; element of the covariance is the uncertainty from the equation,
; which by definition does NOT include anchor geometry or anchor
; measurement uncertainty; the anchor measurement uncertainty is given
; in the corresponding anchor entry for that method and source.  
; The anchor geometry is then added to the diagonal and off-diagonal
; elements for all host equations using that anchor; and the anchor
; method zero point is added to the covariance matrix for all
; host equations using that anchor, method, and source.
; In this approach, there is not a
; separate variable for the anchor distance or for the anchor method
; zero point; they are subsumed into the error calculation.  The
; variable is the host distance.
;
; As far as the supernova calibrators, we can use the same method to
; link the host distance directly to H0.  In this case the equation
; for H0 has an uncertainty based on the supernova photometry and on
; the intrinsic supernova scatter; those will not be covariant with
; any other uncertainties.
;

; Coding note: the interface allows for a configuration file, called
; config.ini by default, which contains all necessary parameters.
; This is done for compatibility with the python code
; written by Emre Ozulker.  Inputs from config.ini are handled by saving
; name/value pairs in two string arrays, cpars and cvals, which are
; then interrogated by the routine input_param.  If a parameter is
; not given explicitly in the call, the value is obtained from cpars
; and cvals (combining multiple entries with the same value of cpars
; into a 1-d array); a default value is passed to input_param for
; each input parameter.  Note that the configuration file option has
; not been thoroughly debugged.

; Note that config must be set to 'none' (or to a value corresponding 
; to a non-existent file) if only parameters at invocation are
; used; otherwise the value in the configuration file could override
; the expected defaults.
  
; Special variables control the output (in part because screen
; output is not compatible with the distributable version of the
; code):
; logfile  : destination of essential results; filename or 'screen' to
;            print at the terminal.  Default: no logging.
; plotfile : destination of the plot; file or 'screen' to plot at the
;            terminal.  Default: no plotting.
; savefile : where to save the DETAILS structure.  Default: no save.
; fitsfile : name of the FITS file to include the relevant arrays.
;            Default: no FITS output.
; Note that make_residual_plot_v3 and print_residuals have been
; designed to work with the DETAILS structure.  print_residuals is not
; called from this program; make_residual_plot_v3 can be, but can also
; be called separately, passing DETAILS for complete information.
;

; Data are expected to be in the current directory or have a specified
; path (absolute or relative) given in the file name for each dataset.
; It is however possible to change the default data location via the
; variable BASEDIR (defaults to './'); this default directory is added to
; each file name that does not already have a path included.  However,
; BASEDIR does not change the starting directory for relative paths;
; all relative pathnames are relative to the current directory when
; the program is executed.  Also, BASEDIR does not affect the location
; of the configuration file; if the configuration file is in a different 
; directory, its name must include the correct path.

; ==========
; = SET UP =
; ==========
  
if (keyword_set(debug)) then do_debug = 1B else do_debug = 0B

if (keyword_set(config) eq 0) then config = 'config.ini'
ff = file_search (config, count=count)
if (config eq '' or config eq 'none' or count le 0) then begin
                                ; No configuration file; parameters
                                ; are either default or from command
                                ; line
   cpars = ''
   cvals = ''
endif else begin
   readcol, 'config.ini', cpars, cvals, format='a,a', comment='#', /silent
endelse

; Set parameters or their defaults
; cpars for control parameters
; cvals for their values

if (keyword_set(logfile) eq 0) then $
   logfile = input_param ('logfile', 'log.txt', cpars, cvals)
if (strlowcase(logfile) eq 'none' or logfile eq '') then $
   doprint=0B else doprint=1B
if (doprint and strmid(logfile,0,1) eq '+') then begin
   print_append = 1B
   logfile=strmid(logfile, 1, strlen(logfile)-1)
endif else print_append = 0B

if (keyword_set(fitsfile) eq 0) then $
   fitsfile = input_param ('fitsfile', 'none', cpars, cvals)
if (strlowcase(fitsfile) eq 'none' or fitsfile eq '') then $
   dofits=0B else dofits=1B

if (keyword_set(savefile) eq 0) then $
   savefile = input_param ('savefile', 'none', cpars, cvals)
if (strlowcase(savefile) eq 'none' or savefile eq '') then $
   dosave=0B else dosave=1B

if (keyword_set(basedir) eq 0) then $
   basedir = input_param ('basedir', './', cpars, cvals)

if (keyword_set(anchor_data) eq 0) then $
   anchor_data = input_param ('anchor_data', 'none', cpars, cvals)
   chg_default_dir, anchor_data, basedir

if (keyword_set(host_data) eq 0) then $
   host_data =   input_param ('host_data', 'host_data.dat', $
                              cpars, cvals)
   chg_default_dir, host_data, basedir

if (keyword_set(sn1a_calib) eq 0) then $
   sn1a_calib =   input_param ('sn1a_calib', 'sn1a_calib.dat', $
                               cpars, cvals)
   chg_default_dir, sn1a_calib, basedir

; sn1a_intrinsic is the extra intrinsic dispersion for SN1A calibrators.
; It is generally included in the data files but it is not included
; for IR.  The recommended values (from Lluis Galbany) are 0.096 mag for H
; and 0.125 mag for J.  This could be included in the file header if
; headers are implemented.
if (keyword_set(sn1a_intrinsic) eq 0) then $
      sn1a_intrinsic = input_param ('sn1a_intrinsic', 0.d0, cpars, cvals)

; The parameter sn1a_ignore_offdiag instructs compute_alpha to ignore
; the off-diagonal terms in the covariance matrix.  This parameter is
; ignored unless a covariance matrix is given for SNe Ia in the Hubble
; flow.  We could issue a warning if the parameter is set and the
; covariance matrix is not provided.  At the moment the parameter is
; silently ignored.
if (keyword_set(sn1a_ignore_offdiag) eq 0) then $
   sn1a_ignore_offdiag = input_param ('sn1a_ignore_offdiag', 0, cpars, cvals)

if (strlowcase(sn1a_calib) ne 'none' and sn1a_calib ne '') then $
   do_sn1a = 1B else do_sn1a = 0B

if (keyword_set(sn1a_hf) eq 0) then $
   sn1a_hf = input_param ('sn1a_hf', 'none', cpars, cvals)
   chg_default_dir, sn1a_hf, basedir
if (keyword_set(sn1a_covar) eq 0) then $
   sn1a_covar = input_param ('sn1a_covar', 'none', cpars, cvals)
   chg_default_dir, sn1a_covar, basedir

; R22: a_b = 0.714                     ; baseline from R22 Table 5
; Assume sigma_a_b = 0.002
; Need to change this behavior to allow for SN1a to not be included.
; Currently do_sn1a is set to zero only if sn1a claibrators are
; excluded.  But we may want to include sn1a calibrators but not
; sn1a in the Hubble flow because of the Coma calibration.  

if (strlowcase(sn1a_hf) eq 'none' or sn1a_hf eq '') then begin
   read_sn1a_hf = 0B
   if (keyword_set(alpha_sn1a_value) eq 0) then alpha_sn1a_value = $
      double (input_param('alpha_sn1a_value', 0.714d0, cpars, cvals))
   if (keyword_set(alpha_sn1a_error) eq 0) then alpha_sn1a_error = $
      double(input_param('alpha_sn1a_error', 0.002d0, cpars, cvals))
   alphadet_sn1a = {read_sn1a_hf:0, alpha_sn1a_error:alpha_sn1a_error}
   name_sn1a_hf = 'none'
endif else read_sn1a_hf = 1B
; NOTE: sn1a_covar is ignored unless read_sn1a_hf is true

if (keyword_set(sn2_calib) eq 0) then $
   sn2_calib = input_param ('sn2_calib', 'none', cpars, cvals)

if (strlowcase(sn2_calib) ne 'none' and sn2_calib ne '') then $
   do_sn2 = 1B else do_sn2 = 0B
if (keyword_set(sn2_hf) eq 0) then $
   sn2_hf = input_param ('sn2_hf', 'none', cpars, cvals)
if (sn2_calib eq '1') then begin
   sn2_calib = 'sn2_calib.dat'
   sn2_hf = 'sn2_hf.dat'
endif

chg_default_dir, sn2_calib, basedir
chg_default_dir, sn2_hf, basedir


; From the message by Thomas de Jaeger:
; a_b = 0.217 +/- 0.008  (check values)
; H0:                        72.904 4.312 
; Mi:                      -16.776 0.121
; -5 ai:                    -1.086 0.039
; and this assumes that magnitudes for SN II are in the I band
if (strlowcase(sn2_hf) eq 'none' or sn2_hf eq '') then begin
   read_sn2_hf = 0B
   if (keyword_set(alpha_sn2_value) eq 0) then alpha_sn2_value = $
      double (input_param('alpha_sn2_value', 0.217, cpars, cvals))
   if (keyword_set(alpha_sn2_error) eq 0) then alpha_sn2_error = $
      double(input_param('alpha_sn2_error', 0.008, cpars, cvals))
endif else read_sn2_hf = 1B

if (keyword_set(min_redshift) eq 0) then $
   min_redshift = input_param ('min_redshift', 0.d0, cpars, cvals)
if (keyword_set(min_redshift) eq 0) then $
   max_redshift = input_param ('max_redshift', 0.d0, cpars, cvals)
redshift_range=[min_redshift, max_redshift]

if (keyword_set(tf_calib) eq 0) then $
   tf_calib = input_param ('tf_calib', 'none', cpars, cvals)
if (strlowcase(tf_calib) ne 'none' and tf_calib ne '') then $
   do_tf = 1B else do_tf = 0B
if (keyword_set(tf_hf) eq 0) then $
   tf_hf = input_param ('tf_hf', 'none', cpars, cvals)
if (tf_calib eq '1') then begin
   tf_calib = 'tf_calib.dat'
   tf_hf = 'tf_hf.dat'
endif
chg_default_dir, tf_calib, basedir
chg_default_dir, tf_hf, basedir

; Values for TF Hubble Flow!!!

if (keyword_set(q0) eq 0) then $
   q0 = double (input_param ('q0', -0.55d0, cpars, cvals))
if (keyword_set(j0) eq 0) then $
   j0 = double (input_param ('j0', 1.d0, cpars, cvals))
if (keyword_set(exclude_list) eq 0) then $
   exclude_list = input_param ('exclude_list', '', cpars, cvals)
if (keyword_set(include_list) eq 0) then $
   include_list = input_param ('include_list', '', cpars, cvals)
; If include_list is set, then exclude_list will be ignored and a
; warning given.  This is handled in the inclusion and exclusion section.

if (keyword_set(sn_exclude) eq 0) then $
   sn_exclude = input_param ('sn_exclude', '', cpars, cvals)
if (keyword_set(sn_include) eq 0) then $
   sn_include = input_param ('sn_include', '', cpars, cvals)

if (keyword_set(coma)) then coma_found = 1B else begin
   coma = input_param ('coma', '', cpars, cvals, coma_found)
endelse
if (coma eq '1') then coma = 'coma_sn.dat'
chg_default_dir, coma, basedir


; Coma-FP handled differentially with respect to Said+ 2024.
; d_coma_nom and h0_coma_nom are their values of distance to Coma and
; the corresponding value of H0.  We rescale the inferred H0 based on
; the distance estimate for Coma from the distance network solution.
; prior_coma is the uncertainty in the SBF distance of Coma in
; magnitudes.  This term can be set to a large value to remove the
; dependence of Coma distance on the SBF measurement.
;
if (keyword_set(h0_coma_nom) eq 0) then $
   h0_coma_nom = input_param ('h0_coma_nom', 76.05d0, cpars, cvals)
if (keyword_set(d_coma_nom) eq 0) then $
   d_coma_nom = input_param ('d_coma_nom', 99.1d0, cpars, cvals)
if (keyword_set(prior_coma) eq 0) then $
   prior_coma = input_param ('prior_coma', 0.1345d0, cpars, cvals)

; 
; FP uncertainty in Hubble flow
if (keyword_set(sigma_mag_fp) eq 0) then $
   sigma_mag_fp = input_param ('sigma_mag_fp', 0.0371, cpars, cvals)
; Uncertainty in the projection of FP calibration to Hubble flow
; Value of 0.0371 from combining in quadrature 0.017 (FP width in
; Hubble flow) and 0.033 (FP uncertainty width in Coma) from Scolnic+ 2025


; Logic for SBF mirrors that for SN2:
; sbf_calib defaults to 'none' (if not given)
; if given with a value other than 'none' or '':
; if it is '1' (the string) of 1 (the number)
;    then the file name defaults to 'sbf_calib.dat'

if (keyword_set(groups_file) eq 0) then $
   groups_file = input_param ('groups_file', 'none', cpars, cvals)
if (strlowcase(groups_file) eq 'none' or groups_file eq '') then $
   do_groups = 0B else do_groups = 1B
; No default for groups_file.
chg_default_dir, groups_file, basedir


if (keyword_set(sbf_calib) eq 0) then $
   sbf_calib = input_param ('sbf_calib', 'none', cpars, cvals)
; test if sbfcalib is a string.  if not, assume it was set to 1 
; (as in /sbf_calib)
; and assign it the default name for the sbf_calib file
sbf_calib_default_filename = 'sbf_calib.dat'

if (size(sbf_calib, /type) ne 7) then begin
   sbf_calib = sbf_calib_default_filename
   do_sbf = 1B
endif else begin
   if (strlowcase(sbf_calib) ne 'none' and sbf_calib ne '') then $
      do_sbf = 1B else do_sbf = 0B
   if (keyword_set(sbf_hf) eq 0) then $
      sbf_hf = input_param ('sbf_hf', 'none', cpars, cvals)
   if (sbf_calib eq '1') then sbf_calib = sbf_calib_default_filename
endelse
if (do_sbf) then begin
   if (keyword_set(sbf_hf) eq 0) then begin
      sbf_hf = input_param ('sbf_hf', 'sbf_hf.dat', cpars, cvals)
   endif else begin
      if (size(sbf_hf, /type) ne 7) then sbf_hf = 'sbf_hf.dat'
   endelse
endif

chg_default_dir, sbf_calib, basedir
chg_default_dir, sbf_hf, basedir

; What to do if do_sbf is set but sbf_hf is not?  We may need a
; check to avoid this situation.

if (keyword_set(mm)) then mm_found = 1B else begin
   mm = input_param ('mm', '', cpars, cvals, mm_found)
endelse
if (mm eq '1') then mm = 'megamaser.dat'
chg_default_dir, mm, basedir

if (keyword_set(epm)) then epm_found = 1B else begin
   epm = input_param ('epm', '', cpars, cvals, epm_found)
endelse
if (epm eq '1') then epm = 'epm.dat'
chg_default_dir, epm, basedir

if (keyword_set(vpec_error) eq 0) then $
   vpec_error = double(input_param('vpec_error', 240.d0, cpars, cvals))

if (keyword_set(sbf_vpec_error) eq 0) then $
   sbf_vpec_error = double(input_param('sbf_vpec_error', vpec_error, $
                                       cpars, cvals))
; sbf_vpec_error defaults to the same value as the rest, but it can be
; set to a diffreent value.  Joe+John recommend a value of 275.

if (keyword_set(vcorr_sbf) eq 0) then $
   vcorr_sbf = input_param('vcorr_sbf', '2M++', cpars, cvals)
if (keyword_set(vcorr_mm) eq 0) then $
   vcorr_mm = input_param('vcorr_mm', '2M++', cpars, cvals)
if (keyword_set(vcorr_epm) eq 0) then $
   vcorr_epm = input_param('vcorr_epm', '2M++', cpars, cvals)
if (keyword_set(vcorr_sn1a) eq 0) then $
   vcorr_sn1a = input_param('vcorr_sn1a', '2M++', cpars, cvals)

; Control use of hosts without an associated calibrator.  By default
; they are included in the matrix solution, the chi^2, and the
; residual plots.  The match between hosts and calibrators is tracked
; by arrays HOST_CLASS and HOST_USED.  HOST_CLASS is a string array
; containing all calibrator categories used in the solution.  For now
; this can have elements 'sn1a', 'sn2', and 'tf'.  HOST_USED tracks
; which of the elements of the HOSTS array is used for each of the
; elements of host_class; it uses bitwise OR to set the lowest bit if
; the host is used for a calibrator in the first eelment of
; HOST_CLASS, the second lowest bit for the second element, and so
; on.  So if HOST_CLASS is [ 'sn1a', 'sn2','tf'], HOST_USED will be 1
; for a host matched only to one or more sn1a calibrators, 2 if it is
; matched to a sn2 calibrator, 4 if it matched to a TF calibrator, 3
; if both 1 and 2 are true, etc.
;

if (keyword_set(nocull) eq 0) then $
   nocull = input_param ('nocull', 0, cpars, cvals)
if (nocull eq 0) then do_cull = 1B else do_cull = 0B


; Control output and logging

if (keyword_set(plotfile) eq 0) then $
   plotfile = input_param('plotfile', 'none', cpars, cvals)
if (plotfile eq '' or strlowcase(plotfile) eq 'none') then $
   doplot = 0B else doplot = 1B

if (keyword_set(plabel) eq 0) then $
   plabel = input_param('plabel', '', cpars, cvals)

if (doprint) then begin
   if (strlowcase(logfile) eq 'screen') then begin
      openw, lun, '/dev/tty', /get_lun
   endif else begin
      if (print_append) then openw, lun, logfile, /get_lun, /append else $
         openw, lun, logfile, /get_lun
   endelse
endif

if (keyword_set(plabel) and doprint) then if (plabel[0] ne '') then $
   printf, lun, ' Run label: ', plabel

forward_function uniquify
if (do_debug) then stop

; First read all host data.  These are all the galaxies used as local
; distance calibrators, whether for SN1a, TF, or SN2.  (May incluide
; SBF as well at a later date; for now SBF are used in a slightly
; different way.)  Any system used as a local calibrator in the
; *_calib files that is not present in the hosts data will be ignored.
; The distance modulus for ALL systems in the hosts data will be
; fitted, but will only result in an H0 constraint if the system
; appears in a *_calib file.

readcol, host_data, host2, mu_host2_value, mu_host2_error, $
         methodname_2, anchorname_2, sourcename_2, $
         comment='#', format='a,f,f,a,a,a', /silent

if (doprint) then printf, lun, ' Hosts read from file ' + $
                          host_data, format='(a)'

; anchorname = 'LMC', 'MW', 'N4258'
; hostname = galaxy designation
; method = 'CEPH', 'TRGB', 'JAGB', etc

is_anchor = host2 eq 'N4258' or host2 eq 'MW' or $
            host2 eq 'LMC' or host2 eq 'SMC' or host2 eq 'M31'
wh = where (is_anchor)
whn = where (is_anchor eq 0)

mas_ref_value = mu_host2_value [wh]
mas_ref_error = mu_host2_error [wh]
mas_ref_anchor = anchorname_2 [wh]
mas_ref_method = methodname_2 [wh]
mas_ref_source = sourcename_2 [wh]

host2 = host2[whn]
mu_host2_value = mu_host2_value[whn]
mu_host2_error = mu_host2_error[whn]
methodname_2 = methodname_2 [whn]
anchorname_2 = anchorname_2 [whn]
sourcename_2 = sourcename_2 [whn]

; Here apply the inclusions and exclusions

; ===================== INCLUSION/EXCLUSION =========================

; Note: the python include/exclude list has a feature that allows
; specifying only match in a specific category, e.g., something
; like anchor:N4258 will match N4258 only when it apears as anchor.
; This can cut down on undesired matches.  This feature is not
; currently available in the IDL version.

; Also note: the exclusion list allows each exclude string to
; contain ampersand-separated substrings.  

; Final note: blanks are currently meaningful within an inclusion or
; exclusion string (not initial or final).  This behavior may change
; in future versions, depending on usage.

if (size(include_list, /type) eq 0) then include_list = ''
if (size(exclude_list, /type) eq 0) then exclude_list = ''
do_include = 0B
if (n_elements(include_list) gt 1) then do_include = 1B
if (include_list[0] ne '' and include_list[0] ne 'none') then do_include = 1B
do_exclude = 0B
if (n_elements(exclude_list) gt 1) then do_exclude = 1B
if (exclude_list[0] ne '' and exclude_list[0] ne 'none') then do_exclude = 1B

if (do_include and do_exclude) then begin
   if (doprint) then printf, lun, ' WARNING: include_list is set; ' + $
                             ' exclude_list will be reset to a blank'
   exclude_list = ''
   do_exclude = 0B
endif

if (do_debug) then stop

n2 = n_elements(host2)
wh = lindgen(n2)     ; accept all by default
nacc = n2            ; unless either include_list or exclude_list are set
if (do_include) then begin
   ; Handle inclusion first
   iincl = replicate(0B, n2)     ; set if item is to be included
   for k = 0, n_elements(include_list)-1 do begin
      if (include_list[k] ne '') then begin
         ; Next line needed to accept multiple inclusions 
         ; separated by '&'.  Recommended method is to use an array instead.
         this_incl = strtrim(strsplit(include_list[k], '&', /extract), 2)
         for jincl = 0, n_elements(this_incl)-1 do begin
            aa = strmatch ([[host2],[methodname_2],[anchorname_2],$
                      [sourcename_2]], this_incl[jincl], /fold_case)
            iincl = iincl or aa[*,0] or aa[*,1] or aa[*,2] or aa[*,3]
            ; include if any matches succeed
         endfor
      endif
   endfor
   wh = where (iincl, nacc)   ; this is the array of included items
endif else if (do_exclude) then begin
   ; Now exclusion
   iexcl = replicate (0B, n2)   ; set if item is to be excluded
   for k = 0, n_elements(exclude_list)-1 do begin
      if (exclude_list[k] ne '') then begin
         this_excl = strtrim(strsplit(exclude_list[k], '&', /extract),2)
         for j = 0, n_elements(this_excl)-1 do begin
            aa = strmatch ([[host2],[methodname_2],[anchorname_2],$
                      [sourcename_2]], this_excl[j], /fold_case)
            iexcl = iexcl or aa[*,0] or aa[*,1] or aa[*,2] or aa[*,3]
            ; exclude if any matches succeed
         endfor
      endif
   endfor
   wh = where (iexcl eq 0, nacc)
endif

if (nacc le 0) then begin
   if (doprint) then printf, lun, ' Too many items excluded, nothing left'
   if (do_debug) then stop
   goto, closeout
endif
if (do_debug) then stop

if (nacc lt n2) then begin
   if (doprint) then printf, lun, 'Excluding ', n2-nacc, ' host measurements'+$
                             ' because of include/exclude lists'
   host2 = host2[wh]
   mu_host2_value = mu_host2_value[wh]
   mu_host2_error = mu_host2_error[wh]
   methodname_2 = methodname_2 [wh]
   anchorname_2 = anchorname_2 [wh]
   sourcename_2 = sourcename_2 [wh]
   n2 = nacc
endif

if (do_debug) then stop

; Check whether the SN host is in the list of local calibrators
; (otherwise it cannot be used)

; Here is the place to verify all calib files and exclude hosts that
; are not used.  However, at this time we have NOT read in all the
; relevant calibration files.  We need to do that first.
; Relevant calibration files may include:
;    - SN1a
;    - SN2
;    - TF
;    - SBF
; We will split the two processes: read in calibrator data (and apply
; exclusions), and host matching.  The latter can result in removal of
; unused hosts.  Variable definition has to happen alter, since no
; variables are defined for hosts that are removed.
;

; Culling will be based on host2 array.  Create hosts array after the culling.


; Read SN1a calibrator file.  NO MATCHING TO HOSTS AT THIS STAGE (done
; after culling).  SN Inclusion/exclusion is applied here (before
; culling, of course).

n_sn1a = 0
if (do_sn1a) then begin
   readcol, sn1a_calib, host_sn1a, name_sn1a, app_sn1a, sigma_sn1a, $
            format='x,a,a,f,f', comment='#', /silent
   if (doprint) then printf, lun, ' SN Ia calibrators read from file ' + $
                             sn1a_calib, format='(a)'
   n_sn1a = n_elements(name_sn1a)
   if (sn1a_intrinsic gt 0) then sigma_sn1a = sqrt(sigma_sn1a^2 + $
                                                   sn1a_intrinsic^2)
   ; Apply SN exclusion
   ; =================== SN INCLUSION AND EXCLUSION =========================
   ; If SN_INCLUDE is not set, then all are included by default.
   ; SN_EXCLUDE overrides SN_INCLUDE.
   sn1a_use = bytarr(n_sn1a)
   ; isn = bytarr(n_sn1a)
   ; Inclusion (with wild card); defaults to include all
   if (n_elements(sn_include) eq 1 and sn_include[0] eq '') then begin
      sn1a_use[*] = 1B
   endif else begin
      for k = 0, n_elements(sn_include)-1 do $
         sn1a_use = sn1a_use or strmatch (name_sn1a, sn_include[k], /fold_case)
   endelse
   ; Exclusion (with wild card; overrides inclusion)
   if (n_elements(sn_exclude) gt 1 or sn_exclude[0] ne '') then begin
      for k = 0, n_elements(sn_exclude)-1 do $
         sn1a_use = sn1a_use and not strmatch(name_sn1a, sn_exclude[k], /fold_case)
   endif
   wh = where (sn1a_use, nwh)
   if (nwh lt n_sn1a) then begin
      if (doprint) then printf, lun, n_sn1a-nwh, ' SN Ia calibrators ' + $
                                'not included or excluded by request'
      host_sn1a = host_sn1a[wh]
      name_sn1a = name_sn1a[wh]
      app_sn1a  = app_sn1a[wh]
      sigma_sn1a = sigma_sn1a[wh]
      n_sn1a = nwh
   endif
endif
if (n_sn1a eq 0) then host_sn1a = ''

; Read sn2_calib file.  NO MATCHING TO HOSTS AT THIS STAGE (done
; after culling.  Inclusion/exclusion (same variables as for SN1a)
; is applied here, before culling.

n_sn2 = 0
if (do_sn2) then begin
   readcol, sn2_calib, host_sn2, name_sn2, app_sn2, sigma_sn2, $
            format='x,a,a,f,f', comment='#', /silent
   if (doprint) then printf, lun, ' Empirical SN II calibrators read ' + $
                             'from file ' + sn2_calib, format='(a)'
   n_sn2 = n_elements(name_sn2)
   ; Inclusion/Exclusion
   sn2_use = bytarr(n_sn2)
   ; Inclusion (with wild card); defaults to include all
   if (n_elements(sn_include) eq 1 and sn_include[0] eq '') then begin
      sn2_use[*] = 1B
   endif else begin
      for k = 0, n_elements(sn_include)-1 do $
         sn2_use = sn2_use or strmatch (name_sn2, sn_include[k], /fold_case)
   endelse
   ; Exclusion (with wild card; overrides inclusion)
   if (n_elements(sn_exclude) gt 1 or sn_exclude[0] ne '') then begin
      for k = 0, n_elements(sn_exclude)-1 do $
         sn2_use = sn2_use and not strmatch(name_sn2, sn_exclude[k], /fold_case)
   endif
   wh = where (sn2_use, nwh)
   if (nwh lt n_sn2) then begin
      if (doprint) then printf, lun, n_sn2-nwh, ' SN Ia calibrators ' + $
                                'not included or excluded by request'
      host_sn2 = host_sn2[wh]
      name_sn2 = name_sn2[wh]
      app_sn2  = app_sn2[wh]
      sigma_sn2 = sigma_sn2[wh]
      n_sn2 = nwh
   endif
endif
if (n_sn2 eq 0) then host_sn2 = ''

; Read sbf_calib file.  NO MATCHING TO HOSTS.
; NO GROUP CONSIDERATIONS (yet).

n_sbf = 0
if (do_sbf) then begin
   readcol, sbf_calib, name_sbf, big_m110_sbf, sigma_big_m110_sbf, $
            small_m110_sbf, sigma_small_m110_sbf, $
             format = 'a,f,f,f,f', comment='#', /silent
   n_sbf = n_elements(name_sbf)
   if (doprint) then if (n_sbf gt 0) then $
      printf, lun, n_sbf, ' SBF calibrators read ' + $
              'from file ' + sbf_calib, format='(2x, i6, a)'
endif
if (n_sbf eq 0) then name_sbf = ''

; Read tf_calib file.  NO MATCHING.

n_tf = 0
if (do_tf) then begin
   readcol, tf_calib, name_tf, small_m_tf, big_m_tf, sigma_big_m_tf, $
            format='a,x,f,f,f', comment='#', /silent
   n_tf = n_elements(name_tf)
   if (doprint) then if (n_tf gt 0) then $
      printf, lun, n_tf, ' Tully-Fisher calibrators read ' + $
              'from file ' + tf_calib, format='(2x, i6, a)'
endif
if (n_tf eq 0) then name_tf = ''

; Read groups file.  No checking or matching.
nghost = 0
if (do_groups) then begin
   readcol, groups_file, gname, g_host_name, g_host_sigma, $
            format='a,a,f', comment='#', /silent
   nghost = n_elements(g_host_name)
   if (doprint) then if (nghost gt 0) then $
      printf, lun, nghost, ' Group members read ', $
              format = '(2x, i6, a)'
endif
if (nghost eq 0) then g_host_name = ''

; Now cull hosts unless nocull is set.  All checks are done on the
; host2 array, and the hosts array is generated again after that.

; Create initial hosts array.  Will be redefined after culling.

n2 = n_elements(host2)
host2_used = lonarr(n2)
host2_used = long64(host2_used)

; Create lists to simplify the code
numbers = [n_sn1a, n_sn2, n_sbf, n_tf, nghost]
h_arrays = list (host_sn1a, host_sn2, name_sbf, name_tf, g_host_name)
bitvalues = [1, 2, 3, 4, 5]
; Find out what host2 are used for
for jcat = 0, n_elements(numbers)-1 do begin
   if (numbers[jcat] gt 0) then begin
      hnames = h_arrays[jcat]
      ivv = 2LL^bitvalues[jcat]
      for k = 0, numbers[jcat]-1 do host2_used = host2_used or $
         (ivv * (host2 eq hnames[k]))
   endif
endfor

if (do_debug) then stop

; If requested, cull host2
if (do_cull) then begin
   wh = where (host2_used ne 0, nwh) ; Hosts that are used for at least one calibration
   if (nwh gt 0) then begin
      if (nwh lt n2) then begin
         host2 = host2[wh]
         host2_used = host2_used [wh]
         mu_host2_value = mu_host2_value[wh]
         mu_host2_error = mu_host2_error[wh]
         methodname_2 = methodname_2 [wh]
         anchorname_2 = anchorname_2 [wh]
         sourcename_2 = sourcename_2 [wh]
         if (doprint) then printf, lun, n2-nwh, ' unused hosts culled', format='(2x,i6,a)'
      endif
   endif else begin
      if (doprint) then printf, lun, ' No hosts remain; abort'  ; should this be legal?
      if (do_debug) then stop
      goto, closeout
   endelse
endif

if (do_debug) then stop

n2 = n_elements(host2)
hosts = uniquify (host2, host_index_2)
nhosts = n_elements(hosts)
host_used = lonarr(nhosts)
host_used = long64(host_used)
for k = 0, n2-1 do host_used[host_index_2[k]] = host_used[host_index_2[k]] or host2_used[k] 

; Note: hosts applies to both host_data and calibrator arrays.
; methods and anchors only apply to host_data.

; uniquify with two arguments returns the index of the first argument
; into the result.  Thus by construction,
; aa = uniquify (bb, cc)
; is defined so that bb equals aa[cc].
;

; The first NHOSTS variables are the distances to hosts in the hosts array.

npars = n_elements (hosts)
par0 = {name:'', value:0.d0, fixed:0b, limits:[0.d0,0.d0], limited:[0B,0B], $
        type:0, objname:''}
parinfo = replicate (par0, npars)
for k = 0, npars-1 do begin
   parinfo[k].name = 'Distance modulus for host/calibrator '+hosts[k]
   parinfo[k].objname = hosts[k]
endfor
ipar = npars-1

; Add groups to hosts.  If do_groups (which means groups is not 'none'
; or ''), then read the groups file, which contains three columns:
; group name, host name, and g_host_sigma.  (In principle each host can
; have a different distance modulus dispersion with respect to the group.)
; hosts that are not in the host_data file are added, but will not
; have values for quantities such as anchor and mas.  We will need to
; check if that creates issues elsewhere; probably not, as these
; quantities are used in setting up distance network equations, and
; these equations do not exist for such hosts.  Hosts that are used only
; for groups are cnsidered used; it is the user's
; responsibility not to add a groups file that is not used in the
; network.
; Additional arrays: igh is the index of each host in groups in the
; overall hosts array; igg is the corresponding group for that host in
; the groups array.  The array igp is the distance modulus parameter
; for each group.

if (do_groups) then begin
   ; Already read
   ; readcol, groups_file, gname, g_host_name, g_host_sigma, format='a,a,f', $
   ;         comment='#', /silent
   groups = uniquify(gname, whg)
   nghost = n_elements(g_host_name)
   igh = replicate (-1, nghost)
   igg = whg
   for k = 0, nghost-1 do begin
      ; Check if the host is already in hosts
      whg = where (hosts eq g_host_name[k], nwhg)
      if (nwhg eq 1) then begin
         igh[k] = whg[0] ; index of matching host
      endif else begin
         if (nwhg gt 1) then begin
            print, ' Found multiple hosts match to host in groups; should not happen'
            stop
         endif
         ; New entry in hosts.  Also add parameter and host_used element
         hosts = [hosts, g_host_name[k]]
         nhosts += 1
         host_used = [host_used, 2LL^5]
         parinfo = [parinfo, par0]
         ipar += 1
         igh[k] = ipar
         parinfo[ipar].name = 'Distance modulus for host/calibrator '+g_host_name[k]
         parinfo[ipar].objname = g_host_name[k]
         npars += 1
      endelse
   endfor
                                ; Now add parameter for each group
   igp = lonarr(n_elements(groups))
   for k = 0, n_elements(groups)-1 do begin
      ipar += 1
      npars += 1
      parinfo = [parinfo, par0]
      parinfo[ipar].name = ' Distance modulus for group '+groups[k]
      parinfo[ipar].objname = groups[k]
      igp[k] = ipar
   endfor
endif

; Will add the relevant equations just after setting up all hosts equations.
;

; At this point there are nhosts + ngroups optimization parameters.
   
; Now decide whether sn1a, sn2, tf, sbf are used of not, and set up
; the appropriate parameters for each.  Also match to hosts array as
; needed.  This was originally done immediately after reading the
; calibration files, but now it needs to be done separately: read the
; calibrator information, cull the hosts, then match calibrators to
; hosts.

; Whether these categories are used or not is based on the variables
; do_sn1a, do_sn2, do_sbf, do_tf.
; For each category used, there is one parameter and n+1 equations.
; The parameter is the reference magnitude calibration or its offset
; from the nominal value..
; The equations are the host distance equations and the Hubble flow
; calibration.
; Equations are all set later.

iabs_sn1a = -1
if (do_sn1a) then begin
   npars = npars+1
   parinfo = [parinfo, par0]
   iabs_sn1a = npars-1
   parinfo[iabs_sn1a].name = 'Reference absolute magnitude of SNe Ia'
endif

iabs_sn2 = -1
if (do_sn2) then begin
   npars = npars+1
   parinfo = [parinfo, par0]
   iabs_sn2 = npars-1
   parinfo[iabs_sn2].name = 'Reference absolute magnitude of SNe II'
endif

iabs_sbf = -1
if (do_sbf) then begin
   npars = npars+1
   parinfo = [parinfo, par0]
   iabs_sbf = npars-1
   parinfo[iabs_sbf].name = 'SBF absolute magnitude offset (m110=M110+M110_offset+mu) '
endif

iabs_tf = -1
if (do_tf) then begin
   npars = npars+1
   parinfo = [parinfo, par0]
   iabs_tf = npars-1
   parinfo[iabs_tf].name = 'TF absolute magnitude offset (m=M+M_offset+mu) '
endif

; And finally the parameter for the logarithm of the Hubble constant

npars = npars+1
parinfo = [parinfo, par0]
ihub = npars-1
parinfo[ihub].name = 'log_10 (Hubble constant) in km/s/Mpc'

; This sets up most parameters.  An additional parameter may be set
; for Coma.  Group parameters for SBF have already been set.


;=========================================================
;===             START UNPACKING HOST DATA             ===
;=========================================================

; Host_Data: for each line in the data file we have an equation
; relating the anchor distance to the distance of that host.
; Anchor distances are given first here.

; And then we have the zero points for each method and anchor
; "mas" refers to method + anchor + source (the assumption is that
; method+anchor+surce errors are uncorrelated if any is different)

; readcol, 'methods.dat', ma_anchor_name, ma_method_name, ma_ref_value, $
;          ma_ref_error, format = 'a,a,f,f'

; readcol, 'sources.dat', etc etc

methods = uniquify (methodname_2, method_index)
anchors = uniquify (anchorname_2, anchor_index)
sources = uniquify (sourcename_2, source_index)

; ad-hoc for anchor

n_anchors = n_elements(anchors)
; Fill in geonetric uncertainty for anchors.  It assumes that anchors
; are only N4258, MW, LMC, or SMC.  Assume the distance modulus
; uncertainty is 0.03.

mu_anchor_error = dblarr (n_anchors)


if (anchor_data ne '' and anchor_data ne 'none') then begin
   readcol, anchor_data, data_anchors, data_mu_value, $
                    data_mu_sigma, format = 'a,d,d', /silent
   if (doprint) then printf, lun, ' Anchor errors read from file '+$
                             'anchors.dat', format='(a)'
endif else begin
   data_anchors = ['N4258','MW','LMC','SMC']
   data_mu_value = [0,0,0,0]
   data_mu_sigma = [0.032d0, 0.020d0, 0.024d0, 0.032d0]
endelse

; data_mu_sigma here is only the geometric component


for k = 0, n_anchors-1 do begin
   wh = where (data_anchors eq anchors[k], nwh)
   if (nwh gt 0) then mu_anchor_error[k] = data_mu_sigma[wh[0]] else $
      mu_anchor_error[k] = 1.d0   ; very large error, no information
endfor

; Note: the lines corresponding to anchor and method should be stripped
; from the host_data input.  This needs to happen at some point.  Also we
; can check if the assumed distance modulus matches that in
; anchor_data, and if not we can adjust all relevant distances accordingly.
; Likely future improvement; currently the distances in the anchors
; file are not used.

mas_name = methodname_2 + '&' + anchorname_2 + '&' + sourcename_2
mas = uniquify (mas_name, mas_index)
n_mas = n_elements(mas)
mas_error = dblarr (n_mas)
mas_ref_name = mas_ref_method + '&' + mas_ref_anchor + '&' + mas_ref_source

for k = 0, n_mas-1 do begin
   wh = where (mas_ref_name eq mas[k], nwh)
   if (nwh gt 0) then mas_error[k] = mas_ref_error[wh[0]] else begin
      if (doprint) then begin
         printf, lun, ' No reference error found for '+mas_name[k]
         printf, lun, 'Bailing out'
      endif
      goto, closeout
   endelse
   if (nwh gt 1 and doprint) then begin
      printf, lun, ' Warning: multiple matches found for ', mas_name[k]
      printf, lun, ' First value taken'
   endif
endfor

; Split the part of the code that reads calibrators from the matching
; to host names.  The reason is to enable exclusion of hosts that do
; not appear in the calibration files.  

; Now match the calibrators for SN1a, SN2, SBF, and TF - if used - to
; the hosts in the hosts array.  The calibrators and their hosts 
; have been read before.  (Check consistency?)

; SN1a calibrators already read, inclusions/exclusions applied
; Match to hosts
; Remove if not matched

if (do_sn1a and n_sn1a gt 0) then begin
   ; Match to a HOST element
   ihost_sn1a = replicate (-1, n_sn1a)
   for k = 0, n_sn1a-1 do begin
      wh = where (hosts eq host_sn1a[k], nwh)
      if (nwh eq 1) then ihost_sn1a[k] = wh[0]
   endfor
   wh = where (ihost_sn1a ge 0, nwh)
   if (n_sn1a gt nwh) then begin
      if (doprint) then begin
         printf, lun, n_sn1a-nwh, ' SN Ia calibrators removed because ' + $
                 'host is not in the HOSTS file'
         host_sn1a = host_sn1a[wh]
         name_sn1a = name_sn1a[wh]
         app_sn1a = app_sn1a[wh]
         sigma_sn1a = sigma_sn1a[wh]
         ihost_sn1a = ihost_sn1a[wh]
         n_sn1a = nwh
      endif
      if (nwh eq 0) then begin
         printf, lun, ' Error: SNe Ia requested, but no host matches. ', $
                 ' Optimization will fail.'
         do_sn1a = 0
         ; need code to avoid degeneracy
      endif
   endif
endif

; SN2

if (do_sn2 and n_sn2 gt 0) then begin
   ; Match to a HOST element
   ihost_sn2 = replicate (-1, n_sn2)
   for k = 0, n_sn2-1 do begin
      wh = where (hosts eq host_sn2[k], nwh)
      if (nwh eq 1) then ihost_sn2[k] = wh[0]
   endfor
   wh = where (ihost_sn2 ge 0, nwh)
   if (n_sn2 gt nwh) then begin
      if (doprint) then begin
         printf, lun, n_sn2-nwh, ' SN II calibrators removed because '+$
                 'host is not in the HOSTS file'
         host_sn2 = host_sn2[wh]
         name_sn2 = name_sn2[wh]
         app_sn2 = app_sn2[wh]
         sigma_sn2 = sigma_sn2[wh]
         ihost_sn2 = ihost_sn2[wh]
         n_sn2 = nwh
      endif
      if (nwh eq 0) then begin
         printf, lun, ' Error: SNe II requested, but no host matches. ', $
                 ' Optimization will fail.'
         do_sn2 = 0
      endif
   endif
endif

; SBF
; Uses a calibrator file with host name, calibrator
; name, m110, error, M110, error.  M110 is based
; on the current calibration (Jensen+ 2025); as for Tully-Fisher, the
; fitted parameter is the calibration offset.  An earlier version
; used the published calibration in either Blakeslee 2021 or Jensen
; 2025, with special code to handle covariances.  

if (do_sbf and n_sbf gt 0) then begin
   host_sbf = name_sbf
   ihost_sbf = replicate (-1, n_sbf)
   for k = 0, n_sbf-1 do begin
      wh = where (hosts eq host_sbf[k], nwh)
      if (nwh eq 1) then ihost_sbf[k] = wh[0]
   endfor
   wh = where (ihost_sbf ge 0, nwh)
   if (n_sbf gt nwh) then begin
      if (doprint) then begin
         printf, lun, n_sbf-nwh, ' SBF calibrators removed because '+$
                 'host is not in the HOSTS or GROUPS file'
         host_sbf = host_sbf[wh]
         name_sbf = name_sbf[wh]
         ihost_sbf = ihost_sbf[wh]
         small_m110_sbf = small_m110_sbf [wh]
         sigma_small_m110_sbf = sigma_small_m110_sbf [wh]
         big_m110_sbf = big_m110_sbf [wh]
         sigma_big_m110_sbf = sigma_big_m110_sbf [wh]
         n_sbf = nwh
      endif
      if (nwh eq 0) then begin
         printf, lun, ' Error: SBF requested, but no host matches. ', $
                 ' Optimization will fail.'
         do_sbf = 0
      endif
   endif
endif

; TF

if (do_tf and n_tf gt 0) then begin
   host_tf = name_tf
   ihost_tf = replicate (-1, n_tf)
   for k = 0, n_tf-1 do begin
      wh = where (hosts eq host_tf[k], nwh)
      if (nwh eq 1) then ihost_tf[k] = wh[0]
   endfor
   wh = where (ihost_tf ge 0, nwh)
   if (n_tf gt nwh) then begin
      if (doprint) then begin
         printf, lun, n_tf-nwh, ' TF calibrators removed because '+$
                 'host is not in the HOSTS file'
         host_tf = host_tf[wh]
         name_tf = name_tf[wh]
         ihost_tf = ihost_tf[wh]
         small_m_tf = small_m_tf [wh]
         big_m_tf = big_m_tf [wh]
         sigma_big_m_tf = sigma_big_m_tf [wh]
         n_tf = nwh
      endif
      if (nwh eq 0) then begin
         printf, lun, ' Error: TF requested, but no host matches. ', $
                 ' Optimization will fail.'
         do_tf = 0
      endif
   endif
endif

; If coma is given read the file with coma SN information and add the
; relevant variables and equations
; Coma does not include host galaxies, only SNe Ia to be calibrated
; Currently these SNe Ia are not subject to include/exclude rules.

n_coma = 0
do_coma = 0
if (coma_found and coma ne 'none' and coma ne '') then do_coma = 1B
if (do_coma) then begin
   if (strmid(coma,0,1) ne '/') then coma_file = coma else $
      coma_file = coma
   readcol, coma_file, coma_sn_name, coma_sn_mag, coma_sn_err, $
            format='a,f,f', comment='#', /silent
   if (doprint) then printf, lun, ' Coma SN Ia calibrators read ' + $
                             'from file ' + coma_file, format='(a)'
   n_coma = n_elements(coma_sn_mag)
   npars = npars+1
   parinfo = [parinfo, par0]
   icoma = npars-1
   parinfo[icoma].name = 'Distance modulus to the Coma cluster'
endif

; Add direct constraints from megamasers
; these have the form:
; Name   pred_value   error
; where the errors are uncorrelated
; We could change to have distance and redshift directly; then we have
; to assume a kinematic model to correct c_z.  Errors must include
; pculiar velocity contributions.  There is a potential problem in
; that peculiar velocity errors can be correlated.  Ignore these
; details for now.

; Note that direct constraints do not add to the number of parameters,
; just equations.

; Read Megamaser data if requested
n_mm = 0
do_mm = 0B
if (mm_found and mm ne 'none' and mm ne '') then do_mm = 1B
mm_file=mm
if (do_mm) then begin
   if (keyword_set(vcorr_mm) eq 0) then vcorr_mm = '2M++'
   fmt = 'a,f,f,f,f,f,x,f'
   case strupcase(vcorr_mm) of 
      'GROUP': fmt = 'a,f,f,f,f,f,f'
      '2M++':  fmt = 'a,f,f,f,f,f,x,f'
      'CF3':   fmt = 'a,f,f,f,f,f,x,x,f'
      'M2000': fmt = 'a,f,f,f,f,f,x,x,f'
      else: if (doprint) then printf, lun, 'Unrecognized VCORR_MM ' + $
                                      'option for megamaser ; using 2M++'
   endcase
   readcol, mm_file, mm_name, mm_distance, $
            mm_d_minus, mm_d_plus, mm_v_obs, mm_v_error, $
            mm_v_corr, format=fmt, comment='#', /silent
   if (doprint) then printf, lun, ' Megamaser information read ' + $
                             'from file ' + mm_file, format='(a)'
   mm_distance_error = 0.5*(abs(mm_d_plus) + abs(mm_d_minus))
   n_mm = n_elements(mm_name)
endif

; Read SN II EPM data if requested; append to megamaser if given.
n_epm = 0
do_epm = 0B
if (epm_found and epm ne 'none' and epm ne '') then do_epm = 1B
epm_file=epm
if (do_epm) then begin
   if (keyword_set(vcorr_epm) eq 0) then vcorr_epm = '2M++'
   fmt = 'x,x,a,f,f,f,x,f'
   case strupcase(vcorr_epm) of 
      'CMB':   fmt = 'x,x,a,f,f,f,x,f'
      '2M++':  fmt = 'x,x,a,f,f,f,x,f'
      else: if (doprint) then printf, lun, 'Unrecognized VCORR_EPM ' + $
                                      'option for SNe II ; using 2M++'
   endcase
   readcol, epm_file, epm_name, epm_distance, $
            epm_distance_error, epm_z_cmb, epm_z_corr, $
            format=fmt, comment='#', /silent
   if (doprint) then printf, lun, ' Astrophysically-calibrated SN II read ' + $
                             'from file ' + epm_file, format='(a)'
   epm_v_corr = epm_z_corr*(!const.c/1.d3)
   epm_v_error = epm_v_corr*0 + 1.d0  ; nominal velocity error
   epm_v_cmb =  epm_z_cmb*(!const.c/1.d3)
   if (strupcase(vcorr_epm) eq 'CMB') then epm_v_corr = epm_v_cmb
   n_epm = n_elements(epm_name)
endif

; ================================================
; ===            EQUATION CREATION             ===
; ================================================

; Create the equations and populate the covariance matrix

; host_data creates n2 equations
; One equation per host_data entry to tie host to anchor
; SN1a: one equation per calibrator + one to tie to H0
; SN2: one equation per calibrator plus one to tie to H0
; TF: One equation per calibrator plus one to tie to H0
; Coma/FP: One equation per calibrator + one for SBF constaint + 1 to
;    tie to H0
; MM: one equation per object tied directly to H0
; EPM: one equation per object tied directly to H0
; SBF: one equation per HF object tied directly to H0

neq = n2
ieq_host_start = 0
ieq_host_end = n2-1

; Set up equation numbers for groups if given
ieq_group_start = -1
ieq_group_end = -1
if (do_groups) then begin
   ieq_group_start = neq
   neq += n_elements(igh)
   ieq_group_end = neq-1
endif

; Set up equation numbers for all sets of data
; SN1a: N_SN1a + 1 (+1 for the Huble flow equation)
ieq_sn1a_start = -1
ieq_sn1a_end = -1
if (do_sn1a) then begin
   ieq_sn1a_start = neq
   neq += n_sn1a + 1
   ieq_sn1a_end = neq-1
endif

ieq_sn2_start = -1
ieq_sn2_end = -1
if (do_sn2) then begin
   ieq_sn2_start = neq
   neq += n_sn2 + 1
   ieq_sn2_end = neq-1
endif

ieq_sbf_start = -1
ieq_sbf_end = -1
if (do_sbf) then begin
   ieq_sbf_start = neq
   neq += n_sbf + 1
   ieq_sbf_end = neq-1
endif

ieq_tf_start = -1
ieq_tf_end = -1
if (do_tf) then begin
   ieq_tf_start = neq
   neq += n_tf + 1
   ieq_tf_end = neq-1
endif

ieq_mm_start = -1
ieq_mm_end = -1
if (do_mm) then begin
   ieq_mm_start = neq
   neq += n_mm
   ieq_mm_end = neq-1
endif

ieq_epm_start = -1
ieq_epm_end = -1
if (do_epm) then begin
   ieq_epm_start = neq
   neq += n_epm
   ieq_epm_end = neq-1
endif

ieq_coma_start = -1
ieq_coma_end = -1
if (do_coma) then begin
   ieq_coma_start = neq
   neq += n_coma + 2
   ieq_coma_end = neq-1
endif

; The two extra equations for Coma are the SBF distance constraint
; and the tie between Coma and H0

; Note that there are typically n_calib + 1 equations for each group.
; The n_calib are the equations that match the absolute calibration
; of the distance indicator to the distance to each calibrator
; The +1 is the Hubble flow equation that uses alpha.  There is none
; for direct distance indicators (mm, epm); there are 2 for Coma
; (includes the sbf constraint).
;

covar = dblarr (neq, neq)
coeffs = dblarr (neq, npars)
yval = dblarr (neq)
yvar = dblarr (neq)
eq_descr = strarr(neq)
eq_shape = strarr(neq)
eq_covar = strarr(neq)

; EQUATIONS: HOSTS

; The arrays below keep track of the anchor and mas error for each
; host_data entry.  These are the terms that get added to all covariances
; with brethren.

mas_host_error = dblarr(n2)
anchor_dist_error = dblarr(n2)

; Add all equations for host_data

for k = 0, n2-1 do begin
   ieq = ieq_host_start + k
   ihost = host_index_2[k]
   coeffs[ieq, ihost] = 1.d0
   ; This equation has the form: 
   ;    mu_host = mu_anchor + (mu_host-mu_anchor) + error
   ; and
   ;    mu_host-mu_anchor = ma_host_method_refvalue - ma_anchor_method_refvalue
   ; with reference appropriate uncertainties.
   ; We subsume all terms into a single constant; no need to track
   ; individual components, although they can be decomposd afterwards.
   yval[ieq] = mu_host2_value[k]
   yvar[ieq] = mu_host2_error[k]^2
   covar[ieq, ieq] = mu_host2_error[k]^2
   eq_descr[ieq] = ' Host equation for ' + host2[k] + $
                   ' with MAS: ' + mas_name[k]
   eq_shape[ieq] = ' mu_host = reported value'
   eq_covar[ieq] += 'Variance of host measurement'
   ; Look at the covariances and ensure they are reasonable.
   ; The total uncertainty of this distance estimate includes:
   ;   - Anchor distance uncertainty
   ;   - Anchor reference value uncertainty (depends on method)
   ;   - Host reference value uncertainty (depends on method)
   ; The last term should be what is left after subtracting
   ; the other two from the total error in quadrature.
   ; So
   y_res_var = mu_host2_error[k]^2 ; - mu_anchor_error[anchor_index[k]]^2 ;
              ; - mas_error[mas_index[k]]^2
              ; This term is not included in mu_host2_error; see below.
   ; Assume that the entries in host_data.dat do NOT include the anchor 
   ; geometric distance uncertainty, and add it to the diagonal
   ; below.
   if (y_res_var le 0) then begin
      if (doprint) then printf, lun, ' Unacceptable value of residual ' + $
                                ' variance for entry ', k
      if (doprint) then printf, lun, ' Bailing out'
      goto, closeout
   endif
endfor

; Now add in all the appropriate covariances between host_data equations
; Every term using the same anchor will have the anchor distance
; uncertainty in the appropriate rows of the covariance matrix

; IMPORTANT NOTE: the uncertainties listed for each host in the
; host_data.dat file DO NOT include the MAS reference uncertainty, which
; is given as the NGC 4258 entry for the same MAS.  So the
; covariance matrix needs to add that term explicitly on the diagonal
; as well.

for k = 0, n_elements(anchors)-1 do begin
   wh = where (anchor_index eq k, nwh)
   if (nwh gt 0) then anchor_dist_error[wh] = mu_anchor_error[k]
   if (nwh gt 0) then for i = 0, nwh-1 do for j = 0, nwh-1 do $
      covar[wh[i],wh[j]] += mu_anchor_error[k]^2
   if (nwh gt 0) then eq_covar [wh] += ' + geometric anchor ' + $
                                       anchors[k]+' distance'
   if (do_debug) then stop
endfor

; Every term using the same method+anchor+source combination will share the
; same reference value uncertainty

for k = 0, n_mas-1 do begin
   wh = where (mas_index eq k, nwh)
   if (nwh gt 0) then mas_host_error[wh] = mas_error[k]
   if (nwh gt 0) then for i = 0, nwh-1 do for j = 0, nwh-1 do $
      covar[wh[i],wh[j]] += mas_error[k]^2
   if (nwh gt 0) then eq_covar[wh] += ' + MAS '+mas[k]+' variance'
   if (do_debug) then stop
endfor

; ====== covariance for same host, method and source

; Must add to covariance off-diagonal for hosts with multiple entries with
; same method and source, but different anchor.  In this case the
; covariance value should be the error in host_data with the reference error
; subtracted in quadrature, BUT we know that that's inconsistent.

; Assume that host_data files will be produced that do NOT include
; the reference error in each host entry.

for k = 0, n_elements(hosts)-1 do begin   ; go over all possible hosts
   wh = where (host2 eq hosts[k], nwh)    ; all equations for that host
   if (nwh gt 1) then begin               ; if only one entry, do nothing
      method_source = methodname_2[wh]+'&'+sourcename_2[wh]
                                          ; create method_source array
      ms_values = uniquify(method_source) ; find possible values
      nms = n_elements(ms_values)
      for ims = 0, nms-1 do begin
         whms = where (method_source eq ms_values[ims], nwhms)
         if (nwhms gt 1) then begin
            orig_idx = wh[whms]   ; original index of elements
            common_ms_covar = min(mu_host2_error[orig_idx]^2)
            ; These SHOULD be the same, but they were not in the
            ;    earlier version of host_data.dat.  They have now been 
            ;    homogeneized.
            ; Ensure that the no diagonal term is smaller than 
            ;    the off-diagonal
            for ii = 0, nwhms-1 do for jj = 0, nwhms-1 do $
               if (ii ne jj) then $
                  covar[orig_idx[ii],orig_idx[jj]] += common_ms_covar
            eq_covar[orig_idx] += ' + method/source host covariance '
         endif
      endfor
   endif
endfor

; EQUATIONS: GROUPS
; There is an equation for each host-group pair.  The equation says
; that a certain host is in a certain group (ie.e, has the dame
; distance modulus) with some tolerance.  These equations are not
; covariant with any others.  (The parameters are affected by other
; equations, as normal.)
; The equations are driven by the arrays igh and igg, set up when
; groups_file is read.  Both arrays have length equal to the number of
; hosts in groups; igh is the host number, igg is the corresponding
; group number.  
; The equation has the form:
; mu_host - mu_group = 0 with variance = sigma_host^2

if (do_groups) then begin
   for k = 0, n_elements(igh)-1 do begin
      ieq = ieq_group_start + k
      ihost = igh[k]
      igroup = igg[k]
      imugroup = igp[igroup]  ; For two groups: igg = 0 or 1; igp is ipar of the corresponding mu_group
      coeffs [ieq, ihost] = 1.d0
      coeffs [ieq, imugroup] = -1.d0
      eq_descr [ieq] = ' Distance modulus match between host ' + hosts[ihost] + ' and group ' + groups[igroup]
      eq_shape [ieq] = ' mu_host - mu_group = 0'
      eq_covar [ieq] = 'None'
      yval [ieq] = 0.d0
      yvar [ieq] = g_host_sigma[k]^2
      covar[ieq,ieq] = g_host_sigma[k]^2
   endfor
endif
      
; Equations for all four (possible) secondary indicators: SN1a, SN2,
; SBF, TF.
; There is an equation for each calibrator, matching it to its
; host/galaxy, as well as one equation per indicator tying it to the
; Hubble constant via the Hubble flow value.

; SN1a: the equations for calibrators have the form:
;    abs_sn1a = app_sn1a - mu_host
; where mu_host is referred to the host for that NS1a, as matched
; above.  (SN1a without a matching host have been excluded.)
; Note that if n_sn1a is zero, then do_sn1a is also false, and
; abs_sn1a is not a parameter.  (all this can be checked in principle)
; iabs_sn1a is the parameter number for abs_sn1a.
; ihost_sn1a is the index of the host for that supernova into the
; hosts array, and thus also the parameter number for mu_host.

;
; Use abs_sn as additional parameter.
; abs_sn is the reference absolute magnitude of a normalized Type Ia supernova.
; In these units:
; abs_sn = sn_app_mag - mu_host
; logh = 0.2*abs_sn + a_b + 5
; iabs_sn1a is the parameter index for abs_sn1a
; ihub is the parameter index for the Hubble constant
;

if (do_sn1a) then begin
   for k = 0, n_sn1a-1 do begin
      ieq = k + ieq_sn1a_start    ; this is the equation tying the SN
                    ; apparent magnitude and the host distance
                    ; to the sn absolute magnitude parameter
      ihost = ihost_sn1a[k]    ; the index of the host -> mu_host
                    ; k is the SN index
      coeffs (ieq, ihost) = 1.d0
      coeffs (ieq, iabs_sn1a) = 1.d0
      yval[ieq] = app_sn1a[k]
      yvar[ieq] = sigma_sn1a[k]^2
      covar[ieq,ieq] = sigma_sn1a[k]^2
      ; app_sn1a = abs_1a + mu_host
      ; where app_sn1a is yval, abs_1a and mu_host are parameters
      ; and sigma_sn1a is the measurement error, which
      ; includes the intrinsic width of the SN distribution
      eq_descr[ieq] = ' Distance for SN ' + name_sn1a[k] + ' in host ' + $
                      hosts[ihost_sn1a[k]]
      eq_shape[ieq] = ' M0 + mu_host = sn_mag'
      eq_covar[ieq] = ' Variance of SN magnitude (incl intrinsic dispersion)'
   endfor

   ; Add the equation relating M0_1a to the Hubble constant
   ; Use a_b_sn1a (given value or computed from sn1a_hf)
   c = !const.c * 1.e-3  ; km/s
   if (read_sn1a_hf) then begin
      if (sn1a_covar ne 'none' and sn1a_covar ne '') then begin
         has_covar = 1B
         ; Read covariance matrix; do not read magnitude error
         readcol, sn1a_covar, value, format='f', /silent ; one-column file
         if (doprint) then printf, lun, ' SN Ia covariance matrix read ' + $
                             'from file ' + sn1a_covar, format='(a)'
         n_covar = nint(value[0])
         lcovar = reform (value[1:*], n_covar, n_covar)
         ; Read HF data in Nils's format
         readcol, sn1a_hf, name_sn1a_hf, app_sn1a_hf, err_sn1a_hf, $
                  zhel_sn1a_hf, zcmb_sn1a_hf, vpec1_sn1a_hf, $
                  vpec2_sn1a_hf, vpec3_sn1a_hf, vpec4_sn1a_hf, $
                  comment='#', format='a,f,f,f,f,f,f,f,f', /silent
         if (doprint) then printf, lun, ' SN Ia in Hubble Flow read ' + $
                                   'from file ' + sn1a_hf, format='(a)'
         case strupcase(vcorr_sn1a) of
            '2M++': vpec_sn1a_hf = vpec4_sn1a_hf ; this is vp_2mpp
            'CMB': vpec_sn1a_hf = replicate (0.d0, n_covar)
            '2MRS': vpec_sn1a_hf = vpec3_sn1a_hf ; this is vp_2mrs
            'SDSS': vpec_sn1a_hf = vpec2_sn1a_hf ; this is vp_2mpp_sdss_6df
            '6DF': vpec_sn1a_hf = vpec2_sn1a_hf  ; this is vp_2mpp_sdss_6df
            'PLAIN': vpec_sn1a_hf = vpec1_sn1a_hf ; this is vp
            else: begin
               if (doprint) then begin
                  printf, lun, ' Unrecognized vcorr_sn1a: ', vcorr_sn1a
                  printf, lun, ' Using 2M++'
                  vpec_sn1a_hf = vpec4_sn1a_hf
               endif
            end
         endcase
         vhel_sn1a_hf = c * ((1+zhel_sn1a_hf)^2-1) / ((1+zhel_sn1a_hf)^2+1)
         vcmb_sn1a_hf = c * ((1+zcmb_sn1a_hf)^2-1) / ((1+zcmb_sn1a_hf)^2+1)
         alpha_sn1a_value = compute_alpha (app_sn1a_hf, vpec=vpec_sn1a_hf, $
            vcmb=vcmb_sn1a_hf, vhel=vhel_sn1a_hf, q0=q0, j0=j0, $
            vdisp=vpec_error, covar=lcovar, details=alphadet_sn1a, $
            ignore_offdiag=sn1a_ignore_offdiag, $
            redshift_range=redshift_range, alpha_err=alpha_sn1a_error, $
            chisq=chisq_sn1a_hf, ndof=ndof_sn1a_hf, fail=fail)
         printf, lun, 'Alpha from SN1a:', alpha_sn1a_value, alpha_sn1a_error, $
                 chisq_sn1a_hf, ndof_sn1a_hf, format='(a,3f12.6,i6)'
      endif else begin
         has_covar = 0B
         ; Read information using the IR format from Lluis
         readcol, sn1a_hf, name_sn1a_hf, app_sn1a_hf, err_sn1a_hf, $
                  zhel_sn1a_hf, zcmb_sn1a_hf, zcorr_sn1a_hf, $
                  format='a,x,f,f,d,d,x,d', comment='#', /silent
         if (doprint) then printf, lun, ' SN Ia in Hubble Flow read ' + $
                             'from file ' + sn1a_hf, format='(a)'
         if (sn1a_intrinsic ne 0) then err_sn1a_hf = $
            sqrt (err_sn1a_hf^2+sn1a_intrinsic^2)
         vhel_sn1a_hf = c * ((1+zhel_sn1a_hf)^2-1) / ((1+zhel_sn1a_hf)^2+1)
         vcmb_sn1a_hf = c * ((1+zcmb_sn1a_hf)^2-1) / ((1+zcmb_sn1a_hf)^2+1)
         case strupcase(vcorr_sn1a) of
            '2M++'  : begin
                  vcorr_sn1a_hf = c * ((1+zcorr_sn1a_hf)^2-1) / $
                                  ((1+zcorr_sn1a_hf)^2+1)
               end
            'CMB' : vcorr_sn1a_hf = vcmb_sn1a_hf
            else: begin
               if (doprint) then begin
                  printf, lun, ' Unknown value for vcorr_sn1a: ', vcorr_sn1a
                  printf, lun, ' Using default correction 2M++'
               endif
               vcorr_sn1a_hf = c * ((1+zcorr_sn1a_hf)^2-1) / $
                              ((1+zcorr_sn1a_hf)^2+1)
            end
         endcase

         alpha_sn1a_value = compute_alpha (app_sn1a_hf, err_sn1a_hf, $
            vhel=vhel_sn1a_hf, vcmb=vcmb_sn1a_hf, vcorr=vcorr_sn1a_hf, $
            vdisp = vpec_error, q0=q0, j0=j0, details=alphadet_sn1a, $
            redshift_range=redshift_range, alpha_err=alpha_sn1a_error, $
            chisq=chisq_sn1a_hf, ndof=ndof_sn1a_hf, fail=fail)
         printf, lun, 'Alpha from SN1a:', alpha_sn1a_value, alpha_sn1a_error, $
                 chisq_sn1a_hf, ndof_sn1a_hf, format='(a,3f12.6,i6)'
      endelse
      ; Here handle the case of fail ne 0
      if (fail) then begin
         if (doprint) then printf, lun, ' Failure in Hubble Flow processing' + $
                                   ' for SNe Ia; alpha error set to 100'
         alpha_sn1a_value = 0.d0
         alpha_sn1a_error = 100.d0
      endif
   endif
   ; alpha_sn1a_value and alpha_sn1a_error have been initialized
   ; Equation for H0 in terms of M0 and a_b:
   ; The basic equation is
   ;    alog10(h_0) == logh = 0.2*(sn_app_mag - mu_host) + a_b + 5
   ; Using the parameter abs_sn, the equations for each SN become:
   ;    abs_sn + mu_host = sn_app_mag  (included above)
   ; The equation tying the Hubble constant and abs_sn is then:
   ;    alog10(h_0) == logh = 0.2*abs_sn + a_b + 5
   ; recast as:
   ;    logh - 0.2*abs_sn = a_b+5
   ; Here a_b is the value of individual distances 
   ieq = ieq_sn1a_end ; The last equation is the Hubble Flow equation
   ieq_sn1a_h0 = ieq
   coeffs[ieq, iabs_sn1a] = -0.2d0
   coeffs[ieq, ihub] = 1.d0
   yval[ieq] = alpha_sn1a_value + 5
   yvar[ieq] = alpha_sn1a_error^2
   covar[ieq, ieq] = alpha_sn1a_error^2
   eq_descr[ieq] = ' Link between M_0 and H_0'
   eq_shape[ieq] = ' alog10(H0) - 0.2 * M_0 = a_b + 5'
   eq_covar[ieq] = ' Variance in Hubble flow magnitude (a_b)'
endif

; SN2 

if (do_sn2) then begin
   for k = 0, n_sn2-1 do begin
      ieq = k + ieq_sn2_start    ; this is the equation tying the SN
                    ; apparent magnitude and the host distance
                    ; to the sn absolute magnitude parameter
      ihost = ihost_sn2[k]    ; the index of the host -> mu_host
                    ; k is the SN index
      coeffs (ieq, ihost) = 1.d0
      coeffs (ieq, iabs_sn2) = 1.d0
      yval[ieq] = app_sn2[k]
      yvar[ieq] = sigma_sn2[k]^2
      covar[ieq,ieq] = sigma_sn2[k]^2
      ; app_sn2 = abs_1a + mu_host
      ; where app_sn2 is yval, abs_1a and mu_host are parameters
      ; and sigma_sn2 is the measurement error, which
      ; includes the intrinsic width of the SN distribution
      eq_descr[ieq] = ' Distance for SN ' + name_sn2[k] + ' in host ' + $
                      hosts[ihost_sn2[k]]
      eq_shape[ieq] = ' M0 + mu_host = sn_mag'
      eq_covar[ieq] = ' Variance of SN magnitude (incl intrinsic dispersion)'
   endfor

   ; Add the equation relating M0_1a to the Hubble constant
   ; Use alpha_sn2 (given value or computed from sn2_hf)

   if (read_sn2_hf) then begin
      readcol, sn2_hf, name_sn2_hf, zcorr_sn2_hf, app_sn2_hf, err_sn2_hf, $
               format='a,f,f,f', comment='#', /silent
      if (doprint) then printf, lun, ' SN II in Hubble Flow read ' + $
                        'from file ' + sn2_hf, format='(a)'
      vcorr_sn2_hf = zcorr_sn2_hf * (!const.c/1.e3)
      alpha_sn2_value = compute_alpha (app_sn2_hf, err_sn2_hf, $
            vcorr=vcorr_sn2_hf, redshift_range=redshift_range, $
            chisq=chisq_sn2_hf, ndof=ndof_sn2_hf, vdisp=vpec_error, $ 
            alpha_err=alpha_sn2_error, details=alphadet_sn2)
      printf, lun, 'Alpha from SN II:', alpha_sn2_value, alpha_sn2_error, $
                 chisq_sn2_hf, ndof_sn2_hf, format='(a,3f12.6,i6)'
   endif
   ; alpha_sn2_value and alpha_sn2_error have been initialized

   ; Equation for H0 in terms of M0 and a_b:
   ; The basic equation is
   ;    alog10(h_0) == logh = 0.2*(sn_app_mag - mu_host) + a_b + 5
   ; Using the parameter abs_sn, the equations for each SN become:
   ;    abs_sn + mu_host = sn_app_mag  (included above)
   ; The equation tying the Hubble constant and abs_sn is then:
   ;    alog10(h_0) == logh = 0.2*abs_sn + a_b + 5
   ; recast as:
   ;    logh - 0.2*abs_sn = a_b+5
   ; Here a_b is the value of individual distances 
   ieq = ieq_sn2_end   ; The last equation is the Hubble Flow equation
   ieq_sn2_h0 = ieq
   coeffs[ieq, iabs_sn2] = -0.2d0
   coeffs[ieq, ihub] = 1.d0
   yval[ieq] = alpha_sn2_value + 5
   yvar[ieq] = alpha_sn2_error^2
   covar[ieq, ieq] = alpha_sn2_error^2
   eq_descr[ieq] = ' Link between M_0 (SNe II) and H_0'
   eq_shape[ieq] = ' alog10(H0) - 0.2 * M_0 = a_b + 5'
   eq_covar[ieq] = ' Variance in Hubble flow magnitude (a_b) for SNe II'
endif

; Surface Brightness Fluctuations (SBF)

if (do_sbf) then begin
   ; For SBF, the host galaxy is also the calibrator.
   ; However, the information is used in
   ; different ways, so a copy of the
   ; vector is created for clarity.
   ; 
   ; The SBF calibration is based on the quantity big_m110 (Mbar_110),
   ; which is the expected absolute magnitude of SBF based on (g-z) color.
   ; The value of big_m110 in this file is based on the original 
   ; calibration in Jensen+ 2025, which can be shifted by an offset
   ; to be determined.  The parameter we fit for is abs_sbf_offset,   
   ; which is defined as the amount to be added to big M110 for the 
   ; final calibration.
   ;
   ; So the SBF calibration looks like:
   ;    m = M + M_offset + mu
   ; and the luminosity distance (modulus) estimate for objects in 
   ; Hubble flow is given by
   ;    mu = m - M - M_offset
   ; which then goes into the luminosity distance-redshift relation
   ;
   ; Each calibrator offers an equation relating the parameters 
   ; iabs_sbf and mu_host with a measured quantity (m-M):
   ;   M_offset + mu_host = m - M  (= yval)
   ;
   for k = 0, n_sbf-1 do begin
      ieq = k + ieq_sbf_start    ; this is the equation tying 
                    ; each SBF calibrator to its host distance
                    ; and to the SBF calibration offset.
                    ; Each equation has the form
                    ;   M_offset + mu_host = m - M  (= yval)
      ihost = ihost_sbf[k]       ; the index of the host -> mu_host
                    ; k is the SBF calibrator index
      coeffs (ieq, ihost) = 1.d0
      coeffs (ieq, iabs_sbf) = 1.d0
      yval[ieq] = small_m110_sbf[k] - big_m110_sbf[k]
      yvar[ieq] = sigma_big_m110_sbf[k]^2 + sigma_small_m110_sbf[k]^2 ; color + intrinsic and phot error
      covar[ieq,ieq] = sigma_big_m110_sbf[k]^2 + sigma_small_m110_sbf[k]^2 ; same
      eq_descr[ieq] = ' SBF calibrator distance for ' + $
                      name_sbf[k] + ' based on host distance from all sources'
      eq_shape[ieq] = ' M110_offset + mu_host = m110 - M110'
      eq_covar[ieq] = ' Variance of m110 and M110 (must include propagated color, intrinsic, ' + $
                      ' and photometry errors)'
   endfor

   ; SBF calibrator equations set.
   ; Now read and set up Hubble flow information, compute alpha, 
   ; and set up equation tying SBF calibration to Hubble flow.
   ; Format notes: the original file has name, color with error, 
   ; mbar110 (observed) with error, distance modulus with error, 
   ; distance with error, vgal, vgrp, vflow, v_2m++, v_rm, N_grp.
   ; The processed file will have: name; bigm110 with error;
   ; m110 with error; vcmb; vcorr.  bigm110 will be computed   
   ; from the color and error-propagated, with an extra term
   ; reflecting theintrinsic width of the distribution (to be
   ; determined).  m110 is directly as given in the original.
   ; vcmb is v_grp, assuming optical convention.  vcorr is 
   ; v_2m++, assuming optical convention.  We defer the issue of
   ; other velocity information.


   readcol, sbf_hf, name_sbf_hf, big_m110_sbf_hf, sigma_big_m110_sbf_hf, $
            m110_sbf_hf, sigma_m110_sbf_hf, $
            sbf_v_cmb, sbf_v_corr, $
            format = 'a,f,f,f,f,f,f', comment='#', /silent
   if (strlowcase(vcorr_sbf) eq 'cmb') then sbf_v_corr=sbf_v_cmb
   if (strlowcase(vcorr_sbf) ne 'cmb' and $
       strlowcase(vcorr_sbf) ne '2m++') then begin
      print, 'Unrecognized code for vcorr_sbf; using 2m++'
      vcorr_sbf = '2M++'
   endif
   
   if (doprint) then printf, lun, ' SBF systems in Hubble Flow read ' + $
                             'from file ' + sbf_hf, format='(a)'
   ; Interprets vgrp as vcmb, and v2m++ as vcorr.  Ignores the 
   ; heliocentric velocity correction
   ; 
   ; Assumes that the bigm110 error includes the measurement uncertainties
   ; for both velocity and parameters, as
   ; well as the intrinsic spread of the relation.
   ; It does NOT include the distance calibration uncertainty.
   ; m110 error is the photometric uncertainty.
   ; Errors are assumed uncorrelated between Hubble Flow objects.
   ; This is a simplification, as some uncertainties related to the
   ; internal SBF calibration parameters (NOT the absolute
   ; magnitude offset) would correlate uncertainties.
   ; Such terms are expected to be subdominant compared to the
   ; absolute magnitude calibration uncertainty.
   ; 
   ; The absolute magnitude calibration uncertainty is taken care of 
   ; directly as the uncertainty in the m_offset parameter.
   ; The above uses 2m++ as the default correction.  Additional
   ; calling parameters can be used to select a different
   ; peculiar velocity correction, redshift range, and so on.

   ; Note that for SBF, the quantity to
   ; form ALPHA is m110 minus bigm110.
   value = m110_sbf_hf - big_m110_sbf_hf
   alpha_sbf_value = compute_alpha (value, sigma_m110_sbf_hf, $
            vcorr = sbf_v_corr, $ ; vcmb = sbf_v_cmb, $
            redshift_range=redshift_range, alpha_err=alpha_sbf_error, $
            chisq=chisq_sbf_hf, ndof=ndof_sbf_hf, vdisp=sbf_vpec_error, $
            /optical, /ignore_cosmology, details=alphadet_sbf)
   printf, lun, 'Alpha from SBF:', alpha_sbf_value, alpha_sbf_error, $
           chisq_sbf_hf, ndof_sbf_hf, format='(a,3f12.6,i6)'
   
   ; Equation for H0 in terms of M0 and a_b:
   ; The basic equation is
   ;    alog10(h_0) == logh = 0.2*(sn_app_mag - mu_host) + a_b + 5
   ; Using the parameter abs_sn, the equations for each SN become:
   ;    abs_sn + mu_host = sn_app_mag  (included above)
   ; The equation tying the Hubble constant and abs_sn is then:
   ;    alog10(h_0) == logh = 0.2*abs_sn + a_b + 5
   ; recast as:
   ;    logh - 0.2*abs_sn = a_b+5
   ; Here a_b is the value of individual distances 
   ieq = ieq_sbf_end ; the last equation is the Hubble flow equation
   ieq_sbf_h0 = ieq
   coeffs[ieq, iabs_sbf] = -0.2d0
   coeffs[ieq, ihub] = 1.d0
   yval[ieq] = alpha_sbf_value + 5
   yvar[ieq] = alpha_sbf_error^2
   covar[ieq, ieq] = alpha_sbf_error^2
   eq_descr[ieq] = ' Link between M_offset(SBF) and H_0'
   eq_shape[ieq] = ' alog10(H0) - 0.2 * M_0 = a_b + 5'
   eq_covar[ieq] = ' Variance in Hubble flow magnitude (a_b)'
endif

; Tully-Fisher (TF)

if (do_tf) then begin
   ; For Tully_Fisher, the host galaxy is also the calibrator.
   ; However, the information is used in
   ; different ways, so a copy of the
   ; vector is created for clarity.
   ; 
   ; The TF calibration is based on the quantity big_m (M),
   ; which is the expected absolute magnitude based on the velocity
   ; width.  The value of big_m in this file is based on a 
   ; notional calibration, which can be shifted by an offset
   ; to be determined.  The parameter we fit for is abs_tf_offset,   
   ; which is defined as the amount to be added to big M for the 
   ; final calibration.
   ;
   ; So the TF calibration looks like:
   ;    m = M + M_offset + mu
   ; and the luminosity distance (modulus) estimate for objects in 
   ; Hubble flow is given by
   ;    mu = m - M - M_offset
   ; which then goes into the luminosity distance-redshift relation
   ;
   ; Each calibrator offers an equation relating the parameters 
   ; iabs_tf and mu_host with a measured quantity (m-M):
   ;   M_offset + mu_host = m - M  (= yval)
   ;
   for k = 0, n_tf-1 do begin
      ieq = k + ieq_tf_start    ; this is the equation tying 
                    ; each TF calibrator to its host distance
                    ; and to the Tf calibration offset.
                    ; Each equation has the form
                    ;   M_offset + mu_host = m - M  (= yval)
      ihost = ihost_tf[k]       ; the index of the host -> mu_host
                    ; k is the TF calibrator index
      coeffs (ieq, ihost) = 1.d0
      coeffs (ieq, iabs_tf) = 1.d0
      yval[ieq] = small_m_tf[k] - big_m_tf[k]
      yvar[ieq] = sigma_big_m_tf[k]^2  ; assume that this includes phot error
      covar[ieq,ieq] = sigma_big_m_tf[k]^2
      eq_descr[ieq] = ' Tully-Fisher calibrator distance for ' + $
                      name_tf[k] + ' based on host distance from all sources'
      eq_shape[ieq] = ' M_offset + mu_host = m - M'
      eq_covar[ieq] = ' Variance of m and M (must include velocity ' + $
                      ' and photometry errors)'
   endfor

   ; Add the equation relating M0_1a to the Hubble constant
   ; Use a_b_tf (given value or computed from tf_hf)

   readcol, tf_hf, name_tf_hf, big_m_tf_hf, small_m_tf_hf, $
            sigma_m_tf_hf, tf_v_cmb, tf_v_corr, $
            format = 'a,x,x,x,f,f,f,f,f', /silent
   if (doprint) then printf, lun, ' Tully-Fisher in Hubble Flow read ' + $
                             'from file ' + tf_hf, format='(a)'
   ; Assumes that the error includes the measurement uncertainties
   ; for both velocity and magnitude, but NOT the calibration uncertainty.
   ; Errors are assumed uncorrelated between Hubble Flow objects.
   ; This is a simplification, as some uncertainties related to the
   ; internal TF calibration parameters (NOT the absolute
   ; magnitude offset) would correlate uncertainties.
   ; Such terms are expected to be subdominant compared to the
   ; absolute magnitude calibration uncertainty.
   ; 
   ; The absolute magnitude calibration uncertainty is taken care of 
   ; directly as the uncertainty in the m_offset parameter.
   ; The above uses 2m++ as the default correction.  Additional
   ; calling parameters can be used to select a different
   ; peculiar velocity correction, redshift range, and so on.

   ; Note that for TF, the quantity to
   ; form ALPHA is small m minus big m.
   value = small_m_tf_hf - big_m_tf_hf
   alpha_tf_value = compute_alpha (value, sigma_m_tf_hf, vcorr=tf_v_corr, $
            redshift_range=redshift_range, alpha_err=alpha_tf_error, $
            chisq=chisq_tf_hf, ndof=ndof_tf_hf, vdisp=vpec_error, $
            details=alphadet_tf)
   printf, lun, 'Alpha from TF:', alpha_tf_value, alpha_tf_error, $
           chisq_tf_hf, ndof_tf_hf, format='(a,3f12.6,i6)'
   ;
   ; Equation for H0 in terms of M0 and a_b:
   ; The basic equation is
   ;    alog10(h_0) == logh = 0.2*(sn_app_mag - mu_host) + a_b + 5
   ; Using the parameter abs_sn, the equations for each SN become:
   ;    abs_sn + mu_host = sn_app_mag  (included above)
   ; The equation tying the Hubble constant and abs_sn is then:
   ;    alog10(h_0) == logh = 0.2*abs_sn + a_b + 5
   ; recast as:
   ;    logh - 0.2*abs_sn = a_b+5
   ; Here a_b is the value of individual distances 
   ieq = ieq_tf_end ; The last equation is the Hubble flow equation
   ieq_tf_h0 = ieq
   coeffs[ieq, iabs_tf] = -0.2d0
   coeffs[ieq, ihub] = 1.d0
   yval[ieq] = alpha_tf_value + 5
   yvar[ieq] = alpha_tf_error^2
   covar[ieq, ieq] = alpha_tf_error^2
   eq_descr[ieq] = ' Link between M_offset(TF) and H_0'
   eq_shape[ieq] = ' alog10(H0) - 0.2 * M_0 = a_b + 5'
   eq_covar[ieq] = ' Variance in Hubble flow magnitude (a_b)'
endif

; Add the equations for direct constraints.  These are in the form
; of redshifts, distaces, and errors.  The redshifts are converted
; into the quantity appearing in eq (4) of Riess+ 2022, which is
; log(redshift) times
; some kinematic terms that account for acceleration and jerk.
; Also, we need distance errors, which we assume uncorelated, and
; redshift errors, which may be correlated via the correction for
; peculiar motions.  The last issue is ignored in this version, but
; may become significant at some level.

if (do_mm) then begin
   ; Notes on MM computations
   ; The MM redshifts are given as velocities using the optical convention
   ; Therefore z = v*c (exactly, by definition)
   ; zhd = vcorr*c (that's how vcorr is defined)
   ; The distances they obtain are *angular diameter* distances, so
   ; df["v_corr"] * (1.0 - 0.5*(q0+3) * df["z"] - 1./6.*(10.0 + 7.0 * q0 + 3.0 * q0**2 - j0) * df["z"]**2)
   ; (from Nils, slack message on July 25)

   mm_v_corr_error = vpec_error
   mm_mu = 5*alog10(mm_distance) + 25  ; distance in Mpc
   mm_mu_error = 5*(mm_distance_error/mm_distance)/alog(10.)
   speed_of_light = !const.c / 1.e3
   mm_z = mm_v_corr / speed_of_light
   mm_cz_term = alog10(mm_v_corr) + $
                alog10 (1.d0 - (3.d0+q0)*mm_z/2.d0 + $
                        (11.d0+7.d0*q0+3.d0*q0^2-j0)*mm_z^2/6.d0)
   mm_v_total_error = sqrt (mm_v_error^2+mm_v_corr_error^2) 
   mm_cz_error = mm_v_total_error / mm_v_corr / alog(10.)
   mm_logh_value = mm_cz_term - 0.2*mm_mu  + 5
   mm_logh_error = sqrt ((0.2*mm_mu_error)^2+mm_cz_error^2)
   ;
   for k = 0, n_mm-1 do begin
      ieq = ieq_mm_start + k
      coeffs[ieq, ihub] = 1.d0
      yval [ieq] = mm_logh_value[k]
      covar[ieq,ieq] = mm_logh_error[k]^2
      eq_descr[ieq] = ' Megamaser constraint from ' + mm_name[k]
      eq_shape[ieq] = ' log10(H0) = measured value'
      eq_covar[ieq] = ' Reported variance of measurement'
   endfor
endif

if (do_epm) then begin
   epm_v_corr_error = vpec_error
   epm_v_corr_error = vpec_error
   epm_mu = 5*alog10(epm_distance) + 25  ; distance in Mpc
   epm_mu_error = 5*(epm_distance_error/epm_distance)/alog(10.)
   speed_of_light = !const.c / 1.e3
   epm_z = epm_v_corr / speed_of_light
   epm_cz_term = alog10(epm_v_corr * (1.d0 + (1.d0-q0)*epm_z/2 - $
                        (1.d0-q0-3*q0^2+j0)/6.d0*epm_z^2))
   epm_v_total_error = sqrt (epm_v_error^2+epm_v_corr_error^2) 
   epm_cz_error = epm_v_total_error / epm_v_corr / alog(10.)
   epm_logh_value = epm_cz_term - 0.2*epm_mu  + 5
   epm_logh_error = sqrt ((0.2*epm_mu_error)^2+epm_cz_error^2)
   ;
   for k = 0, n_epm-1 do begin
      ieq = ieq_epm_start + k
      coeffs[ieq, ihub] = 1.d0
      yval [ieq] = epm_logh_value[k]
      covar[ieq,ieq] = epm_logh_error[k]^2
      eq_descr[ieq] = ' EPM constraint from SN II named' + epm_name[k]
      eq_shape[ieq] = ' log10(H0) = measured value'
      eq_covar[ieq] = ' Reported variance of measurement'
   endfor
endif

; Future upgrade: add any covariance between distance or velocity
; measurements in MM or EPM groups.
; For megamasers, covariance should only come from velocity corrections.
; For EPM, it could come from systematic errors in the calibration.
; For now these are deemed subdominant to statistical uncertainties.  

; Now add the equations for Coma.  These include:
; One equation for each SN hosted in Coma, listed in the file 'coma'
;    if given
; The external (SBF) constraint
; The equation relating mu_coma to H0
;    which comes from the FP:
;    H0=(76.05 +/- 1.3)* 99.1/D_Coma  where D_Coma is in Mpc

; if (keyword_set(prior_coma) eq 0) then prior_coma = 0.13d0
; Nominal SBF prior on Coma distance; default value set in input.
; Can be set to a large value to remove Coma distance dependence 
; on TRGB-calibrated SBF.

if (n_coma gt 0) then begin
   for k = 0, n_coma-1 do begin
      ieq = ieq_coma_start+k  ; this is the equation constraining mu_coma from 
            ; the apparent magnitude of the k-th Coma supernova
            ; and abs_sn
            ; mu_coma = coma_sn_mag[k] - abs_sn
            ; or mu_coma + abs_sn = coma_sn_mag[k] with 
            ;   sigma = coma_sn_err[k]
      coeffs[ieq, icoma] = 1.d0
      coeffs[ieq, iabs_sn1a] = 1.d0
      yval[ieq] = coma_sn_mag[k]
      covar[ieq,ieq] = coma_sn_err[k]^2
      eq_descr[ieq] = ' Distance estimate to Coma SN '+coma_sn_name[k]
      eq_shape[ieq] = ' mu_coma - M_0 = observed SN magnitude'
      eq_covar[ieq] = ' Variance in apparent magnitude (incl intrinsic spread)'
   endfor
   ; Now the constraint from SBF : mu_coma = value + error
   ; This states that the value of 99.1 has an uncertainty of 0.13 mag
   ieq = ieq+1
   coeffs[ieq, icoma] = 1.d0
   mu_coma_nom = 25 + 5*alog10(d_coma_nom)
   yval[ieq] = mu_coma_nom
   covar[ieq,ieq] = prior_coma^2
   eq_descr[ieq] = ' SBF distance to Coma '
   eq_shape[ieq] = ' mu_coma = observed SBF value'
   eq_covar[ieq] = ' Variance in SBF distance modulus to Coma'
   ; And finally the equation for H0
   ; h0 = h0_coma_nom + d_coma_nom / d_coma
   ; h0_coma_nom = 76.05 is the default value set in input parameters
   logh_nom = alog10(h0_coma_nom)
   ; d_coma_nom = 99.1d0 default value set in input parameters
   mu_coma_nom = 5*alog10(d_coma_nom)+25
   ; These equations simply state that the FP
   ; value of h0 is 76.05 if D_coma is 99.1 Mpc
   ; 
   ; The corresponding equations are:
   ; alog10(h0) == logh = alog10(h0_coma_nom) + alog10(d_coma_nom/d_coma)
   ; logh = logh_nom + 0.2*(mu_coma_nom-mu_coma)
   ; logh + 0.2*mu_coma = logh_nom + 0.2_mu_coma_nom
   ieq = ieq+1
   coeffs[ieq,ihub] = 1.d0
   coeffs[ieq,icoma] = 0.2d0
   yval[ieq] = logh_nom + 0.2*mu_coma_nom
   covar[ieq,ieq] = (sigma_mag_fp/5.)^2 ; covariance in log10(h0); 
                 ; magnitude error divided by 5
   eq_descr[ieq] = ' H0 from Coma distance and Fundamental Plane'
   eq_shape[ieq] = ' alog10(H0) + 0.2*mu_coma = expected FP value'
   eq_covar[ieq] = ' Variance in cosmological FP value (assumed)'
endif

; This completes the matrices.  Now simply solve
; the linear equation.

; inv_covar = la_invert(covar) ; also called "precision" matrix
; (weight matrix when diagonal)

; FOR TEST PURPOSES ONLY: if FORCE_DIAG is set,
; then simply zero all off-diagonal elements of COVAR
; NOT TO BE USED FOR ANY PRODUCTION RUNS

if (keyword_set(force_diag)) then begin
   if (force_diag) then begin
      for k = 0, neq-2 do covar[k, k+1:neq-1]=0.d0
      for k = 0, neq-2 do covar[k+1:neq-1, k]=0.d0
   endif
endif
; stop

inv_covar = mpinv(covar, rank=rank) ; also called "precision" matrix (weight matrix when diagonal)

if (rank ne neq) then begin
   if (doprint) then printf, lun, 'Rank reduction: ', neq-rank
   cond1 = 0
   endif else cond1 = cond(covar)

solmat = transpose(coeffs) # inv_covar # coeffs 
invsolmat = la_invert (solmat)   ; also the covariance matrix of the solution
cond2 = cond(solmat)
solval = transpose(coeffs) # inv_covar # yval

params = invsolmat # solval

logh0_value = params[ihub]
logh0_var = invsolmat[ihub,ihub]

residuals = yval - coeffs#params
h0_value = 10^(logh0_value)
h0_error = 10^(logh0_value+sqrt(logh0_var))-h0_value
; ndof = neq - npars
ndof = rank - npars
chi2 = transpose(residuals) # inv_covar # residuals

; Output values of specific parameters
if (do_sn1a) then begin
   abs_sn1a_value = params[iabs_sn1a]
   abs_sn1a_error = sqrt(invsolmat[iabs_sn1a,iabs_sn1a])
   if (doprint) then printf, lun, ' Absolute magnitude offset for SNe Ia: ', $
                             abs_sn1a_value, ' +/- ', abs_sn1a_error, $
                             format='(2(a,f10.5))'
endif else begin
   abs_sn1a_value = 0.d0
   abs_sn1a_error = -1.d0
endelse   

if (do_sbf) then begin
   abs_sbf_value = params[iabs_sbf]
   abs_sbf_error = sqrt(invsolmat[iabs_sbf,iabs_sbf])
   if (doprint) then printf, lun, ' Absolute magnitude offset for SBF:  ', $
                             abs_sbf_value, ' +/- ', abs_sbf_error, $
                             format='(2(a,f10.5))'
endif else begin
   abs_sbf_value = 0.d0
   abs_sbf_error = -1.d0
endelse

if (do_sn2) then begin
   abs_sn2_value = params[iabs_sn2]
   abs_sn2_error = sqrt(invsolmat[iabs_sn2,iabs_sn2])
   if (doprint) then printf, lun, ' Absolute magnitude of SNe II: '+$
         '(empirically calibrated): ', abs_sn2_value, ' +/- ', abs_sn2_error, $
         format='(2(a,f14.5))'
endif else begin
   abs_sn2_value = 0.d0
   abs_sn2_error = -1.d0
endelse

if (do_tf) then begin
   abs_tf_value = params[iabs_tf]
   abs_tf_error = sqrt(invsolmat[iabs_tf,iabs_tf])
   if (doprint) then printf, lun, ' Absolute magnitude offset for TF: ', $
         abs_tf_value, ' +/- ', abs_tf_error, format='(2(a,f10.5))'
endif else begin
   abs_tf_value = 0.d0
   abs_tf_error = -1.d0
endelse

if (n_coma gt 0) then begin
   mu_coma_value = params[icoma]
   mu_coma_error = sqrt (invsolmat[icoma,icoma])
   if (doprint) then printf, lun, 'Distance modulus to Coma:   ', $
                             mu_coma_value, ' +/- ', mu_coma_error, $
                             format='(2(a,f10.5))'
endif
if (doprint) then printf, lun, ' H0 Value and uncertainty:  ', $
                          h0_value, ' +/- ', h0_error, h0_error/h0_value, $
                          ' %', 5*alog10(1+h0_error/h0_value), ' mag', $
                          format='(2(a,f10.5), 2(2x, f10.5, a))'

if (doprint) then printf, lun, ' Chi^2, # DOF            :  ', $
                          chi2, ndof, chi2/ndof, $
                          format='(a,f14.5,2x,i6,2x,f12.6)'

anchor_details = {n_anchors:n_anchors, anchors:anchors, $
                  mu_anchor_error:mu_anchor_error, $
                  data_anchors:data_anchors, data_mu_value:data_mu_value, $
                  data_mu_sigma:data_mu_sigma}
; 'anchors' are the anchors used in the run, and mu_anchor_error the 
; assumed geometric uncertainties.  'data_anchors' etc are the full
; anchor list, regardless of whether they are used or not.
; 'mu_anchor_value' (from 'data_mu_value') is not used yet; it could be
; used to assess the sensitivity of the process to the assumed
; anchor distances, but that's not a priority.

mas_details = {n_mas:n_mas, mas:mas, mas_error:mas_error, $
               mas_ref_name:mas_ref_name, mas_ref_value:mas_ref_value, $
               mas_ref_error:mas_ref_error, mas_ref_anchor:mas_ref_anchor, $
               mas_ref_method:mas_ref_method, mas_ref_source:mas_ref_source}
; These are the names of the method+anchor+source array items, and
; their assumed uncertainties.  They are currently extracted from the
; relevant lines in host_data.dat.  The index vector in host_data_details
; index mas/mas_error, which are matched to mas_ref_name etc.

host_data_details = {host2:host2, n2:n2, methodname_2:methodname_2, $
                     anchorname_2:anchorname_2, sourcename_2:sourcename_2, $
                     mu_host2_value:mu_host2_value, $
                     mu_host2_error:mu_host2_error, $
                     mas_name:mas_name, mas_index:mas_index, $
                     anchor_dist_error:anchor_dist_error, $
                     mas_host_error:mas_host_error, hosts:hosts, $
                     nhosts:nhosts, host_index_2:host_index_2}

if (do_sn1a) then begin
   sn1a_details = {do_sn1a:do_sn1a, n_sn1a:n_sn1a, host_sn1a:host_sn1a, $
                   name_sn1a:name_sn1a, app_sn1a:app_sn1a, $
                   ieq_sn1a_start:ieq_sn1a_start, ieq_sn1a_end:ieq_sn1a_end, $
                   sigma_sn1a:sigma_sn1a, ihost_sn1a:ihost_sn1a, $
                   read_sn1a_hf:read_sn1a_hf, vcorr_sn1a:vcorr_sn1a, $
                   name_sn1a_hf:name_sn1a_hf, $
                   alpha_sn1a_value:alpha_sn1a_value, $
                   alpha_sn1a_error:alpha_sn1a_error, $
                   alphadet_sn1a:alphadet_sn1a}
endif else sn1a_details = {do_sn1a:0B}


if (do_sn2) then begin
   sn2_details = {do_sn2:do_sn2, n_sn2:n_sn2, host_sn2:host_sn2, $
                  ieq_sn2_start:ieq_sn2_start, ieq_sn2_end:ieq_sn2_end, $
                  name_sn2:name_sn2, app_sn2:app_sn2, $
                  sigma_sn2:sigma_sn2, ihost_sn2:ihost_sn2, $
                  name_sn2_hf:name_sn2_hf, $
                  alpha_sn2_value:alpha_sn2_value, $
                  alpha_sn2_error:alpha_sn2_error, $
                  alphadet_sn2:alphadet_sn2}
endif else sn2_details = {do_sn2:0B}

if (do_tf) then begin
   tf_details = {do_tf:do_tf, n_tf:n_tf, host_tf:host_tf, $
                 ieq_tf_start:ieq_tf_start, ieq_tf_end:ieq_tf_end, $
                 name_tf:name_tf, small_m_tf:small_m_tf, $
                 big_m_tf:big_m_tf, sigma_big_m_tf:sigma_big_m_tf, $
                 ihost_tf:ihost_tf, $
                 name_tf_hf:name_tf_hf, $
                 alpha_tf_value:alpha_tf_value, $
                 alpha_tf_error:alpha_tf_error, $
                 alphadet_tf:alphadet_tf}
endif else tf_details = {do_tf:0B}

if (n_coma eq 0) then begin
   icoma=-1
   ieq_coma_start = -1
   coma_sn_name = ['']
   coma_sn_mag = [0]
   coma_sn_err = [0]
   mu_coma_value = 0
   mu_coma_error = 0
endif

coma_details = {ieq_coma_start:ieq_coma_start, n_coma:n_coma, icoma:icoma, $
                coma_sn_name:coma_sn_name, $
                coma_sn_mag:coma_sn_mag, coma_sn_err:coma_sn_err, $
                mu_coma_value:mu_coma_value, mu_coma_error:mu_coma_error, $
                prior_coma:prior_coma, h0_coma_nom:h0_coma_nom, $
                d_coma_nom:d_coma_nom, sigma_mag_fp:sigma_mag_fp}

if (n_mm eq 0) then begin
   mm_details = {n_mm:0, ieq_mm_start:-1}
endif else begin
   mm_details = {n_mm:n_mm, ieq_mm_start:ieq_mm_start, $
                 mm_name:mm_name, $
                 mm_distance:mm_distance, $
                 mm_distance_error:mm_distance_error, $
                 vcorr_mm:vcorr_mm, $
                 mm_v_obs:mm_v_obs, $
                 mm_v_corr:mm_v_corr, $
                 mm_v_error:mm_v_error, $
                 mm_mu:mm_mu, $
                 mm_mu_error:mm_mu_error, $
                 mm_z:mm_z, $
                 mm_v_corr_error:mm_v_corr_error, $
                 mm_v_total_error:mm_v_total_error, $
                 mm_cz_term:mm_cz_term, $
                 mm_cz_error:mm_cz_error, $
                 mm_logh_value:mm_logh_value, $
                 mm_logh_error:mm_logh_error}
endelse

if (n_epm eq 0) then begin
   epm_details = {n_epm:0, ieq_epm_start:-1}
endif else begin
   epm_details = {n_epm:n_epm, ieq_epm_start:ieq_epm_start, $
                 epm_name:epm_name, $
                 epm_distance:epm_distance, $
                 epm_distance_error:epm_distance_error, $
                 vcorr_epm:vcorr_epm, $
                 ; epm_v_obs:epm_v_obs, $
                 epm_v_corr:epm_v_corr, $
                 epm_v_error:epm_v_error, $
                 epm_mu:epm_mu, $
                 epm_mu_error:epm_mu_error, $
                 epm_z:epm_z, $
                 epm_v_corr_error:epm_v_corr_error, $
                 epm_v_total_error:epm_v_total_error, $
                 epm_cz_term:epm_cz_term, $
                 epm_cz_error:epm_cz_error, $
                 epm_logh_value:epm_logh_value, $
                 epm_logh_error:epm_logh_error}
endelse

if (n_sbf eq 0) then sbf_details = {do_sbf:0B, n_sbf:0, ieq_sbf_start:-1} else begin
   sbf_details = {do_sbf:do_sbf, n_sbf:n_sbf, ieq_sbf_start:ieq_sbf_start, $
                  ieq_sbf_end:ieq_sbf_end, iabs_sbf:iabs_sbf, $
                  name_sbf:name_sbf, host_sbf:host_sbf, ihost_sbf:ihost_sbf, $
                  big_m110_sbf:big_m110_sbf, $
                  sigma_big_m110_sbf:sigma_big_m110_sbf, $
                  small_m110_sbf:small_m110_sbf, $
                  sigma_small_m110_sbf:sigma_small_m110_sbf, $
                  name_sbf_hf:name_sbf_hf, $
                  big_m110_sbf_hf:big_m110_sbf_hf, $
                  sigma_big_m110_sbf_hf:sigma_big_m110_sbf_hf, $
                  m110_sbf_hf:m110_sbf_hf, $
                  sigma_m110_sbf_hf:sigma_m110_sbf_hf, $
                  sbf_v_cmb:sbf_v_cmb, sbf_v_corr:sbf_v_corr, $
                  vcorr_sbf:vcorr_sbf, alpha_sbf_value:alpha_sbf_value, $
                  alpha_sbf_error:alpha_sbf_error, alphadet_sbf:alphadet_sbf}
endelse
; logh0est is the value for h0 estimated for that specific Hubble flow object.
; It can be computed as follows:
; log(h0) - 0.2*M0 ~ alphavec + 5
; where in this case M0 is M0+offset
; or big_m110_sbf_hf + params_value[iabs_sbf]
; so that
; logh0 = 0.2*big_m110_sbf_hf + 0.2*params_value[iabs_sbf] + alphavec
; + 5

if (do_groups eq 0) then group_details = {ngroups:0} else begin
   group_details = {groups:groups, member_name:hosts[igh], $
                    group_name:groups[igg], mu_index:igp}
endelse


; Additional arrays: igh is the index of each host in groups in the
; overall hosts array; igg is the corresponding group for that host in
; the groups array.  The array igp is the distance modulus parameter
; for each group.
   
if (n_elements(host_class) eq 0) then host_class = ''
details = {h0_value:h0_value, h0_error:h0_error, chi2:chi2, ndof:ndof, $
           abs_sn1a_value:abs_sn1a_value, abs_sn1a_error:abs_sn1a_error, $
           abs_sbf_value:abs_sbf_value, abs_sbf_error:abs_sbf_error, $
           abs_sn2_value:abs_sn2_value, abs_sn2_error:abs_sn2_error, $
           abs_tf_value:abs_tf_value, abs_tf_error:abs_tf_error, $
           mu_coma_value:mu_coma_value, mu_coma_error:mu_coma_error, $
           npars:npars, parinfo:parinfo, params_name:parinfo[*].name, $
           params_value:params, params_var:invsolmat, $
           host_class:host_class, host_used:host_used, do_cull:do_cull, $
           neq:neq, covar:covar, inv_covar:inv_covar, coeffs:coeffs, $
           eq_descr:eq_descr, eq_shape:eq_shape, eq_covar:eq_covar, $
           cond1:cond1, cond2:cond2, yval:yval, residuals:residuals, $
           exclude_list:exclude_list, sn_exclude:sn_exclude, $
           ihub:ihub, icoma:icoma, iabs_sn1a:iabs_sn1a, $
           iabs_sbf:iabs_sbf, iabs_sn2:iabs_sn2, iabs_tf:iabs_tf, $
           anchor_details:anchor_details, mas_details:mas_details, $
           host_data_details:host_data_details, sn1a_details:sn1a_details, $
           sn2_details:sn2_details, tf_details:tf_details, $
           coma_details:coma_details, mm_details:mm_details, $
           epm_details:epm_details, sbf_details:sbf_details, $
           group_details:group_details, q0:q0, j0:j0}

if (dosave) then save, file=savefile, details

if (dofits) then begin
   mkhdr, header, '', /extend
   sxaddpar, header, 'H0_VALUE', h0_value
   sxaddpar, header, 'HO_ERROR', h0_error
   sxaddpar, header, 'HOST_DAT', host_data
   sxaddpar, header, 'SN1A_CAL', sn1a_calib
   if (do_coma) then sxaddpar, header, 'COMA', coma
   if (do_mm) then sxaddpar, header, 'MM', mm
   if (do_sbf) then sxaddpar, header, 'SBF', sbf_calib
   writefits, fitsfile, '', header
   ; mwrfits, ' ', fitsfile, header, /create
   ; Parameter list
   mkhdr, head1, params, /image
   sxaddpar, head1, 'xtension', 'PARAMS'
   for k = 0,npars-1 do begin
      sxaddpar, head1, 'PARAM'+string(k, format='(i03)'), parinfo[k].name
   endfor
   writefits, fitsfile, params, head1, /append
   ; mwrfits, params, fitsfile, head1
   ; Coefficients
   mkhdr, head2, coeffs, /image
   sxaddpar, head2, 'xtension', 'COEFFS'
   for k = 0,neq-1 do begin
      sxaddpar, head2, 'EQN'+string(k, format='(i05)'), eq_descr[k]
   endfor
   writefits, fitsfile, coeffs, head2, /append
   ; mwrfits, coeffs, fitsfile, head2
   ; Values
   mkhdr, head3, yval, /image
   sxaddpar, head3, 'xtension', 'DATA_VAL'
   writefits, fitsfile, yval, head3, /append
   ; mwrfits, yval, fitsfile, head3
   ; Covar
   mkhdr, head4, covar, /image
   sxaddpar, head4, 'xtension', 'COVAR'
   writefits, fitsfile, covar, head4, /append
   ; Sigma2
   mkhdr, head5, invsolmat, /image
   sxaddpar, head5, 'xtension', 'SIGMA2'
   writefits, fitsfile, invsolmat, head5, /append
endif

if (keyword_set(plabel) eq 0) then plabel=''
if (doplot) then make_residual_plot_v3, plotfile, details, plabel=plabel, /rota

closeout:
if (doprint) then begin
   close, lun
   free_lun, lun
endif
return

end
