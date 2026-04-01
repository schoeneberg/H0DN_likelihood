Pro vars, descr=descr_list, details=det_list, savefile=savefile, $
                     basedir=basedir, resid=resid, outdir=outdir

hvaluecmb = 67.24d0
herrorcmb =  0.35d0

; Some variants upon special request
refvar_array = !null
id_array = !null

; Done by hand in part to show the choices more clearly.
; if RESID is set, print residuals for baseline, kitchen sink with
; prefix RESID

descr_list = list()
det_list = list()
latexlabel = !null
inlatex = !null
if (keyword_set(outdir) eq 0) then outdir = './' else begin
   ; make the subdirectory if it does not exist
   nn = file_search (outdir, /test_directory, count=count)
   if (count eq 0) then spawn, 'mkdir ' + outdir
endelse
outdir = outdir+'/'  ; sometimes forgotten; no harm if doubled

if (keyword_set(basedir) eq 0) then basedir = './' else basedir=basedir+'/'

; In this version we switch to computing a_b via the alpha method.
; Since the defaults have not changed for v3p1, the file names for
; SN1a in the Hubble flow need to be given explicitly.  We expect that
; the default behavior will change for v3p2.
;
; Update 20250627: remove b25 to avoid using old Cepheid calibrations.

; V00 - baseline
; Baseline includes Cepheids, DEB, TRGB, NGC4258, Gaia, SNe Ia, SBF, Masers.
; Baseline does not include JAGB, Miras, SNe II calibrated astrophysically, 
; SNe II calibrated empirically, DESI FP, Tully-Fisher, SMC as anchor.
; These are included in variants.
;

; BASELINE

id = 'V00'
llab = 'base'
latexlabel = [latexlabel, llab]
label = id + ': Baseline '
descr_list.add, label
print, ' Processing '+label

org_v3p9, config='none', anchor_data='anchors.dat', basedir=basedir, $
          details=dd, /sbf_calib, groups_file='groups.dat', /mm, $
          sn1a_calib='sn1a_calib.dat', $
          sn1a_hf='sn1a_hf_pp.dat', sn1a_covar='sn1a_covar_pp.dat', $
          exclude_list=['SMC','*jagb*','*mira*','b25'], $
          log=outdir+id+'.txt' , plot=outdir+id+'.ps', plabel=label
det_list.add, dd
inlatex = [inlatex, 1B]
id_array = [id_array, id]
refvar = 'none'
refvar_array = [refvar_array, refvar]

if (keyword_set(resid)) then print_residuals, dd, $
   prefix=outdir+'baseline_'+resid

; 01-07: Add JAGB, Miras, NSFP, SNII, EPM, TF, SMC

; V01: add JAGB
id = 'V01'
llab = 'jagb'
latexlabel = [latexlabel, llab]
label = id + ': Baseline + JAGB'
descr_list.add, label
print, ' Processing '+label
org_v3p9, config='none', anchor_data='anchors.dat', basedir=basedir, $
          details=dd, /sbf_calib, groups_file='groups.dat', /mm, $
          sn1a_calib='sn1a_calib.dat', sn1a_hf='sn1a_hf_pp.dat', $
          sn1a_covar='sn1a_covar_pp.dat', $
          exclude_list=['SMC','*mira*','b25'], $
          log=outdir+id+'.txt', plot=outdir+id+'.ps', $
          plabel=label
det_list.add, dd
inlatex = [inlatex, 1B]
id_array = [id_array, id]
refvar = 'V00'
refvar_array = [refvar_array, refvar]

; V02: add Miras
id = 'V02'
llab = 'mira'
latexlabel = [latexlabel, llab]
label = id + ': Baseline + Miras'
descr_list.add, label
print, ' Processing '+label
org_v3p9, config='none', anchor_data='anchors.dat', basedir=basedir, $
          details=dd, /sbf_calib, groups_file='groups.dat', /mm, $
          sn1a_calib='sn1a_calib.dat', sn1a_hf='sn1a_hf_pp.dat', $
          sn1a_covar='sn1a_covar_pp.dat', $
          exclude_list=['SMC','*jagb*','b25'], $
          log=outdir+id+'.txt', plot=outdir+id+'.ps', $
          plabel=label
inlatex = [inlatex, 1B]
det_list.add, dd
id_array = [id_array, id]
refvar = 'V00'
refvar_array = [refvar_array, refvar]

; V03: add DESI FP
id = 'V03'
llab = 'fp'
latexlabel = [latexlabel, llab]
label = id + ': Baseline + DESI FP calibrated to Coma'
descr_list.add, label
print, ' Processing '+label
org_v3p9, config='none', anchor_data='anchors.dat', basedir=basedir, $
          details=dd, /sbf_calib, groups_file='groups.dat', /mm, $
          sn1a_calib='sn1a_calib.dat', sn1a_hf='sn1a_hf_pp.dat', $
          sn1a_covar='sn1a_covar_pp.dat', /coma, $
          exclude_list=['SMC','*jagb*','*mira*','b25'], $
          log=outdir+id+'.txt', plot=outdir+id+'.ps', $
          plabel=label
det_list.add, dd
inlatex = [inlatex, 1B]
id_array = [id_array, id]
refvar = 'V00'
refvar_array = [refvar_array, refvar]

; V04: Include empirically calibrated SNe II
id = 'V04'
llab = 'snii'
latexlabel = [latexlabel, llab]
label = id + ': Baseline + empirically calibrated SNe II'
descr_list.add, label
print, ' Processing '+label
org_v3p9, config='none', anchor_data='anchors.dat', basedir=basedir, $
          details=dd, /sbf_calib, groups_file='groups.dat', /mm, $
          sn1a_calib='sn1a_calib.dat', sn1a_hf='sn1a_hf_pp.dat', $
          sn1a_covar='sn1a_covar_pp.dat', /sn2_calib, $
          exclude_list=['SMC','*jagb*', '*mira*','b25'], $
          log=outdir+id+'.txt', plot=outdir+id+'.ps', $
          plabel=label
det_list.add, dd
inlatex = [inlatex, 1B]
id_array = [id_array, id]
refvar = 'V00'
refvar_array = [refvar_array, refvar]

; V04 -> V05: add astrophysically calibrated SNe II
id = 'V05'
llab = 'epm'
latexlabel = [latexlabel, llab]
label = id + ': Baseline + SNe II w. Expanding Photosphere method'
descr_list.add, label
print, ' Processing '+label
org_v3p9, config='none', anchor_data='anchors.dat', basedir=basedir, $
          details=dd, /sbf_calib, groups_file='groups.dat', /mm, /epm, $
          sn1a_calib='sn1a_calib.dat', sn1a_hf='sn1a_hf_pp.dat', $
          sn1a_covar='sn1a_covar_pp.dat', $
          exclude_list=['SMC','*jagb*','*mira*','b25'], $
          log=outdir+id+'.txt', plot=outdir+id+'.ps', $
          plabel=label
det_list.add, dd
inlatex = [inlatex, 1B]
id_array = [id_array, id]
refvar = 'V00'
refvar_array = [refvar_array, refvar]

; V06: include Tully-Fisher from CF4. Results are preliminary, error
; estimates need update.

id = 'V06'
llab = 'tf'
latexlabel = [latexlabel, llab]
label = id + ': Baseline + Tully-Fisher'
descr_list.add, label
print, ' Processing '+label
org_v3p9, config='none', anchor_data='anchors.dat', basedir=basedir, $
          details=dd, /sbf_calib, groups_file='groups.dat', /mm, /tf_calib, $
          sn1a_calib='sn1a_calib.dat', sn1a_hf='sn1a_hf_pp.dat', $
          sn1a_covar='sn1a_covar_pp.dat', $
          exclude_list=['SMC','*jagb*','*mira*','b25'], $
          log=outdir+id+'.txt', plot=outdir+id+'.ps', $
          plabel=label
det_list.add, dd
inlatex = [inlatex, 1B]
id_array = [id_array, id]
refvar = 'V00'
refvar_array = [refvar_array, refvar]


; V07: Include SMC as anchor
id = 'V07'
llab = 'smc'
latexlabel = [latexlabel, llab]
label = id + ': Baseline + SMC'
descr_list.add, label
print, ' Processing '+label
org_v3p9, config='none', anchor_data='anchors.dat', basedir=basedir, $
          details=dd, /sbf_calib, groups_file='groups.dat', /mm, $
          sn1a_calib='sn1a_calib.dat', sn1a_hf='sn1a_hf_pp.dat', $
          sn1a_covar='sn1a_covar_pp.dat', $
          exclude_list=['*jagb*','*mira*','b25'], $
          log=outdir+id+'.txt', plot=outdir+id+'.ps', $
          plabel=label
det_list.add, dd
inlatex = [inlatex, 1B]
id_array = [id_array, id]
refvar = 'V00'
refvar_array = [refvar_array, refvar]

; V08: remove all Cepheids.  Keep SBF (2025 uses TRGB calibration)

id = 'V08'
llab = 'noceph'
latexlabel = [latexlabel, llab]
label = id + ': Baseline w/o Cepheids'
descr_list.add, label
print, ' Processing '+label
org_v3p9, config='none', anchor_data='anchors.dat', basedir=basedir, $
          details=dd, /mm, /sbf_calib, groups_file='groups.dat', $
          sn1a_calib='sn1a_calib.dat', sn1a_hf='sn1a_hf_pp.dat', $
          sn1a_covar='sn1a_covar_pp.dat', $
          exclude_list=['SMC','*ceph*', '*jagb*','*mira*','b25'], $
          log=outdir+id+'.txt', plot=outdir+id+'.ps', $
          plabel=label
det_list.add, dd
names_08 = dd.sn1a_details.name_sn1a
inlatex = [inlatex, 1B]
id_array = [id_array, id]
refvar = 'V08B'
refvar_array = [refvar_array, refvar]

; V08B: Custom baseline for V08.  Exclude SN calibrators only measured
; via Cepheids; remaining SNe are as listed in array names_08 created under V08 above.
id = 'V08B'
llab = 'custbase_noceph'
latexlabel = [latexlabel, llab]
label = id + ': Custom baseline for V08'
descr_list.add, label
print, ' Processing '+label
org_v3p9, config='none', anchor_data='anchors.dat', basedir=basedir, $
          details=dd, /mm, /sbf_calib, groups_file='groups.dat', $
          sn1a_calib='sn1a_calib.dat', $
          sn1a_hf='sn1a_hf_pp.dat', sn1a_covar='sn1a_covar_pp.dat', $
          exclude_list=['SMC', '*jagb*','*mira*','b25'], $
          sn_include=names_08, log=outdir+id+'.txt', $
          plot=outdir+id+'.ps', plabel=label
det_list.add, dd
inlatex = [inlatex, 1B]
id_array = [id_array, id]
refvar = 'none'
refvar_array = [refvar_array, refvar]

; V09: remove all TRGB. Remove SBF (2025 uses TRGB calibration)
id = 'V09'
llab = 'notrgb'
latexlabel = [latexlabel, llab]
label = id + ': Baseline w/o TRGB, SBF'
descr_list.add, label
print, ' Processing '+label
org_v3p9, config='none', anchor_data='anchors.dat', basedir=basedir, $
          details=dd, /mm, sn1a_calib='sn1a_calib.dat', $
          sn1a_hf='sn1a_hf_pp.dat', sn1a_covar='sn1a_covar_pp.dat', $
          exclude_list=['SMC','*trgb*', '*jagb*', '*mira*','b25'], $
          log=outdir+id+'.txt', plot=outdir+id+'.ps', $
          plabel=label
det_list.add, dd
names_09 = dd.sn1a_details.name_sn1a
inlatex = [inlatex, 1B]
id_array = [id_array, id]
refvar = 'V09B'
refvar_array = [refvar_array, refvar]

; V09B: custom baseline for V09 (no TRGB, SBF).  Excludes SNe with
; only TRGB distance, but does not exclude SBF.  This is achieved
; by including only SN calibrators with data in V09, as listed in
; array names_09 created there
id = 'V09B'
llab = 'custbase_notrgb'
latexlabel = [latexlabel, llab]
label = id + ': Custom baseline for V09'
descr_list.add, label
print, ' Processing '+label
org_v3p9, config='none', anchor_data='anchors.dat', basedir=basedir, $
          details=dd, /mm, /sbf_calib, $
          groups_file='groups.dat', sn1a_calib='sn1a_calib.dat', $
          sn1a_hf='sn1a_hf_pp.dat', sn1a_covar='sn1a_covar_pp.dat', $
          exclude_list=['SMC', '*jagb*', '*mira*','b25'], $
          sn_include = names_09, $
          log=outdir+id+'.txt', plot=outdir+id+'.ps', $
          plabel=label
det_list.add, dd
inlatex = [inlatex, 1B]
id_array = [id_array, id]
refvar = 'none'
refvar_array = [refvar_array, refvar]

; V10: remove Gaia parallaxes (thus MW)
id = 'V10'
llab = 'nomw'
latexlabel = [latexlabel, llab]
label = id + ': Baseline w/o Gaia parallaxes, MW'
descr_list.add, label
print, ' Processing '+label
org_v3p9, config='none', anchor_data='anchors.dat', basedir=basedir, $
          details=dd, /sbf_calib, groups_file='groups.dat', /mm, $
          sn1a_calib='sn1a_calib.dat', sn1a_hf='sn1a_hf_pp.dat', $
          sn1a_covar='sn1a_covar_pp.dat', $
          exclude_list=['SMC','*MW*', '*jagb*', '*mira*','b25'], $
          log=outdir+id+'.txt', plot=outdir+id+'.ps', $
          plabel=label
det_list.add, dd
inlatex = [inlatex, 1B]
id_array = [id_array, id]
refvar = 'V00'
refvar_array = [refvar_array, refvar]

; V11: remove DEB (thus LMC+SMC)
id = 'V11'
llab = 'nolmc'
latexlabel = [latexlabel, llab]
label = id + ': Baseline w/o DEB, LMC'
descr_list.add, label
print, ' Processing '+label
org_v3p9, config='none', anchor_data='anchors.dat', basedir=basedir, $
          details=dd, /sbf_calib, groups_file='groups.dat', /mm, $
          sn1a_calib='sn1a_calib.dat', sn1a_hf='sn1a_hf_pp.dat', $
          sn1a_covar='sn1a_covar_pp.dat', $
          exclude_list=['SMC','*LMC*', '*jagb*', '*mira*','b25'], $
          log=outdir+id+'.txt', plot=outdir+id+'.ps', $
          plabel=label
det_list.add, dd
inlatex = [inlatex, 1B]
id_array = [id_array, id]
refvar = 'V00'
refvar_array = [refvar_array, refvar]

; V12: remove NGC 4258 (thus all TRGB)
id = 'V12'
llab = 'noftfe'
latexlabel = [latexlabel, llab]
label = id + ': Baseline w/o NGC 4258'
descr_list.add, label
print, ' Processing '+label
org_v3p9, config='none', anchor_data='anchors.dat', basedir=basedir, $
          details=dd, /mm, sn1a_calib='sn1a_calib.dat', $
          sn1a_hf='sn1a_hf_pp.dat', sn1a_covar='sn1a_covar_pp.dat', $
          exclude_list=['SMC','N4258', '*jagb*', '*mira*','b25'], $
          log=outdir+id+'.txt', plot=outdir+id+'.ps', $
          plabel=label
det_list.add, dd
names_12 = dd.sn1a_details.name_sn1a
inlatex = [inlatex, 1B]
id_array = [id_array, id]
refvar = 'V12B'
refvar_array = [refvar_array, refvar]

;  V12B: custom baseline for V12 (no NGC 4258).  This custom baseline
;  also excludes SBF, unlike V09B.
id = 'V12B'
llab = 'custbase_noftfe'
latexlabel = [latexlabel, llab]
label = id + ': Custom baseline for V12'
descr_list.add, label
print, ' Processing '+label
org_v3p9, config='none', anchor_data='anchors.dat', basedir=basedir, $
          details=dd, /mm, sn1a_calib='sn1a_calib.dat', $
          sn1a_hf='sn1a_hf_pp.dat', sn1a_covar='sn1a_covar_pp.dat', $
          exclude_list=['SMC', '*jagb*', '*mira*','b25'], $
          sn_include = names_12, $
          log=outdir+id+'.txt', plot=outdir+id+'.ps', $
          plabel=label
det_list.add, dd
inlatex = [inlatex, 1B]
id_array = [id_array, id]
refvar = 'none'
refvar_array = [refvar_array, refvar]


; V13: remove SNe Ia
id = 'V13'
llab = 'nosnia'
latexlabel = [latexlabel, llab]
label = id + ': Baseline w/o SNe Ia'
descr_list.add, label
print, ' Processing '+label
org_v3p9, config='none', anchor_data='anchors.dat', basedir=basedir, $
          details=dd, /sbf_calib, groups_file='groups.dat', /mm, $
          exclude_list=['SMC','*jagb*', '*mira*','b25'], $
          sn1a_calib='none', log=outdir+id+'.txt', plot=outdir+id+'.ps', $
          plabel=label
det_list.add, dd
inlatex = [inlatex, 1B]
id_array = [id_array, id]
refvar = 'V00'
refvar_array = [refvar_array, refvar]

; V14: remove SBF
id = 'V14'
llab = 'nosbf'
latexlabel = [latexlabel, llab]
label = id + ': Baseline w/o SBF'
descr_list.add, label
print, ' Processing '+label
org_v3p9, config='none', anchor_data='anchors.dat', basedir=basedir, $
          details=dd, /mm, sn1a_calib='sn1a_calib.dat', $
          sn1a_hf='sn1a_hf_pp.dat', sn1a_covar='sn1a_covar_pp.dat', $
          exclude_list=['SMC','*jagb*', '*mira*','b25'], $
          log=outdir+id+'.txt', plot=outdir+id+'.ps', $
          plabel=label
det_list.add, dd
inlatex = [inlatex, 1B]
id_array = [id_array, id]
refvar = 'V00'
refvar_array = [refvar_array, refvar]

; V15: remove Megamasers
id = 'V15'
llab = 'nomm'
latexlabel = [latexlabel, llab]
label = id + ': Baseline w/o masers in the Hubble flow'
descr_list.add, label
print, ' Processing '+label
org_v3p9, config='none', anchor_data='anchors.dat', basedir=basedir, $
          details=dd, /sbf_calib, $
          groups_file='groups.dat', sn1a_calib='sn1a_calib.dat', $
          sn1a_hf='sn1a_hf_pp.dat', sn1a_covar='sn1a_covar_pp.dat', $
          exclude_list=['SMC','*jagb*', '*mira*','b25'], $
          log=outdir+id+'.txt', plot=outdir+id+'.ps', $
          plabel=label
det_list.add, dd
inlatex = [inlatex, 1B]
id_array = [id_array, id]
refvar = 'V00'
refvar_array = [refvar_array, refvar]

; V16: Exclude HST data
id = 'V16'
llab = 'nohst'
latexlabel = [latexlabel, llab]
label = id + ': Exclude HST data, LMC, MW'
descr_list.add, label
print, ' Processing '+label
org_v3p9, config='none', anchor_data='anchors.dat', basedir=basedir, $
          details=dd, /sbf_calib, groups_file='groups.dat', /mm, $
          sn1a_calib='sn1a_calib.dat', sn1a_hf='sn1a_hf_pp.dat', $
          sn1a_covar='sn1a_covar_pp.dat', $
          exclude_list=['SMC','*jagb*', '*mira*','b25','*hst*','MW','LMC'], $
          log=outdir+id+'.txt', plot=outdir+id+'.ps', plabel=label
det_list.add, dd
names_16 = dd.sn1a_details.name_sn1a
inlatex = [inlatex, 1B]
id_array = [id_array, id]
refvar = 'V16B'
refvar_array = [refvar_array, refvar]

; V16B: Custom baseline for V16.  Excludes SNe with only HST measurements
; (formally includes only SNe with JWST measurements, listed in array names_16
; defined in variant V16 above)
id = 'V16B'
llab = 'custbase_nohst'
latexlabel = [latexlabel, llab]
label = id + ': Custom baseline for V16'
descr_list.add, label
print, ' Processing '+label
org_v3p9, config='none', anchor_data='anchors.dat', basedir=basedir, $
          details=dd, /sbf_calib, groups_file='groups.dat', /mm, $
          sn1a_calib='sn1a_calib.dat', sn1a_hf='sn1a_hf_pp.dat', $
          sn1a_covar='sn1a_covar_pp.dat', $
          exclude_list=['SMC','*jagb*', '*mira*','b25'], $
          sn_include = names_16, $
          log=outdir+id+'.txt', plot=outdir+id+'.ps', plabel=label
det_list.add, dd
inlatex = [inlatex, 1B]
id_array = [id_array, id]
refvar = 'none'
refvar_array = [refvar_array, refvar]

; V17: Exclude JWST data
id = 'V17'
llab = 'nojwst'
latexlabel = [latexlabel, llab]
label = id + ': Exclude JWST data'
descr_list.add, label
print, ' Processing '+label
org_v3p9, config='none', anchor_data='anchors.dat', basedir=basedir, $
          details=dd, /mm, sn1a_calib='sn1a_calib.dat', $
          sn1a_hf='sn1a_hf_pp.dat', sn1a_covar='sn1a_covar_pp.dat', $
          exclude_list=['SMC','*jagb*', '*mira*','b25','*jwst*'], $
          log=outdir+id+'.txt', plot=outdir+id+'.ps', $
          plabel=label
; /sbf_calib, groups_file='groups.dat', 
det_list.add, dd
inlatex = [inlatex, 1B]
id_array = [id_array, id]
refvar = 'V00'
refvar_array = [refvar_array, refvar]

; V18: Remove old SNe

id = 'V18'
llab = 'nooldsne'
latexlabel = [latexlabel, llab]
label = id + ': Exclude SN 1994D and earlier'
descr_list.add, label
print, ' Processing '+label
org_v3p9, config='none', anchor_data='anchors.dat', basedir=basedir, $
          details=dd, /sbf_calib, groups_file='groups.dat', /mm, $
          sn1a_calib='sn1a_calib.dat', sn1a_hf='sn1a_hf_pp.dat', $
          sn1a_covar='sn1a_covar_pp.dat', $
          exclude_list=['SMC','*jagb*', '*mira*','b25'], $
          sn_exclude=['198*','199[0-3]*','1994D*'], $
          log=outdir+id+'.txt', plot=outdir+id+'.ps', $
          plabel=label
det_list.add, dd
inlatex = [inlatex, 1B]
id_array = [id_array, id]
refvar = 'V00'
refvar_array = [refvar_array, refvar]

; No SNe before 2000.  Not an official variant, not included in the
; LaTeX tables.
id = 'V18A'
llab = 'nooldsnevar'
latexlabel = [latexlabel, llab]
label = id + ': Exclude SNe before 2000'
descr_list.add, label
print, ' Processing '+label
org_v3p9, config='none', anchor_data='anchors.dat', basedir=basedir, $
          details=dd, /sbf_calib, groups_file='groups.dat', /mm, $
          sn1a_calib='sn1a_calib.dat', sn1a_hf='sn1a_hf_pp.dat', $
          sn1a_covar='sn1a_covar_pp.dat', $
          exclude_list=['SMC','*jagb*', '*mira*','b25'], $
          sn_exclude=['198*','199*'], $
          log=outdir+id+'.txt', plot=outdir+id+'.ps', $
          plabel=label
det_list.add, dd
inlatex = [inlatex, 0B]
id_array = [id_array, id]
refvar = 'V00'
refvar_array = [refvar_array, refvar]

; V15 -> V19: Remove peculiar velocity corrections
; Cannot be done fully yet; requires appropriate value of a_B, error
; Partial version uses SBF, MM w/o corrections
; Check that these are in the CMB frame

id = 'V19'
llab = 'usecmb'
latexlabel = [latexlabel, llab]
label = id + ': Baseline w. SNe Ia, SBF, Masers in CMB frame'
descr_list.add, label
print, ' Processing '+label
org_v3p9, config='none', anchor_data='anchors.dat', basedir=basedir, $
          details=dd, /sbf_calib, groups_file='groups.dat', /mm, $
          sn1a_calib='sn1a_calib.dat', sn1a_hf='sn1a_hf_pp.dat', $
          sn1a_covar='sn1a_covar_pp.dat', vcorr_sn1a='CMB', $
          vcorr_mm='group', vcorr_sbf='cmb', $
          exclude_list=['SMC','*jagb*', '*mira*','b25'], $
          log=outdir+id+'.txt', plot=outdir+id+'.ps', $
          plabel=label
det_list.add, dd
inlatex = [inlatex, 1B]
id_array = [id_array, id]
refvar = 'none'
refvar_array = [refvar_array, refvar]

; V20: remove Hubble flow objects w. z < 0.06
; Also exclude SBF
id = 'V20'
llab = 'nohighz'
latexlabel = [latexlabel, llab]
label = id + ': Hubble flow SNe Ia w. z < 0.06'
descr_list.add, label
print, ' Processing '+label
org_v3p9, config='none', anchor_data='anchors.dat', basedir=basedir, $
          details=dd, /mm, sn1a_calib='sn1a_calib.dat', $
          sn1a_hf='sn1a_hf_pp.dat', sn1a_covar='sn1a_covar_pp.dat', $
          exclude_list=['SMC','*jagb*', '*mira*','b25'], $          
          min_redshift=0.023d0, max_redshift=0.06d0, $
          log=outdir+id+'.txt', plot=outdir+id+'.ps', $
          plabel=label
det_list.add, dd
inlatex = [inlatex, 1B]
id_array = [id_array, id]
refvar = 'none'
refvar_array = [refvar_array, refvar]

; V19 -> none: TRGB w. OGLE calibration
; Cannot be done yet


; V28 -> V21 : Pantheon+ SNe w. direct a_B computation and redshift selection
; Exclude SBF
id = 'V21'
llab = 'nolowz'
latexlabel = [latexlabel, llab]
label = id + ': SNe Ia in redshift range 0.03--0.10 '
descr_list.add, label
print, ' Processing '+label
org_v3p9, config='none', anchor_data='anchors.dat', basedir=basedir, $
          details=dd, /mm, sn1a_calib='sn1a_calib.dat', $
          sn1a_hf='sn1a_hf_pp.dat', sn1a_covar='sn1a_covar_pp.dat', $
          exclude_list=['SMC','*jagb*', '*mira*','b25'], $
          min_redshift = 0.03, max_redshift = 0.1, $
          log=outdir+id+'.txt', plot=outdir+id+'.ps', $
          plabel=label
det_list.add, dd
inlatex = [inlatex, 1B]
id_array = [id_array, id]
refvar = 'none'
refvar_array = [refvar_array, refvar]

; V22 : CSP/Snoopy instead of Panheon+
id = 'V22'
llab = 'usecsp'
latexlabel = [latexlabel, llab]
label = id + ': SNe Ia from CSP, all 55 calibrators'
descr_list.add, label
print, ' Processing '+label
org_v3p9, config='none', anchor_data='anchors.dat', basedir=basedir, $
          details=dd, /sbf_calib, groups_file='groups.dat', /mm, $
          exclude_list=['SMC','*jagb*', '*mira*','b25'], $
          sn1a_calib='sn1a_cal_csp_burns.dat', $
          sn1a_hf='hf_csp_lluis_v2.dat', $
          sn1a_covar='hf_csp_lluis_v2_cov.dat', $
          log=outdir+id+'.txt', plot=outdir+id+'.ps', $
          plabel=label
det_list.add, dd
inlatex = [inlatex, 1B]
id_array = [id_array, id]
refvar = 'none'
refvar_array = [refvar_array, refvar]

; V23 : Bayes SN instead of Panheon+
id = 'V23'
llab = 'usebayes'
latexlabel = [latexlabel, llab]
label = id + ': SNe Ia from BayesSN' 
descr_list.add, label
print, ' Processing '+label
org_v3p9, config='none', anchor_data='anchors.dat', basedir=basedir, $
          details=dd, /sbf_calib, groups_file='groups.dat', /mm, $
          sn1a_calib='sn1a_calib_bayessn.dat', sn1a_hf='sn1a_hf_bayessn.dat', $
          sn1a_covar='sn1a_covar_bayessn.dat', sn1a_intrinsic=0.15d0, $
          exclude_list=['SMC','*jagb*', '*mira*','b25'], $
          log=outdir+id+'.txt', plot=outdir+id+'.ps', plabel=label
det_list.add, dd
inlatex = [inlatex, 1B]
id_array = [id_array, id]
refvar = 'none'
refvar_array = [refvar_array, refvar]

; V24 : SALT3
id = 'V24'
llab = 'usesaltthree'
latexlabel = [latexlabel, llab]
label = id + ': SNe Ia from PP processed w. SALT3 fitter' 
descr_list.add, label
print, ' Processing '+label
org_v3p9, config='none', anchor_data='anchors.dat', basedir=basedir, $
          details=dd, /sbf_calib, groups_file='groups.dat', /mm, $
          sn1a_calib='rung3_salt3.dat', sn1a_hf='hf_salt3.dat', $
          sn1a_covar='hf_salt3_cov.dat', $
          exclude_list=['SMC','*jagb*', '*mira*','b25'], $
          log=outdir+id+'.txt', plot=outdir+id+'.ps', $
          plabel=label
det_list.add, dd
inlatex = [inlatex, 1B]
id_array = [id_array, id]
refvar = 'none'
refvar_array = [refvar_array, refvar]

; V25 : NIR/SNe (Galbany) instead of Panheon+.  Default H band
id = 'V25'
llab = 'usehband'
latexlabel = [latexlabel, llab]
label = id + ': SNe Ia in H-band'
descr_list.add, label
print, ' Processing '+label
org_v3p9, config='none', anchor_data='anchors.dat', basedir=basedir, $
          details=dd, /sbf_calib, groups_file='groups.dat', /mm, $
          sn1a_calib='sn1a_calib_H.dat', sn1a_hf='sn1a_hf_H.dat', $
          sn1a_intrinsic=0.096d0, $
          exclude_list=['SMC','*jagb*', '*mira*','b25'], $
          log=outdir+id+'.txt', plot=outdir+id+'.ps', $
          plabel=label
det_list.add, dd
inlatex = [inlatex, 1B]
id_array = [id_array, id]
refvar = 'none'
refvar_array = [refvar_array, refvar]

; V26 : NIR/SNe (Galbany) instead of Panheon+.  Optional J band
id = 'V26'
llab = 'usejband'
latexlabel = [latexlabel, llab]
label = id + ': SNe Ia in J-band'
descr_list.add, label
print, ' Processing '+label
org_v3p9, config='none', anchor_data='anchors.dat', basedir=basedir, $
          details=dd, /sbf_calib, groups_file='groups.dat', /mm, $
          sn1a_calib='sn1a_calib_J.dat', sn1a_hf='sn1a_hf_J.dat', $
          sn1a_intrinsic=0.125d0, $
          exclude_list=['SMC','*jagb*', '*mira*','b25'], $
          log=outdir+id+'.txt', plot=outdir+id+'.ps', $
          plabel=label
det_list.add, dd
inlatex = [inlatex, 1B]
id_array = [id_array, id]
refvar = 'none'
refvar_array = [refvar_array, refvar]

; V27: Pantheon+ SNe ignoring off-diagonal covariance in Hubble flow
id = 'V27'
llab = 'nocovar'
latexlabel = [latexlabel, llab]
label = id + ': No off-diagonal covariance for SNe Ia in Hubble flow'
descr_list.add, label
print, ' Processing '+label
org_v3p9, config='none', anchor_data='anchors.dat', basedir=basedir, $
          details=dd, /sbf_calib, groups_file='groups.dat', /mm, $
          sn1a_calib='sn1a_calib.dat', sn1a_hf='sn1a_hf_pp.dat', $
          sn1a_covar='sn1a_covar_pp.dat', /sn1a_ignore_offdiag, $
          exclude_list=['SMC','*jagb*', '*mira*','b25'], $
          log=outdir+id+'.txt', plot=outdir+id+'.ps', $
          plabel=label
det_list.add, dd
inlatex = [inlatex, 1B]
id_array = [id_array, id]
refvar = 'none'
refvar_array = [refvar_array, refvar]

; V28 No metallicity correction
id = 'V28'
llab = 'nometcorr'
latexlabel = [latexlabel, llab]
label = id + ': No metallicity correction for Cepheid PL'
descr_list.add, label
print, ' Processing '+label
org_v3p9, config='none', anchor_data='anchors.dat', basedir=basedir, $
          details=dd, /sbf_calib, groups_file='groups.dat', /mm, $
          host_data='host_data_nomet.dat', sn1a_calib='sn1a_calib.dat', $
          sn1a_hf='sn1a_hf_pp.dat', sn1a_covar='sn1a_covar_pp.dat', $
          exclude_list=['SMC','*jagb*','*mira*','b25'], $
          log=outdir+id+'.txt' , plot=outdir+id+'.ps', plabel=label ; , $
          ; fitsfile = id+'.fits'
det_list.add, dd
inlatex = [inlatex, 1B]
id_array = [id_array, id]
refvar = 'none'
refvar_array = [refvar_array, refvar]

; V97: All SN1a, nothing else
id = 'V97'
llab = 'onlysnia'
latexlabel = [latexlabel, llab]
label = id + ': Only SNe Ia in the Hubble flow '
descr_list.add, label
print, ' Processing '+label
org_v3p9, config='none', anchor_data='anchors.dat', basedir=basedir, $
          details=dd, sn1a_calib='sn1a_calib.dat', $
          sn1a_hf='sn1a_hf_pp.dat', sn1a_covar='sn1a_covar_pp.dat', $
          exclude_list=['smc','b25'], log=outdir+id+'.txt', $
          plot=outdir+id+'.ps', plabel=label
det_list.add, dd
inlatex = [inlatex, 0B]
id_array = [id_array, id]
refvar = 'none'
refvar_array = [refvar_array, refvar]

; V98: Kitchen sink minus SNe Ia
; Note: as a stop gap, this is achieved by setting the parameter
; alpha_sn1a_error to a very large number.  For SN1A, simply setting
; the hubble flow file to zero reverts to the default values of
; apha_sn1a and its uncertainty, unless overridden.  This behavior
; needs to change, and alpha be treated the same way as other Hubble
; Flow objects.
id = 'V98a'
llab = 'hfnosnia'
latexlabel = [latexlabel, llab]
label = id + ': All except SNe Ia in Hubble flow'
descr_list.add, label
print, ' Processing '+label
org_v3p9, config='none', anchor_data='anchors.dat', basedir=basedir, $
          details=dd, sn1a_hf='none', alpha_sn1a_error=100., $
          /sbf_calib, groups_file='groups.dat', $
          /mm, /epm, /sn2_calib, /coma, /tf_calib, $
          exclude_list='b25', log=outdir+id+'.txt', $
          plot=outdir+id+'.ps', plabel=label
          
det_list.add, dd
inlatex = [inlatex, 0B]
id_array = [id_array, id]
refvar = 'none'
refvar_array = [refvar_array, refvar]

; V98b: Kitchen sink minus SNe Ia and TF
; Note: as a stop gap, this is achieved by setting the parameter
; alpha_sn1a_error to a very large number.  For SN1A, simply setting
; the Hubble flow file to zero reverts to the default values of
; apha_sn1a and its uncertainty, unless overridden.  Setting the sn1a
; calibrator file to none excludes sn1a altogehter, which means no
; SN constraint for Coma and has other impacts in the full distance
; network.  This behavior will  be changed in future versions, and alpha 
; will be treated the same way as other Hubble Flow objects.

id = 'V98b'
llab = 'hftwo'
latexlabel = [latexlabel, llab]
label = id + ': All except SNe Ia and TF in Hubble flow'
descr_list.add, label
print, ' Processing '+label
org_v3p9, config='none', anchor_data='anchors.dat', basedir=basedir, $
          details=dd, sn1a_hf='none',  alpha_sn1a_error=100., $
          /sbf_calib, groups_file='groups.dat', $
          /mm, /epm, /sn2_calib, /coma, exclude_list='b25', $
          log=outdir+id+'.txt', plot=outdir+id+'.ps',  plabel=label
det_list.add, dd
inlatex = [inlatex, 0B]
id_array = [id_array, id]
refvar = 'none'
refvar_array = [refvar_array, refvar]

; V98c: Only TF in the Hubble flow
; Note: as a stop gap, this is achieved by setting the parameter
; alpha_sn1a_error to a very large number.  For SN1A, simply setting
; the hubble flow file to zero reverts to the default values of
; apha_sn1a and its uncertainty, unless overridden.  This behavior
; will be changed in future versions, and alpha will treated 
; the same way as other Hubble Flow objects.
id = 'V98c'
llab = 'hfonlytf'
latexlabel = [latexlabel, llab]
label = id + ': Only TF in Hubble flow'
descr_list.add, label
print, ' Processing '+label
org_v3p9, config='none', anchor_data='anchors.dat', basedir=basedir, $
          details=dd, sn1a_hf='none', alpha_sn1a_error=100., $
          /tf_calib, exclude_list='b25', $
          log=outdir+id+'.txt', plot=outdir+id+'.ps', plabel=label
det_list.add, dd
inlatex = [inlatex, 0B]
id_array = [id_array, id]
refvar = 'none'
refvar_array = [refvar_array, refvar]


; V99: Kitchen sink: all available data, but only one set of SN Ia
id = 'V99'
llab = 'everything'
latexlabel = [latexlabel, llab]
label = id + ': Everything available'
descr_list.add, label
print, ' Processing '+label
org_v3p9, config='none', anchor_data='anchors.dat', basedir=basedir,  $
          details=dd, /sbf_calib, $
          groups_file='groups.dat', /mm, /epm, /sn2_calib, /coma, $
          sn1a_calib='sn1a_calib.dat', sn1a_hf='sn1a_hf_pp.dat', $
          sn1a_covar='sn1a_covar_pp.dat', /tf_calib, exclude_list='b25', $
          log=outdir+id+'.txt', plot=outdir+id+'.ps', plabel=label
det_list.add, dd
inlatex = [inlatex, 1B]
if (keyword_set(resid)) then print_residuals, dd, $
   prefix=outdir+'everything_'+resid
id_array = [id_array, id]
refvar = 'V00'
refvar_array = [refvar_array, refvar]


; V99A: Kitchen sink minus TF
id = 'V99A'
llab = 'everythingnotf'
latexlabel = [latexlabel, llab]
label = id + ': Everything available except TF'
descr_list.add, label
print, ' Processing '+label
org_v3p9, config='none', anchor_data='anchors.dat', basedir=basedir, $
            details=dd, /sbf_calib, $
          groups_file='groups.dat', /mm, /epm, /sn2_calib, /coma, $
          sn1a_calib='sn1a_calib.dat', sn1a_hf='sn1a_hf_pp.dat', $
          sn1a_covar='sn1a_covar_pp.dat', exclude_list='b25', $
          log=outdir+id+'.txt', plot=outdir+id+'.ps', plabel=label
det_list.add, dd
inlatex = [inlatex, 1B]
id_array = [id_array, id]
refvar = 'V00'
refvar_array = [refvar_array, refvar]

;=================================================
; Variants complete.  Create LaTeX and log outputs
;=================================================

; The following sections produce LaTeX files for use in the paper.
; Two are tables that are inserted directly into the text.
; The third is a set of variable definitions that allows the values
; of H0, irts error, and other parameters to be called as variables
; in the LaTeX source, and thus be updated automatically.  Note that
; numbers are not allowed in variable names with the default TeX catcodes;
; so the variable names are all letters.  Also note that in at least
; one case, variables do not have a unique abbreviation; (La)TeX does not
; allow abbreviations, so this is not an issue per se, but it could trip
; up searches and other tools.

nvar = n_elements(det_list)

if (keyword_set(savefile)) then save, file=outdir+savefile, det_list, $
                         descr_list, latexlabel, inlatex

openw, lun, outdir+'vars_out.txt', /get_lun
openw, lun2, outdir+'vartable.tex', /get_lun
openw, lun3, outdir+'vardef.tex', /get_lun
; Table header
; printf, lun2, '\newcommand{\ptt}[1]{\parbox[t]{4cm}{\raggedright #1}}'
printf, lun2, '\begin{table*}'
printf, lun2, '\caption{Main variants for $ H_0 $ calculation ' + $
              '\label{tab:variants}}'  
printf, lun2, '\footnotesize\centering'
printf, lun2, '\begin{tabular}{lccrrcrl}'
printf, lun2, '\hline\hline'
printf, lun2, '\noalign{\smallskip}'
printf, lun2, '\# & $ H_0 $ & 1--$\sigma$ & $ N_{\rm dof} $ '+$
              '& Reduced &  PTE$^a$ & $ \Delta_\sigma^b $ & Description \\'
printf, lun2, '           & \multicolumn{2}{c}{\Hunit} & & '+$
        '$ \chi^2 $ \\'
; printf, lun2, '\# & $ H_0 $ & 1--$\sigma$ &  $ \chi^2 $ & $ N_{\rm dof} $ '+$
;               '& Reduced & $ M_0 $ & $ \alpha_{SN Ia} $ & '+$
;               '$ N_{\rm calib} $ & Description \\'
; printf, lun2, '           & \multicolumn{2}{c}{$ \rm km/s/Mpc $} & & & '+$
;         '$ \chi^2 $ & & & (SNe Ia)   \\'
printf, lun2, '\noalign{\smallskip}'
printf, lun2, '\hline'
printf, lun2, '\noalign{\smallskip}'
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

   ; Significance of difference.  Only valid when comparing variants that 
   ; are a strict extension of one another, i.e., such that ove variant 
   ; completely encompasses the other.  Then the significance of the 
   ; change can be expressed as delta_h0 / sqrt (max(err^2)-min(err^2))
   ; The reference variant to compare to is in refvar_array[k]; if
   ; 'none' or '', no significance is computed.   
   sigdiff = 0.d0
   sigdiffstr = '  ---  '
   if (refvar_array[k] ne 'none' and refvar_array[k] ne '') then begin
      whref = where (id_array eq refvar_array[k], nref)
      if (nref eq 1) then begin
         delta = dd.h0_value - det_list[whref[0]].h0_value
         ea = dd.h0_error
         eb = det_list[whref[0]].h0_error
         sigdiff = delta / sqrt (max([ea^2,eb^2])-min([ea^2,eb^2]))
         sigdiffstr = string(sigdiff, format='(1x,f5.2)')+' '
      endif
   endif
   pte = pte_chi2 (dd.chi2[0], dd.ndof[0])
   if (pte gt 0.05) then begin
      pte_str = string(pte, format='(1x,f5.3)')+' '
   endif else begin
      pteconv = pte
      ndex = 0
      while (pteconv lt 1) do begin
         pteconv = pteconv * 10
         ndex = ndex+1
      endwhile
      pte_str = string(pteconv, 'e-', ndex, format='(1x,f3.1,a,i0)')
   endelse
   ;
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
   printf, lun, k+1, dd.h0_value, dd.h0_error, dd.chi2[0], dd.ndof, $
           dd.chi2[0]/dd.ndof, pte_str, sigdiffstr, $
           absmag_sn1a, absmag_sn1a, alpha_sn1a, $
           ealpha_sn1a, n_sn1a, $
           absmag_sn2, eabsmag_sn2, alpha_sn2, ealpha_sn2, n_sn2, label, $
           format='(i2, ") ", 2f10.4, 2x, f10.4, i7, f10.4, 2a, 3x, ' + $
           '2(2f10.3, 2f10.4, i4), 2x, a)'
   ; stop
   print, k+1, dd.h0_value, dd.h0_error, dd.chi2[0], dd.ndof, $
          dd.chi2[0]/dd.ndof, pte_str, sigdiffstr, $
          absmag_sn1a, eabsmag_sn1a, alpha_sn1a, $
          ealpha_sn1a, n_sn1a, absmag_sn2, eabsmag_sn2, $
          alpha_sn2, ealpha_sn2, n_sn2, label, $
          format='(i2, ") ", 2f10.4, 2x, f10.4, i7, f10.4, 2a, 3x, ' + $
          '2(2f10.3, 2f10.4, i4), 2x, a)'
   ; stop
   ; For LaTeX
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
           dd.h0_error, ' & ', dd.ndof, ' & ', $
           dd.chi2[0]/dd.ndof, ' & ', pte_str, ' & ', sigdiffstr, ' & ', $
           label2, '\\', format = '(a5, a3, f6.3, a3, f5.3, a3, i3, ' + $
           ' a3, f6.4, 2(a3, a7), a3, a, a3)'

   ; Old print statement; includes chi2 and some SN1a parameters,
   ; does not include PTE, sigma_diff
   ; if (inlatex[k]) then printf, lun2, label1, ' & ', dd.h0_value, ' & ', $
           ; dd.h0_error, ' & ', dd.ndof, ' & ', dd.chi2[0]/dd.ndof, $
           ; ' & ', s_absmag, ' & ', s_eabsmag, ' & ', $
           ; s_alpha, ' & ', s_ealpha, ' & ', $
           ; n_sn1a, ' & ', label2, ' \\', $
           ; format='(a5, a3, f6.3, a3, f5.3, a3, f8.4, a3, i3, a3, f7.5, ' + $
           ; ' 4(a3, a7), 2(a3, a8), a3, i3, a3, a, a3)'
           ; dd.h0_error, ' & ', dd.chi2[0], ' & ', dd.ndof, ' & ', $
           ; dd.chi2[0]/dd.ndof, ' & ', pte_str, ' & ', sigdiffstr, $
           ; ' & ', s_absmag, ' & ', s_eabsmag, ' & ', $
           ; s_alpha, ' & ', s_ealpha, ' & ', $
           ; n_sn1a, ' & ', label2, ' \\', $
           ; format='(a5, a3, f6.3, a3, f5.3, a3, f8.4, a3, i3, a3, f7.5, ' + $
           ; ' 4(a3, a7), 2(a3, a8), a3, i3, a3, a, a3)'
   ; LaTeX variable definitions.  Use 2 decimal places for H0 and error
   ; stop
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

printf, lun2, '\noalign{\smallskip}'
printf, lun2, '\hline'
printf, lun2, '\noalign{\smallskip}'
printf, lun2, 'O1 & 73.110 & 0.920 & 57 & 0.9360 & & & Orth. path 1 MW+LMC/SMC+Ceph+SNIa+FP \\' 
printf, lun2, 'O2 & 73.451 & 1.777 & 23 & 0.4853 & & & Orth. path 2 N4258+TRGB+SBF+MM (V00 equivalent) \\' 
printf, lun2, 'O2\_V99a & 74.083 & 1.249 & 27 & 0.6016 & & & Orth. path 2 O2+SNII+EPM (V99a equivalent) \\'
printf, lun2, 'O2\_V99 & 74.780 & 1.133 & 151 & 1.4613 & & & Orth. path 2 O2+SNII+EPM+TF(V99 equivalent) \\'
printf, lun2, '\noalign{\smallskip}'
printf, lun2, '\hline\hline'
printf, lun2, '\end{tabular}'
printf, lun2, '\tablefoot{{\\'
printf, lun2, 'a) PTE, or Probability to Exceed, is the probability'
printf, lun2, 'of obtaining a value of $ \chi^2 $ larger than measured for that variant.\\'
printf, lun2, 'b) $ \Delta_\sigma $ is the significance of the change in {\Hcst}'
printf, lun2, 'for the current variant in units of the relative uncertainty, determined as the'
printf, lun2, 'quadrature difference of the total uncertainties.  See text for details.\\}'
printf, lun2, '\end{table}'

close, lun
free_lun, lun
close, lun2
free_lun, lun2
close, lun3
free_lun, lun3

print_extrapars, det_list, descr_list, latexlabel=latexlabel, $
                 inlatex=inlatex, outdir=outdir

return
end
