Pro make_host_distance_table, dd_base=dd_base, dd_all=dd_all, $
                              outfile=outfile, debug=debug, order=order
;
; Extract the hosts result information from the two details
; structures.  Output a LaTeX table with hosts and their distance
; moduli and errors.
;
; if ORDER, reorder output in correct order without regard to
; data availability
  
h_base = dd_base.host_data_details
h_all = dd_all.host_data_details

host_base = h_base.hosts
nbase = h_base.nhosts

host_all = h_all.hosts
nall = h_all.nhosts

; Edit the host name and order to be more publication-worthy
; Assume:
; N for NGC
; M for Messier - but check M1337
; M1337 = Mrk 1337
; P for ?
; U for UGC
; UA for UGCA
catalog = strarr(nall)
catkey = strarr(nall)
idstring = strarr(nall)
idnumber = lonarr(nall)
lexall = strarr(nall)
host_form = strarr(nall)

for k = 0, nall-1 do begin
   catkey[k] = strmid(host_all[k],0,1)
   idstring[k] = strmid(host_all[k],1,strlen(host_all[k])-1)
   if (host_all[k] eq 'M1337') then catkey[k] = 'Mrk'
   if (strmid(host_all[k],0,2) eq 'UA') then begin
      catkey[k] = 'UA'
      idstring[k] = strmid(host_all[k],2,strlen(host_all[k])-2)
   endif
   if (strmid(host_all[k],0,2) eq 'IC') then begin
      catkey[k] = 'IC'
      idstring[k] = strmid(host_all[k],2,strlen(host_all[k])-2)
   endif
endfor
; stop
idnumber = long(idstring)
for k = 0, nall-1 do begin
   case catkey[k] of
      'IC' : catalog[k] = 'IC'
      'M': catalog[k] = 'M'
      'N': catalog[k] = 'NGC'
      'P' : catalog[k] = 'PGC'
      'Mrk': catalog[k] = 'Mrk'
      'U': catalog[k] = 'UGC'
      'UA' : catalog[k] = 'UGCA'
      else: catalog[k] = catkey[k]
   endcase
   lexall[k] = string (catalog[k], idnumber[k], format='(a,i09)')
   host_form[k] = string (catalog[k], idnumber[k], format='(a, 1x, i0)')
endfor

; M1337 from R22
; UA319 from A21 = UGCA 391
; P* from A21



ibase = replicate (-1, nall)
for k = 0, nall-1 do begin
   wh = where (host_base eq host_all[k], nwh)
   if (nwh eq 1) then ibase[k] = wh[0]
endfor

whboth = where (ibase ge 0, nboth)
whonly = where (ibase lt 0, nonly)

muall = dd_all.params_value
emuall = dblarr (nall)
for k = 0,nall-1 do emuall[k] = sqrt(dd_all.params_var[k,k])

mubase = dd_base.params_value
emubase = dblarr (nbase)
for k = 0,nbase-1 do emubase[k] = sqrt(dd_base.params_var[k,k])


if (keyword_set(outfile) eq 0) then outfile = 'distance_table.tex'
openw, lun, outfile, /get_lun

if (keyword_set(order)) then begin
   iorder = sort(lexall)
   for k = 0, nall-1 do begin
      i = iorder[k]
      j = ibase[i]
      if (j ne -1) then begin
         printf, lun, host_form[i], mubase[j], emubase[j], $
                 muall[i], emuall[i], $
                 format='(a14, 4x, 2(" & ",f8.3), 4x, 2(" & ",f8.3), " \\")'
      endif else begin          ; no baseline value
         printf, lun, host_form[i], '   ---  ', '   ---  ', $
                 muall[i], emuall[i], $
                 format='(a14, 4x, 2(" & ",a8), 4x, 2(" & ",f8.3), " \\")'
      endelse
   endfor
endif else begin
   for k = 0, nboth-1 do begin
      i = whboth[k]
      j = ibase[i]
      printf, lun, host_all[i], mubase[j], emubase[j], muall[i], emuall[i], $
          format='(a14, 4x, 2(" & ",f8.3), 4x, 2(" & ",f8.3), " \\")'
   endfor
   ;
   for k = 0, nonly-1 do begin
      i = whonly[k]
      printf, lun, host_all[i], '   ---  ', '   ---  ', muall[i], emuall[i], $
              format='(a14, 4x, 2(" & ",a8), 4x, 2(" & ",f8.3), " \\")'
   endfor
endelse

close, lun
free_lun, lun

if (keyword_set(debug)) then stop
return
end
