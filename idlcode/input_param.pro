Function input_param, name, default, plist, pvalues, was_found
;
; Extract the parameter whose name in PLIST matches NAME from
; the array PVALUES.  If not found, return DEFAULT.  Parameter is
; not modified; if PVALUES is a string array, then it needs to be
; converted in the calling routine.
; WAS_FOUND is set to 1 if the parameter was actually found in PLIST,
; or 0 if it was set to the default.
; Note that this routine will return arrays if multiple elements in
; PLIST match NAME.  This may not always be what is desired; an
; earlier version actually returned an error if there were multiple
; matches.  But in this way, short arrays can be returned.
;

wh = where (strmatch (plist, name, /fold_case) ne 0, nwh)
if (nwh eq 0) then begin
   result = default
   was_found = 0B
endif else begin
   if (nwh eq 1) then begin
      result = pvalues[wh[0]]
      was_found = 1B
   endif else begin
      result = pvalues[wh]
      ; print, ' Error: duplicate parameter found for ', name
      ; stop
   endelse
endelse

return, result
end
