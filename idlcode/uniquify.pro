Function uniquify, input_name, index, sort_index=sindex, robust=robust
;
; Returns a subarray containing the unique elements of input_name in
; sorted order.  Optionally INDEX returns the index of each element of
; input_name in the uniquified array.  Thus result[index[i]] equals
; input_name[i].
;

temp = sort(input_name)
name_temp = input_name[temp]
utemp = uniq(name_temp)
uname_temp = name_temp[utemp]

iii = lindgen (n_elements(temp))
sindex = iii[sort(temp)]  ; sindex[i] is the position of input_name[i] in the sorted array
                          ; for example, sindex[0] is the order
                          ; position of the first element of input_name
                          ; i.e., name_temp[sindex[i]] equals input_name[i]
; index requires a loop; only execute if requested
if (n_params() gt 1) then begin
   index = lonarr(n_elements(input_name))
   if (keyword_set(robust) eq 0) then begin
      for j = 0L, n_elements(utemp)-2L do index[utemp[j]+1:*] = index[utemp[j]+1:*]+1  ; this gives the element of uname_temp
      index = index[sindex]
   endif else begin
      for i = 0L, n_elements(uname_temp)-1L do begin
         wh = where (input_name eq uname_temp[i], nwh)
         if (nwh gt 0) then index[wh] = i
      endfor
   endelse
   ; The regular version increments all
   ; elements of index after each element
   ; of uniq(name_temp).  This index
   ; array needs to be mapped into the
   ; original elements of input_name
   ; using sindex.  Fo clarity, the index computed
   ; in the for loop is such that
   ; uname_temp[index[i]] equals name_temp[i]
   ; thus  uname_temp[index[sindex[i]]]
   ;       equals name_temp[sindex[i]] equals
   ;       input_name[i]
   ; An earlier version incorrectly assumed that UNIQ returned the FIRST
   ; index for equal elements.  In fact
   ; it returns the LAST such index.
   ; The robust version just uses where
   ; to determine the mapping.  This is
   ; much cleaner if potentially slower.
   ;
   ; Test using a 100K element array with 63K distinct elements:
   ;    14s for the "normal" version
   ;    12s for the "robust" version
   ;    results are identical
   ; stop
endif
return, uname_temp
end
