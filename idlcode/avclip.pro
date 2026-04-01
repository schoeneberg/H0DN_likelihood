function avclip, array, disp=disp, accept=accept, niter=niter, kfactor=kfactor
;
; Simple average + sigma clip routine.  Does NITER [3] iterations with
; KFACTOR [3] -sigma clip.  Returns average as function value,
; dispersion and "pass" array as arguments.
;

if (keyword_set(niter) eq 0) then niter = 3
if (keyword_set(kfactor) eq 0) then kfactor = 3.0

mean = avg (array)
disp = sigma (array)

for i = 0, niter-1 do begin
   accept = ( abs(array - mean) le kfactor*disp )
   ww = where (accept)
   mean = avg (array[ww])
   disp = sigma (array[ww])
endfor

return, mean
end
