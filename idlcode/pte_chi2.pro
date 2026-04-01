Function pte_chi2, value, df, verbose=verbose
;
; Compute the probability to exceed the given value of CHI2
; for a chi2 distribution with DF degrees of freedom.
;
; The chi2 distribution is
; f(chi2, df) = 1/(2^df/2) 1/gamma(df/2) sqrt(chi2)^(df-2) exp(-chi2/2)
;
; and the probability to exceed the value z is
; pte(z, df) = int_z^infty f(t,df) dt
;
; The cdf (integral from 0 to z) is computed as:
;
; int_0^z f(t, df) dt = ligamma(df/2. z/2) / gamma(df/2)
; where ligamma is the lower incomplete gamma function.
;
; IDL implements the *regularized* incomplete gamma function, defined as
; igamma(a,z) = int_0^z exp(-t) t^(a-1) dt / int_0^infty exp(-t) t^(a-1) dt
;
; Thus the PTE value is the complement to 1 of igamma(df/2, z/2)
;

; Note that this method of computation can be inaccurate for large values,
; since the precision of the igamma approximation may not be optimal and
; cancellation can occur.  A direct computation of the integral may be
; preferable for large values of the first argument.
  
result = 1.d0 - igamma (df/2, value/2)

if (keyword_set(verbose)) then print, gamma(df/2), igamma(df/2, value/2) 

return, result
end
