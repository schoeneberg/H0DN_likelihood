Function mpinv, a, tol=tol, rank=rank
;
; Compute the Moore-Penrose inverse M of A using singular value decomposition.
;
; M is defined by the Moore-Penrose conditions:
; A M A = A
; M A M = M
; (A M)^* = A M
; (M A)^* = M A
;
; M can be computed by singular value decomposition.  If
; A = U S V^*
; with U = m*m orthogonal matrix, V = n*n orthogonal matrix, and
; S = m*n diagonal materix with non-negative real elements (A is m*n), then
;
; M = V S+ U^*
; where S+ is the pseudoinverse of S: a diagonal n*m matrix whose diagonal
; elements are the reciprocal of those of S, or 0 if the corresponding
; element of S is 0.
;
; Numerically, elements of S that are close to 0 need to be set to 0.
; The TOL parameter defines the range of elements that are et to 0;
; it defaults to epsilon * (m > n) * max (S), where epsilon is the
; machine precision and thedefault value is the same used in the
; function pinv in GNU Octave.
;

sz = size(a)
if (sz[0] ne 2) then begin
   print, ' Error: argument must be 2-dimensional array'
   return, 0
endif
m = sz[1]
n = sz[2]
is_complex = 0B
is_double = 0B
itype = size(a, /type)
if (itype eq 6 or itype eq 9) then is_complex = 1B
if (itype eq 5 or itype eq 9) then is_double = 1B

la_svd, a, w, u, v

if (is_double eq 0) then q = machar () else q = machar(/double)
epsilon = q.eps

if (keyword_set(tol) eq 0) then tol = max(w) * (m > n) * epsilon

whsmall = where (abs(w) le tol, nsmall)
whlarge = where (abs(w) gt tol, nlarge)
winv = w
if (nsmall gt 0) then winv[whsmall] = 0.d0
if (nlarge gt 0) then winv[whlarge] = 1.d0 / w[whlarge]

if (is_complex) then uu = conj(transpose(u)) else uu = transpose(u)
result = v ## diag_matrix(winv) ## uu

rank = nlarge

return, result
end
