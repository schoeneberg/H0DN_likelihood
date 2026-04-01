Pro chg_default_dir, filename, basedir, notrim=notrim

; Adds the default directory BASEDIR to FILENAME unless FILENAME already
; has a directory designation (established by searching for
; a forward slash in the filename).  Do not add the default directory
; if FILENAME has the special values 'none' or '' (after changing to
; lowercase and trimming leading and trailing blanks).
;
; Leading blanks on filename are trimmed on output
; (to avoid embedded blanks) unless NOTRIM is set.

has_dir = strpos(filename, '/') ge 0

fcompare = strlowcase(strtrim(filename, 2))
if (keyword_set(notrim) eq 0) then notrim = 0
if (not notrim) then filename = strtrim(filename, 1)
if (not has_dir and fcompare ne 'none' and fcompare ne '') then $
    filename = basedir + filename

return
end
