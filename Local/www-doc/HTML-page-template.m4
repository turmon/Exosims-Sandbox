dnl
dnl Webpage template as an M4 document
dnl macros below are defined in the command-line 
dnl of the m4 invocation
dnl
_startpage(META_Title)
<!-- begin included content -->
dnl undivert() is a gnu extension that does
dnl not process the named file with m4
dnl (otherwise, occurrence of m4 constructs
dnl in HTML_CONTENT would be interpreted)
undivert(HTML_CONTENT)
<!-- end included content -->
_endpage
