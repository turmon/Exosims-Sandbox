dnl macros for HTML header and footer
changequote(`{{',`}}')dnl
define({{_startpage}}, {{<html>
  <!-- This HTML file was auto-generated by a macro system -->
  <!-- To update, please change the original markdown (.md) file -->
  <!-- and re-run the templating system. -->
  <head>
    <title>$1</title>
    <link rel="stylesheet" href="ensemble-reports.css">
  </head>
  <body>
}})dnl
dnl
define({{_endpage}}, {{
  </body>
</html>}})dnl
dnl