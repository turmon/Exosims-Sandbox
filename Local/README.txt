Information that supports EXOSIMS for one specific mission application.

This was intended to consist only of files that were customized for a 
single mission.  It was thought that most plots would be mission-specific, 
and that some EXOSIMS source code specializations would pertain to only
one mission.  

This did not prove to be an important distinction, so many of these
files could go in ../util as well.

Files here:

EXOSIMS Drivers
  ipcluster_ensemble_jpl_old.py -- obsolete version of the below
  ipcluster_ensemble_jpl_driver.py -- python driver for Exosims runs in Sandbox
  EXOSIMS_local/* -- classes supporting DRM output to Sandbox structures

Support for EXOSIMS runs under ipython-parallel
  00-path.py
  10-threads.py

Matlab-based plots
  Matlab/*

Web server for results viewing
  www-resources/  -- CSS and JS and other static pieces supporting 
                     the results website, and accessed by httpd.
  www-service/    -- Invocation and logs for the web-server, apache httpd.
