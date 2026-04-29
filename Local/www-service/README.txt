
Contents:

[1]
config/*
httpd.conf files for apache httpd server

[2]
var/*
Log and process ID (pid) files for running HTTP servers.

For instance, the files for the (default) server on port 8088 are:

  http-mustang2.8090.log       -- server log file
  http-mustang2.8090.pid       -- process ID of server

Old logs for servers that were run earlier, but have now exited,
are kept around and truncated when a new server on that port
is run again.

For more, see util/html-serve.py.
