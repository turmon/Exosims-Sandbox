#!/usr/bin/env bash
#
# html-serve.sh: start/stop HTTP server
#
# This starts or stops the HTTP server that serves up the Exosims
# table of contents and ensemble-summary files.
#
# We run the server on a non-standard port - 8090 by default.
# To see the table of contents, navigate to:
#   http://mustang2.jpl.nasa.gov:8090/sims/
#
# Usage:
#   html-serve.sh [-h] [-p PORT] MODE
#
# where MODE is one of:
#   start: start the server on the given port
#   stop: stop the server on the given port
#   status: list the servers running, if any
#   init: check that server status files work: should be unneeded.
#
# and:
#  -p PORT => gives the HTTP port number.  Default is 8090.
#  -h      => shows this help text.
#
# Implementation:
# -- Status files are kept in these files:
#   Local/www-service/http-HOST.PORT.*
# where PORT is typically 8090.  (Formerly 8088.)
# There is a server error log (.log) and a file containing the process ID (PID)
# of the server (.pid).
# The PID is used to signal the running server to exit.  If the basic "stop" 
# MODE of this command does not work, you can kill the server by sending TERM
# to process ID (pid) named in that file.
# 
# turmon sep 2018, dec 2021
#
## [end comment block]

# exit-on-error
set -euo pipefail

PROGNAME=$(basename $0)

# attempt to give group-write to created files
umask 002

# some defaults
#  group allowed to kill the server
SERVER_GROUP=exosims
#  directory for log files
SERVER_DIR=Local/www-service
# current working directory - httpd needs absolute pathnames
CURR_DIR=$(pwd)
# default port -- formerly 8088
DEFAULT_PORT=8090


port=$DEFAULT_PORT
while getopts "hp:" opt; do
    case $opt in
	p)
	    # port number
	    port="$OPTARG"
            echo "${PROGNAME}: Operating on port $port"
	    ;;
	h)
	    # help text
	    sed 's/^#//' $(which $0) | awk '/^#/{exit};NR>1{print}'
	    exit 2
	    ;;
	\?)
	    echo "${PROGNAME}: Invalid option, exiting.  Try -h." >&2
	    exit 2
	    ;;
    esac
done
shift $((OPTIND-1))

# check port number
if ! [[ $port =~ ^[0-9]+$ ]] ; then
   echo "${PROGNAME}: Error: Supplied port must be a number." >&2
   exit 1
elif (( $port < 1024 )); then
   echo "${PROGNAME}: Error: Supplied port must be 1024 or greater." >&2
   exit 1
fi

# enforce 1 argument
if [ $# -ne 1 ]; then
   echo "${PROGNAME}: Error: Need exactly one argument" >&2
   exit 1
fi

mode="$1"

host=$(hostname -s)

# specific filenames we need
SERVER_LOG=$SERVER_DIR/http-$host.$port.log
SERVER_PID=$SERVER_DIR/http-$host.$port.pid


# set up the server state dir
function setup {
    # (assume it is already created and has right mode)
    ## mkdir -m 775 -p "$SERVER_DIR"
    ## chgrp $SERVER_GROUP "$SERVER_DIR"
    # copy and reset the server log
    [ -f "$SERVER_LOG" ] && cp -pf "$SERVER_LOG" "$SERVER_LOG".old
    rm -f "$SERVER_LOG" && touch "$SERVER_LOG"
    echo "${PROGNAME}: HTTP file setup OK." >> "$SERVER_LOG"
    chgrp $SERVER_GROUP "$SERVER_LOG"
    chmod 664 "$SERVER_LOG"
    # remove server PID file, if any
    rm -f "$SERVER_PID"
}


if [ "$mode" == start ]; then
    #
    # start the HTTP server
    #
    # Try not to have two going on same machine + port!
    if [ -r $SERVER_PID ]; then
	echo "${PROGNAME}: Server may be running on $port already." >&2
	echo "${PROGNAME}: Attempting to shut down cleanly and restart." >&2
	kill -TERM $(cat $SERVER_PID)
	sleep 0.1
	# if the PID file was removed, that's a good sign
	if [ -r "$SERVER_PID" ]; then
	    echo "${PROGNAME}: Could not shut down cleanly." >&2
	else
	    echo "${PROGNAME}: Shut down cleanly.  Restarting." >&2
	fi
    fi
    # reset the file status
    setup
    # start the server and grab its PID so we can kill it later
    if true; then
	# apache httpd - multi-threaded, serves movies properly
	# minimal config, PID and logging configured on command line
	httpd -f "$CURR_DIR/$SERVER_DIR/config/httpd-apache2.4.conf" \
	      -c "PidFile $CURR_DIR/$SERVER_PID" \
	      -c "ErrorLog $CURR_DIR/$SERVER_LOG" \
	      -c "Listen $port"
	# capture server PID -- httpd forks, so the pid changes
	server_pid_var=$(sleep 0.1 && cat "$CURR_DIR/$SERVER_PID")
    else
	# python SimpleHTTPServer: single-threaded, does not support returning
	# byte ranges, which makes .mp4s not render in Safari.
	# mostly-daemonized: nohup, stdout to /dev/null, sterrr to logfile
	# (12/2021: this has not been tested recently)
	nohup python -m SimpleHTTPServer "$port" 1>/dev/null 2>>$SERVER_LOG &
	# capture server PID - last backgrounded command
	server_pid_var=$!
    fi
    # status update to stdout
    [ $port != $DEFAULT_PORT ] && p_text="-p $port" || p_text=""
    echo "${PROGNAME}: Started the server.  Stop with: \`${PROGNAME} $p_text stop'" 

elif [ "$mode" == stop ]; then
    #
    # stop the HTTP server
    #
    if [ ! -r $SERVER_PID ]; then
	echo "${PROGNAME}: No server is running on $host on port number $port." >&2
	echo "${PROGNAME}: Try \`${PROGNAME} status', or check $SERVER_DIR for another port number?" >&2
	exit 1
    fi
    echo "${PROGNAME}: Killing server at $SERVER_PID."
    kill -TERM $(cat $SERVER_PID)

elif [ "$mode" == status ]; then
    #
    # which servers are running?
    #
    filepat="$SERVER_DIR/http-*.pid"
    if ls $filepat &> /dev/null; then
	echo "Found server(s):" 
	for pfile in $filepat; do
	    echo "$pfile"
	    p_info=$(stat -c 'started by %U on %y' $pfile)
	    x_hopo=$(echo $pfile | sed -e 's|.*/http-||' -e 's|\.pid||' -e 's|\.|:|')
	    x_host=$(echo $x_hopo | sed -e 's|:.*||')
	    x_port=$(echo $x_hopo | sed -e 's|.*:||')
	    echo "  " "HOST:PORT = $x_hopo"
	    echo "  " "$p_info"
	    echo "  " "View at this URL: http://${x_host}.jpl.nasa.gov:${x_port}/sims"
	    if wget -q -t 1 -O /dev/null "http://${x_host}.jpl.nasa.gov:${x_port}/sims/"; then
		echo "  " "Server seems to respond OK."
	    else
		echo "  " "Server does not respond."
	    fi
	done
    else
	echo "No HTTP servers running." 
    fi
fi

