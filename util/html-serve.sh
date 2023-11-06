#!/usr/bin/env bash
#
# html-serve.sh: start/stop/status HTTP server
#
# This starts, stops, or status-queries the HTTP server that 
# serves up the EXOSIMS table of contents and ensemble-summary files.
#
# For security, the server only serves out to connections coming
# from localhost, i.e., the system the server is running on.
#
# We run the server on a non-standard port - 8090 by default.
# To see the table of contents, navigate to:
#   http://localhost:8090/sims/
# To connect your own system to localhost, use an ssh tunnel:
#   ssh -fnNT -L8090:localhost:8090 mustang2.jpl.nasa.gov
# and visit the above URL on your own system (e.g., laptop).
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
#   Local/www-service/var/http-HOST.PORT.*
# where PORT is typically 8090.  (Formerly 8088.)
# There is a server error log (.log) and a file containing the process ID (PID)
# of the server (.pid).
# The PID is used to signal the running server to exit.  If the basic "stop" 
# MODE of this command does not work, you can kill the server by sending TERM
# to process ID (pid) named in that file.
# 
# turmon sep 2018, dec 2021, jan 2022
#
## [end comment block]

# exit-on-error
set -euo pipefail

PROGNAME=$(basename $0)

# attempt to give group-write to created files
umask 002

# runtime context (s383/JPL supercomputer)
if [ -d /proj/exep ]; then
    CONTEXT=s383
elif [ -d /projects/exo_yield ]; then
    CONTEXT=jplsc
else
    echo "${PROGNAME}: Fatal: Could not determine runtime context." >&2
    exit 1
fi

## defaults
# SERVER_GROUP: whole group is allowed to kill the server
# DOC_ROOT: httpd needs to know the rooted dir of the content
if [ $CONTEXT = s383 ]; then
    SERVER_GROUP=exosims
    DOC_ROOT=/proj/exep/rhonda/Sandbox/HabEx
elif [ $CONTEXT = jplsc ]; then
    SERVER_GROUP=exo-yield
    DOC_ROOT=/projects/exo_yield/Sandbox/hwo
fi
# httpd config file
SERVER_CONFIG=Local/www-service/config/httpd-apache2.4.conf
# directory for log and PID files produced by httpd
SERVER_VARDIR=Local/www-service/var

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
   echo "${PROGNAME}: Error: Need exactly one argument, use -h for help." >&2
   exit 1
fi

mode="$1"

host=$(hostname -s)

# specific filenames, now that we know $port
SERVER_LOG=$SERVER_VARDIR/http-$host.$port.log
SERVER_PID=$SERVER_VARDIR/http-$host.$port.pid


# set up the server state files
function setup {
    # copy and reset the server log
    [ -f "$SERVER_LOG" ] && cp -pf "$SERVER_LOG" "$SERVER_LOG".old
    echo > "$SERVER_LOG"
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
	kill -TERM $(cat $SERVER_PID) || true
	sleep 0.1
	# if the PID file was not removed, the server is apparently not 
        # responding ... proceed anyway (e.g., host machine restart)
	if [ -r "$SERVER_PID" ]; then
	    echo "${PROGNAME}: Server process not responding or not present." >&2
	    echo "${PROGNAME}: Proceeding with restart anyway." >&2
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
	# 1: formerly: "ServerName $(hostname)", "Listen $port", but
	# now only listen for connections from localhost - must make an
	# ssh tunnel from your remote->localhost to get a server connection
	# 2: for use across s383/jpl-sc, we export a variable for
	# httpd.conf ($SERVER_CONFIG) to pick up
	DOCUMENT_ROOT="$DOC_ROOT" \
	httpd -f "$CURR_DIR/$SERVER_CONFIG" \
	      -c "PidFile $CURR_DIR/$SERVER_PID" \
	      -c "ErrorLog $CURR_DIR/$SERVER_LOG" \
	      -c "ServerName localhost" \
	      -c "Listen localhost:$port"
    else
	# python SimpleHTTPServer: single-threaded, does not support returning
	# byte ranges, which makes .mp4s not render in Safari.
	# mostly-daemonized: nohup, stdout to /dev/null, sterrr to logfile
	# (12/2021: this has not been tested recently)
	nohup python -m http.server "$port" 1>/dev/null 2>>$SERVER_LOG &
	# capture server PID - last backgrounded command
	server_pid_var=$!
    fi
    # status update to stdout
    [ $port != $DEFAULT_PORT ] && p_text=" -p $port " || p_text=" "
    echo "${PROGNAME}: Started the server.  Stop with: \`${PROGNAME}${p_text}stop'" 

elif [ "$mode" == stop ]; then
    #
    # stop the HTTP server
    #
    if [ ! -r $SERVER_PID ]; then
	echo "${PROGNAME}: No server is running on $host on port number $port." >&2
	echo "${PROGNAME}: Try \`${PROGNAME} status', or check $SERVER_VARDIR for another port number?" >&2
	exit 1
    fi
    echo "${PROGNAME}: Killing server named in $SERVER_PID."
    kill -TERM $(cat $SERVER_PID)

elif [ "$mode" == status ]; then
    #
    # which servers are running?
    #
    filepat="$SERVER_VARDIR/http-*.pid"
    if ls $filepat &> /dev/null; then
	echo "Looking in: $SERVER_VARDIR"
	echo "Found server(s):" 
	for pfile in $filepat; do
	    echo "$pfile"
	    p_info=$(stat -c 'started by %U on %y' $pfile)
	    x_hopo=$(echo $pfile | sed -e 's|.*/http-||' -e 's|\.pid||' -e 's|\.|:|')
	    x_host=$(echo $x_hopo | sed -e 's|:.*||')
	    x_port=$(echo $x_hopo | sed -e 's|.*:||')
	    # we only listen on localhost since 2023-October
	    # (the above filename manipulations find the host the server is running on)
	    x_host=localhost
	    echo "  " "HOST:PORT = ${x_host}:${x_port}"
	    echo "  " "$p_info"
	    echo "  " "View at this URL: http://${x_host}:${x_port}/sims"
	    if wget -q -t 1 -O /dev/null "http://${x_host}:${x_port}/sims/"; then
		echo "  " "Server seems to respond OK to requests over HTTP."
	    else
		echo "  " "Server does not respond to requests over HTTP."
	    fi
	done
    else
	echo "No HTTP servers running." 
    fi
fi

