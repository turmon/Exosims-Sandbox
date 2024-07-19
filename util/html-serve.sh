#!/usr/bin/env bash
#
# html-serve.sh: start/stop/status HTTP server
#
# This starts, stops, or status-queries the HTTP server that 
# serves up the EXOSIMS index and ensemble-summary files.
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
#   html-serve.sh [-h] [-p PORT] [-s apache|simple] MODE
#
# where MODE is one of:
#   start: start the server on the given port
#   stop: stop the server on the given port
#   status: list the servers running, if any
#   init: check that server status files work: should be unneeded.
#
# and:
#  -p PORT   => gives the HTTP port number.  Default is 8090.
#  -s apache => use apache2/httpd server (default)
#  -s simple => use python SimpleHTTPServer (less performant)
#  -h        => shows this help text.
#
# Implementation:
# -- Status files are kept in these files:
#    Local/www-service/var/http-HOST.PORT.*
# where PORT is typically 8090.
# There is a server error log (.log) and a file containing the process ID (PID)
# of the server (.pid).
# The PID is used to signal the running server to exit.  If the basic "stop" 
# MODE of this command does not work, you can kill the server by sending TERM
# to process ID (pid) named in that file. "ps ufx" can also be helpful to
# get status of the server (if you're the process owner).
# 
## [end comment block]

# turmon sep 2018, dec 2021, jan 2022, jul 2024

# exit-on-error
set -euo pipefail

PROGNAME=$(basename "$0")

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
# SERVER: apache server binary
# SERVER_GROUP: whole group is allowed to kill the server
# DOC_ROOT: httpd needs to know the rooted dir of the content
if [[ $CONTEXT = s383 ]]; then
    SERVER=apache2
    SERVER_GROUP=exosims
    DOC_ROOT=/proj/exep/rhonda/Sandbox/HabEx
    # (need the slash)
    MOD_ROOT=/usr/lib/apache2/
elif [ $CONTEXT = jplsc ]; then
    SERVER=httpd
    SERVER_GROUP=exo-yield
    DOC_ROOT=/projects/exo_yield/Sandbox/hwo
    # (empty)
    MOD_ROOT=
fi
# httpd config file
SERVER_CONFIG=Local/www-service/config/httpd-apache2.4.conf
# directory for log and PID files produced by httpd
SERVER_VARDIR=Local/www-service/var

# current working directory - httpd needs absolute pathnames
CURR_DIR=$(pwd)
# default port
DEFAULT_PORT=8090
# default server
DEFAULT_SERVER=apache


port=$DEFAULT_PORT
server=$DEFAULT_SERVER
while getopts "hp:s:" opt; do
    case $opt in
	s)
	    # server flavor
	    server="$OPTARG"
            echo "${PROGNAME}: Server type set to: $server"
	    ;;
	p)
	    # port number
	    port="$OPTARG"
            echo "${PROGNAME}: Operating port set to: $port"
	    ;;
	h)
	    # help text
	    sed 's/^# \?//' "$(which "$0")" | awk '/^#/{exit};NR>1{print}'
	    exit 2
	    ;;
	\?)
	    echo "${PROGNAME}: Invalid option, exiting.  Try -h." >&2
	    exit 2
	    ;;
    esac
done
shift $((OPTIND-1))

# enforce 1 argument
if [ $# -ne 1 ]; then
   echo "${PROGNAME}: Error: Need exactly one argument, use -h for help." >&2
   exit 1
fi

# extract the mode
mode="$1"

# check port number
if ! [[ $port =~ ^[0-9]+$ ]] ; then
   echo "${PROGNAME}: Error: Supplied port must be a number." >&2
   exit 1
elif (( $port < 1024 )); then
   echo "${PROGNAME}: Error: Supplied port must be 1024 or greater." >&2
   exit 1
fi

# check server, message only for "start" mode
if [[ "$server" == simple ]] ; then
    if [[ "$mode" == start ]]; then
	echo "${PROGNAME}: Using server: $server."
    fi
elif [[ "$server" == apache ]] ; then
    if [[ "$mode" == start ]]; then
	echo "${PROGNAME}: Using server: $server."
    fi
else
   echo "${PROGNAME}: Error: Supplied server ($server) not recognized." >&2
   exit 1
fi

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
    chgrp "$SERVER_GROUP" "$SERVER_LOG"
    chmod 664 "$SERVER_LOG"
    # remove server PID file, if any
    rm -f "$SERVER_PID"
}


if [[ "$mode" == start ]]; then
    #
    # start the HTTP server
    #
    # Try not to have two going on same machine + port!
    if [ -r "$SERVER_PID" ]; then
	echo "${PROGNAME}: Server may be running on $port already." >&2
	echo "${PROGNAME}: Attempting to shut down cleanly and restart." >&2
	kill -TERM "$(cat "$SERVER_PID")" || true
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
    if [[ "$server" == "apache" ]]; then
	# apache httpd - multi-threaded, serves movies properly
	# minimal config, PID and logging configured on command line
	# 1: formerly: "ServerName $(hostname)", "Listen $port", but
	# now only listen for connections from localhost - must make an
	# ssh tunnel from your remote->localhost to get a server connection
	# 2: for use across s383/jpl-sc, we export two variables for
	# httpd.conf ($SERVER_CONFIG) to pick up: the first 2 lines below.
	MODULE_ROOT="$MOD_ROOT" \
	DOCUMENT_ROOT="$DOC_ROOT" \
	$SERVER -f "$CURR_DIR/$SERVER_CONFIG" \
	      -c "DefaultRuntimeDir $(pwd)" \
	      -c "PidFile $CURR_DIR/$SERVER_PID" \
	      -c "ErrorLog $CURR_DIR/$SERVER_LOG" \
	      -c "ServerName localhost" \
	      -c "Listen localhost:$port"
    elif [[ "$server" == "simple" ]]; then
	# python SimpleHTTPServer: single-threaded, does not support returning
	# byte ranges, which makes .mp4s not render in Safari.
	# mostly-daemonized: nohup, stdout to /dev/null, sterrr to logfile
	# (07/2024: this mostly works)
	nohup python -m http.server "$port" -b localhost 1>/dev/null 2>> "$SERVER_LOG" &
	# capture server PID - last backgrounded command
	server_pid_var=$!
        echo "$server_pid_var" > $SERVER_PID
    else
	echo "${PROGNAME}: Unreachable." >&2
	exit 2
    fi
    # status update to stdout
    [ "$port" != $DEFAULT_PORT ] && p_text=" -p $port " || p_text=" "
    echo "${PROGNAME}: Started the server.  Stop with: \`${PROGNAME}${p_text}stop'" 

elif [[ "$mode" == stop ]]; then
    #
    # stop the HTTP server
    #
    if [ ! -r "$SERVER_PID" ]; then
	echo "${PROGNAME}: No server is running on $host on port number $port." >&2
	echo "${PROGNAME}: Try \`${PROGNAME} status', or check $SERVER_VARDIR for another port number?" >&2
	exit 1
    fi
    echo "${PROGNAME}: Killing server named in $SERVER_PID."
    kill -TERM "$(cat "$SERVER_PID")"

elif [[ "$mode" == status ]]; then
    #
    # which servers are running on this host?
    # (we can only status-test servers on localhost)
    #
    filepat="$SERVER_VARDIR/http-$host.*.pid"
    if ls $filepat &> /dev/null; then
	echo "Looking in: $SERVER_VARDIR"
	echo "Found possible server(s):" 
	for pfile in $filepat; do
	    echo "$pfile"
	    p_info=$(stat -c 'Started by %U on %y' "$pfile")
	    x_hopo=$(echo "$pfile" | sed -e 's|.*/http-||' -e 's|\.pid||' -e 's|\.|:|')
	    # we listen on localhost, x_host is the server hostname
	    x_host=$(echo "$x_hopo" | sed -e 's|:.*||')
	    x_serv=localhost
	    x_port=$(echo "$x_hopo" | sed -e 's|.*:||')
	    # (the above filename manipulations find the host/port the server is running on)
	    echo "  " "Running on $x_host"
	    echo "  " "$p_info"
	    echo "  " "HOST:PORT = ${x_serv}:${x_port}"
	    echo "  " "View at this URL: http://${x_serv}:${x_port}/sims"
            echo "  " "Making trial request to server..."
	    if wget -q -t 1 -O /dev/null "http://${x_serv}:${x_port}/sims/"; then
		echo "  " "Server seems to respond OK to requests over HTTP."
	    else
		echo "  " "Server does not respond to requests over HTTP."
	    fi
	done
    else
	echo "No HTTP servers running on this host." 
    fi

else
    echo "${PROGNAME}: Unrecognized mode input $mode." >&2
    exit 1
fi

