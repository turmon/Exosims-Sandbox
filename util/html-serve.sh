#!/bin/bash
#
# html-serve.sh: start/stop HTTP server
#
# This starts or stops the HTTP server that serves up the Exosims
# table of contents and ensemble-summary files.
# We run the server on a non-standard port - 8088 by default.  So,
# to see the table of contents, navigate to:
#  http://aftac1.jpl.nasa.gov:8088/sims/
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
#  -p PORT => gives the HTTP port number.  Default is 8088.
#  -h      => shows this help text.
#
# Implementation note: Status files are kept in files called:
#   Local/www-service/http.PORT.*
# where PORT is typically 8088.
# There is a server log, a process ID of the server, and a fifo that
# is used to signal the running server to exit.  If the basic "stop" 
# MODE of this command does not work, you can kill the server by the
# process ID (pid) in that directory.
# 
# turmon sep 2018
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
# default port
DEFAULT_PORT=8088


port=$DEFAULT_PORT
while getopts "hp:" opt; do
    case $opt in
	p)
	    # port number
	    port="$OPTARG"
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

# specific filenames we need
SERVER_FIFO=$SERVER_DIR/http.$port.fifo
SERVER_LOG=$SERVER_DIR/http.$port.log
SERVER_PID=$SERVER_DIR/http.$port.pid


# set up the server state dir
function setup {
    # (assume it is already created and has right mode)
    # mkdir -m 775 -p "$SERVER_DIR"
    # chgrp $SERVER_GROUP "$SERVER_DIR"
    # re-create fifo
    rm -f "$SERVER_FIFO"
    mkfifo -m 664 "$SERVER_FIFO"
    chgrp $SERVER_GROUP "$SERVER_FIFO"
    # truncate and re-create server log
    rm -f "$SERVER_LOG" && touch "$SERVER_LOG"
    echo "${PROGNAME}: HTTP file setup OK." >> "$SERVER_LOG"
    chgrp $SERVER_GROUP "$SERVER_LOG"
    chmod 664 "$SERVER_LOG"
    # truncate server PID file
    rm -f "$SERVER_PID" && touch "$SERVER_PID"
    chgrp $SERVER_GROUP "$SERVER_PID"
    chmod 664 "$SERVER_PID"
}


if [ "$mode" == start ]; then
    #
    # start the HTTP server
    #
    # Try not to have two going on same port!
    if [ -r $SERVER_FIFO ]; then
	echo "${PROGNAME}: Server may be running on $port already." >&2
	echo "${PROGNAME}: Attempting to shut down cleanly and restart." >&2
	echo > $SERVER_FIFO &
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
	httpd -f "$CURR_DIR/$SERVER_DIR/httpd.conf" \
	      -c "PidFile $CURR_DIR/$SERVER_PID" \
	      -c "ErrorLog $CURR_DIR/$SERVER_LOG" \
	      -c "Listen $port"
	# capture server PID -- httpd forks, so the pid changes
	server_pid_var=$(sleep 0.1 && cat "$CURR_DIR/$SERVER_PID")
    else
	# python SimpleHTTPServer: single-threaded, does not support returning
	# byte ranges, which makes .mp4s not render in Safari.
	# mostly-daemonized: nohup, stdout to /dev/null, sterrr to logfile
	nohup python -m SimpleHTTPServer "$port" 1>/dev/null 2>>$SERVER_LOG &
	# capture server PID - last backgrounded command
	server_pid_var=$!
    fi
    # last status update to stdout
    [ $port != $DEFAULT_PORT ] && p_text="-p $port" || p_text=""
    echo "${PROGNAME}: Started the server.  Stop with: \`${PROGNAME} $p_text stop'" 
    # this is just informational, the correct way to kill is with the FIFO
    echo $server_pid_var >> $SERVER_PID
    # mostly-daemonize ourselves
    trap "" SIGHUP # ignore exit of parent shell
    exec >/dev/null
    exec 2>/dev/null
    exec 0</dev/null
    # wait for a message on the FIFO, kill the server when received, clean up.
    # all done in background, so we can terminate this command.
    { read < $SERVER_FIFO; 
      echo "${PROGNAME}: Killing the HTTP server at $port" >> $SERVER_LOG;
      kill "$server_pid_var";
      rm -f $SERVER_FIFO $SERVER_PID; } &
elif [ "$mode" == new ]; then
    #
    # start the HTTP server
    #
    # Try not to have two going on same port!
    if [ -r $SERVER_FIFO ]; then
	echo "${PROGNAME}: Server may be running on $port already." >&2
	echo "${PROGNAME}: Attempting to shut down cleanly and restart." >&2
	echo > $SERVER_FIFO &
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
    # apache httpd, minimal config, PID and logging configured on command line
    httpd -f "$CURR_DIR/httpd.conf" \
	  -c "PidFile $CURR_DIR/$SERVER_PID" \
	  -c "ErrorLog $CURR_DIR/$SERVER_LOG" \
	  -c "Listen $port"
    # capture server PID -- httpd forks, so the pid changes
    server_pid_var=$(sleep 0.1 && cat "$CURR_DIR/$SERVER_PID")
    # last status update to stdout
    [ $port != $DEFAULT_PORT ] && p_text="-p $port" || p_text=""
    echo "${PROGNAME}: Started the server.  Stop with: \`${PROGNAME} $p_text stop'" 
    # this is just informational, the correct way to kill is with the FIFO
    echo $server_pid_var >> $SERVER_PID
    # mostly-daemonize ourselves
    trap "" SIGHUP # ignore exit of parent shell
    exec >/dev/null
    exec 2>/dev/null
    exec 0</dev/null
    # wait for a message on the FIFO, kill the server when received, clean up.
    # all done in background, so we can terminate this command.
    { read < $SERVER_FIFO; 
      kill "$server_pid_var";
      echo "${PROGNAME}: Killing the HTTP server at $port" >> $SERVER_LOG;
      rm -f $SERVER_FIFO $SERVER_PID; } &
elif [ "$mode" == stop ]; then
    #
    # stop the HTTP server
    #
    if [ ! -r $SERVER_FIFO ]; then
	echo "${PROGNAME}: No server is running on $port." >&2
	echo "${PROGNAME}: Try \`${PROGNAME} status', or check $SERVER_DIR for another port number?" >&2
	exit 1
    fi
    echo "${PROGNAME}: Killing server at $SERVER_FIFO." 
    echo > $SERVER_FIFO
elif [ "$mode" == init ]; then
    #
    # ensure that files can be set up -- not necessary
    #
    setup
    # don't want the files hanging round...
    # but this means they could be created OK
    rm -f $SERVER_LOG $SERVER_FIFO $SERVER_PID
    echo "${PROGNAME}: Initialized OK (see $SERVER_DIR)." 
elif [ "$mode" == status ]; then
    #
    # which servers are running?
    #
    filepat="$SERVER_DIR/*.pid"
    if ls $filepat &> /dev/null; then
	echo "Found server(s):" 
	for pfile in $filepat; do
	    p_info=$(stat -c 'started by %U on %y' $pfile)
	    p_num=$(echo $pfile | sed -e 's|.*/http\.||' -e 's|\..*||')
	    echo "  " "PORT=$p_num" " -- " "$p_info"
	done
	echo "Should be able to view at http://$HOST:PORT/sims"
    else
	echo "No HTTP servers running." 
    fi
fi

