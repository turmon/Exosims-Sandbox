#!/usr/bin/env bash
# 
# Top-level sh driver for certain "experiment", or multi-ensemble drm plots
#
# usage:
#   plot_experiment.sh [-d] in_experiment dest_tmpl
# where:
#   in_experiment: filename giving experiment summary
#       e.g.: Experiments/HabEx_4m_TSDD_pop100DD_revwt.json
#   dest_tmpl: the destination template for graphical output
#       e.g.: exps/HabEx_4m_TSDD_pop100DD_revwt/gfx/plot-%s.%s
# and:
#   -d is a flag which signals to use the developer MATLABPATH (~turmon) rather
#      than the regular production path.
#
# For example,
#
#  plot_experiment.sh Experiments/HabEx_4m_TSDD_pop100DD_revwt.json /tmp/revwt-plot-%s.%s
#
# See also: 
#
# Michael Turmon, JPL, Jun 2018

# exit-on-error
set -euo pipefail

progname=`basename $0`
USAGE="Usage: $progname [ -d ] in_experiment dest_tmpl"

# get options
developer_path=0
while getopts "d" o; do
    case "$o" in
	d)    developer_path=1;;
	[?])  echo "$USAGE" 1>&2
	      exit 2;;
    esac
done
shift `expr $OPTIND - 1`

# get arguments
if [ "$#" -ne 2 ]; then
    echo "$USAGE" 1>&2
    exit 2
fi

echo "${progname}: Beginning."
# ok to get the args
in_experiment="$1"
dest_tmpl="$2"

# make created files group-writable
umask 0002

# try to create the destination dir
mkdir -p "$(dirname $dest_tmpl)"
# attempt to make it group-writable (and ignore errors)
chmod g+sw "$(dirname $dest_tmpl)" 2> /dev/null || true

# set up matlab path - has to be rooted, because Matlab binary chdir's internally
rootDir=$PWD/Local/Matlab
export MATLABPATH="$rootDir/mfile:$rootDir/export_fig"

# Matlab version
matlab_binary=/usr/local/matlab-9.2/bin/matlab
# needed for R2010b (see above)
OSCHECK_ENFORCE_LIMITS=0; export OSCHECK_ENFORCE_LIMITS

# the basic matlab call
# notes:
#  the try/catch is good matlab to catch trouble
#    but, there is no way to return error status from matlab,
#    so we resort to looking for lines with ERROR
#    Thus, to play nice, the script has to start its error messages this
#    way.
tmpnam="/tmp/$progname.$$.err"
rm -f "$tmpnam"
t0=`date +%s`
$matlab_binary -nodesktop -nosplash -nodisplay -logfile "$tmpnam" <<EOF
% use outer try/catch to catch *all* exceptions
try,
  in_experiment='$in_experiment';
  dest_tmpl='$dest_tmpl';
  load_experiment_script;
catch ME,
  % ensure these errors start on a new line
  fprintf(1, '\\nERROR: %s in %s (Line %d)\\n', ME.message, ME.stack(1).name, ME.stack(1).line);
  fprintf(2, '\\nERROR: %s in %s (Line %d)\\n', ME.message, ME.stack(1).name, ME.stack(1).line);
  rethrow(ME);
end;
quit
EOF
t1=`date +%s`
tdiff=$(( $t1 - $t0 ))
echo "${progname}: Matlab call took $tdiff seconds."

# generate exit status
if grep -q "^ERROR" "$tmpnam" ; then
    # grep above matched an error
    echo "${progname}: Exiting with error." 
    echo "${progname}: Exiting with error." 1>&2
    exit 1
else
    # no error found
    echo "${progname}: Done."
    rm "$tmpnam"
    exit 0
fi

