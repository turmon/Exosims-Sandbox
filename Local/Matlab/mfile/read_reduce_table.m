function t = read_reduce_table(in_tmpl, name)
%read_reduce_table(in_tmpl, name)	read one reduce-*.csv table, gzipped or not
% 
% t = read_reduce_table(in_tmpl, name)
% * Fill in_tmpl (a filename template with two %s's) with NAME and the file
% extension, and read the resulting table.
% * The larger reduction outputs -- the earth-char list, the planet
% population -- are written by reduce_drms.py as .csv.gz, which readtable
% cannot open, so such a file is unpacked to a temporary directory first.
% * A plain .csv wins when both are present, matching the python readers
% (see common_style.resolve_csv_path).  If neither exists we read the .csv
% anyway, so the caller sees readtable's usual error, naming the file it
% would have expected.
%
% Inputs:
%   string in_tmpl -- template with two %s's (name, extension)
%   string name -- table name, e.g. 'earth-char-list'
% 
% Outputs:
%   table t
% 
% See Also:  plot_drms_script

% 
% Error checking
% 
if all(nargin  ~= [2]), error ('Bad input arg number'); end

fn = sprintf(in_tmpl, name, 'csv');
fn_gz = sprintf(in_tmpl, name, 'csv.gz');

if exist(fn, 'file') || ~exist(fn_gz, 'file'),
    t = readtable(fn);
    return;
end;

% unpack the .gz to a scratch directory, read it, and clean up
tmp_dir = tempname();
mkdir(tmp_dir);
cleanup = onCleanup(@() rmdir(tmp_dir, 's'));
fn_out = gunzip(fn_gz, tmp_dir);
t = readtable(fn_out{1});
