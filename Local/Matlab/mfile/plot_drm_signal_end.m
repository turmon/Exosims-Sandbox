function plot_drm_signal_end(dest_tmpl, mode, t_info)
%plot_drm_signal_end	write a file to signal end of processing
% 
% plot_drm_signal_end(dest_tmpl, mode, t_info)
% * Write a sentinel file whose main purpose is to signal that the
% processing ended OK.
% 
% Inputs:
%   string dest_tmpl
%   string mode
%   table t_info
% 
% Outputs:
%   (to disk)
% 
% See Also:  

% 
% Error checking
% 
% if all(nargin  ~= [3]), error ('Bad input arg number'); end
% if all(nargout ~= [0 1]), error ('Bad output arg number'); end  
% 


%
% Signal success
%
fn_output = sprintf(dest_tmpl, 'info', 'txt');
fp = fopen(fn_output, 'w');
fprintf(fp, 'Graphics written by %s\n', getenv('USER'));
fprintf(fp, 'From a run that was reduced on %s\n', t_info{1,'runtime'}{1});
fclose(fp);

return
