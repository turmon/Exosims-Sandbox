function titl = plot_make_title(t_info)
%plot_make_title(t_info)	make title string from t_info structure
% 
% titl = plot_make_title(t_info)
% * Format a textual experiment title from the overall info structure.
%
% Inputs:
%   table t_info
% 
% Outputs:
%   string titl
% 
% See Also:  

% 
% Error checking
% 
if all(nargin  ~= [1]), error ('Bad input arg number'); end

% allow passing an empty metadata table as []
if ismatrix(t_info) && isempty(t_info),
    t_info = cell2table(cell(0));
end;

% format the title
if any(strcmp('experiment', fieldnames(t_info))),
    exp_name = t_info{1,'experiment'}{1};
    if length(exp_name) > 0 && exp_name(1) == ' ',
        % auto-generated title, derived from sim dir name
        titl = sprintf('%s, Ensemble Size %d', ...
                         exp_name(2:end), t_info{1,'ensemble_size'});
    else,
        % manually-set title, from EnsembleName.txt file -- will not
        % begin with space
        titl = exp_name;
    end;
else,
    % should never happen
    titl = '';
end

return
end
