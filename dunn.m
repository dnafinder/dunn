function [stats, results] = dunn(varargin)
%DUNN Stepdown Dunn procedure for nonparametric multiple comparisons.
%
%   DUNN(X1, X2, ..., XK)
%   DUNN(X1, X2, ..., XK, 'Ctrl', CTRL, 'Display', DISPLAY)
%
%   Performs Dunn's stepdown procedure for multiple comparisons between
%   K independent groups using rank-based statistics.
%
%   Each Xi is the data vector for the i-th group (row or column, with
%   possibly different lengths).
%
%   Name–Value pair arguments:
%     'Ctrl'    - logical-like flag (true/false, 0/1, 'on'/'off',
%                 'true'/'false', 'yes'/'no')
%                 true  -> the first group (X1) is a control group
%                 false -> no control group (default)
%
%     'Display' - logical-like flag controlling command-window output:
%                 true  -> print tables and messages (default)
%                 false -> run silently, only return outputs
%
%   Outputs:
%     stats   - structure with summary statistics for groups and procedure
%     results - table with pairwise comparison results
%
%   If no output arguments are requested, the function behaves in a
%   MATLAB-style way: it prints to the command window (unless Display is
%   set to false) and does not return stats or results.
%
%   Example:
%
%   x1 = [7.68 7.69 7.70 7.70 7.72 7.73 7.73 7.76];
%   x2 = [7.71 7.73 7.74 7.74 7.78 7.78 7.80 7.81];
%   x3 = [7.74 7.75 7.77 7.78 7.80 7.81 7.84];
%   x4 = [7.71 7.71 7.74 7.79 7.81 7.85 7.87 7.91];
%
%   % No control group, with printed output:
%   dunn(x1, x2, x3, x4)
%
%   % First sample is the control group, no printed output:
%   [stats, results] = dunn(x1, x2, x3, x4, 'Ctrl', true, 'Display', false);
%
%   Created by Giuseppe Cardillo
%   giuseppe.cardillo.75@gmail.com
%
%   To cite this file, this would be an appropriate format:
%   Cardillo G. (2006). Dunn's Test: a procedure for multiple
%   nonparametric comparisons.
%   Available on GitHub: https://github.com/dnafinder/dunn

% -------------------------------------------------------------------------
% Input parsing: data arguments and Name–Value pairs
% -------------------------------------------------------------------------

narginchk(2, Inf);

% Separate data arguments (groups) from Name–Value pairs
isNameArg = @(c) ischar(c) || (isstring(c) && isscalar(c));
firstNV   = find(cellfun(isNameArg, varargin), 1, 'first');

if isempty(firstNV)
    dataArgs = varargin;
    nvArgs   = {};
else
    dataArgs = varargin(1:firstNV-1);
    nvArgs   = varargin(firstNV:end);
end

k = numel(dataArgs);
assert(k >= 2, 'DUNN requires at least two groups.');

% Parse Name–Value pairs
p = inputParser;
p.FunctionName = mfilename;
addParameter(p, 'Ctrl', false, @validateLogicalLike);
addParameter(p, 'Display', true, @validateLogicalLike);
parse(p, nvArgs{:});

ctrl        = logical(normalizeLogicalLike(p.Results.Ctrl));
displayFlag = logical(normalizeLogicalLike(p.Results.Display));

% Collect data and build x and g
xCell = cell(1,k);
N     = zeros(k,1); % sample size per group

for ii = 1:k
    Xi = dataArgs{ii};
    validateattributes(Xi, {'numeric'}, ...
        {'vector','real','finite','nonnan','nonempty'}, ...
        mfilename, sprintf('X%d',ii), ii);
    Xi        = Xi(:);          % column vector
    xCell{ii} = Xi;
    N(ii)     = numel(Xi);
end

% Concatenate data (row vector)
x = vertcat(xCell{:}).';

% Group labels 1..k, ordered like the inputs
groups = (1:k).';
gCell  = arrayfun(@(ii) ii*ones(N(ii),1), 1:k, 'UniformOutput', false);
g      = vertcat(gCell{:}).';

% -------------------------------------------------------------------------
% Main computation
% -------------------------------------------------------------------------

tot = sum(N);                 % total number of observations
[r, tieadj] = tiedrank(x);    % ranks and tie adjustment
f  = (tot*(tot+1)/12) - (tieadj/(6*(tot-1))); % variance factor with ties

% Sum and mean of ranks for each group
sumRanks  = zeros(k,1);
for ii = 1:k
    sumRanks(ii) = sum(r(g==ii));
end
meanRanks = sumRanks ./ N;

% Matrix for comparisons: col 1 = mean ranks, col 2 = group index
Mr = [meanRanks, groups];

% Group summary table
Tgroups = table(groups, N, sumRanks, meanRanks, ...
    'VariableNames', {'Group','N','Sum_of_ranks','Mean_rank'});

alphaGlobal = 0.05;
count       = 0;

% Stepdown procedure
if ~ctrl
    % No control group: all pairwise comparisons with stepdown
    kstar    = 0.5*k*(k-1);
    alphaEff = 1 - (1-alphaGlobal)^(1/kstar);
    vcrit    = -sqrt(2).*erfcinv(2-alphaEff);
    
    % Sort mean ranks in descending order
    Mr = sortrows(Mr, -1);
    
    pb = cell(kstar,4); % preallocate results
    
    for I = 1:k-1
        comp = true; % comparison checker
        for J = k:-1:(I+1)
            count = count + 1;
            if comp % comparison is necessary
                Qvalue();
            else % comparison is not necessary
                pb(count,:) = { ...
                    sprintf('%d-%d', Mr(I,2), Mr(J,2)), ...
                    'No comparison made', ...
                    vcrit, ...
                    'Accept H0'};
            end
        end
    end
else
    % With control group: first group is the control
    alphaEff = 1 - (1-alphaGlobal)^(1/(k-1));
    vcrit    = -sqrt(2).*erfcinv(2-alphaEff);
    
    pb = cell(k-1,4);
    I  = 1; % control group index in Mr (first row, unsorted)
    for J = 2:k
        count = count + 1;
        Qvalue();
    end
end

% Pairwise comparison results table
results = cell2table(pb, 'VariableNames', {'Comparison','Q_value','Crit_Q','Comment'});

% -------------------------------------------------------------------------
% Build outputs
% -------------------------------------------------------------------------

stats = struct();
stats.groups      = groups;
stats.N           = N;
stats.sumRanks    = sumRanks;
stats.meanRanks   = meanRanks;
stats.tiesFactor  = 2*tieadj;
stats.varFactor   = f;
stats.totalN      = tot;
stats.alphaGlobal = alphaGlobal;
stats.alphaEff    = alphaEff;
stats.critQ       = vcrit;
stats.k           = k;
stats.ctrl        = ctrl;

% Display (if requested)
if displayFlag
    disp('STEPDOWN DUNN TEST FOR NON PARAMETRIC MULTIPLE COMPARISONS');
    disp(' ');
    disp(Tgroups);
    fprintf('Ties factor: %.4g\n', 2*tieadj);
    disp(' ');
    disp(results);
end

% If no output is requested, clear outputs (MATLAB-style)
if nargout == 0
    clear stats results
end

% -------------------------------------------------------------------------
% Nested function for Q calculation
% -------------------------------------------------------------------------
    function Qvalue
        % Numerator: absolute difference of mean ranks between groups I and J
        num = abs(diff(Mr([I J],1)));
        
        % Denominator with ties correction: depends on sample sizes
        denom = sqrt(f * sum(1./N(Mr([I J],2))));
        
        % Q value
        Q = num/denom;
        
        % Store comparison
        pb(count,1:3) = { ...
            sprintf('%d-%d', Mr(I,2), Mr(J,2)), ...
            Q, ...
            vcrit};
        
        % Decision rule
        if Q > vcrit
            pb(count,4) = {'Reject H0'};
            comp = true;  %#ok<NASGU> % more comparisons are required
        else
            pb(count,4) = {'Fail to reject H0'};
            comp = false; %#ok<NASGU> % no more intermediate comparisons are required
        end
    end

end

% -------------------------------------------------------------------------
% Local helper functions
% -------------------------------------------------------------------------

function tf = validateLogicalLike(x)
%VALIDATELOGICALLIKE Helper for inputParser: check logical-like values.
    try
        normalizeLogicalLike(x);
        tf = true;
    catch
        tf = false;
    end
end

function y = normalizeLogicalLike(x)
%NORMALIZELOGICALLIKE Convert various logical-like inputs to true/false.
    if islogical(x)
        y = x;
    elseif isnumeric(x) && isscalar(x)
        y = (x ~= 0);
    elseif ischar(x) || (isstring(x) && isscalar(x))
        s = lower(char(x));
        if any(strcmp(s, {'true','on','yes'}))
            y = true;
        elseif any(strcmp(s, {'false','off','no'}))
            y = false;
        else
            error('Invalid logical-like value: %s', s);
        end
    else
        error('Invalid type for logical-like option.');
    end
end
