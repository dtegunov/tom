function [ind,t0] = tom_crossing(S,t,level,par)
% CROSSING find the crossings of a given level of a signal
%   ind = CROSSING(S) returns an index vector ind, the signal
%   S crosses zero at ind or at between ind and ind+1
%   [ind,t0] = CROSSING(S,t) additionally returns a time
%   vector t0 of the zero crossings of the signal S. The crossing
%   times are linearly interpolated between the given times t
%   [ind,t0] = CROSSING(S,t,level) returns the crossings of the
%   given level instead of the zero crossings
%   ind = CROSSING(S,[],level) as above but without time interpolation
%   [ind,t] = CROSSING(S,t,level,par) allows additional parameters
%       par.interpolation = {'none'|'linear'}
%

% Steffen Brueckner, 2002-09-25
% Copyright (c) Steffen Brueckner, 2002
% brueckner@sbrs.net

% check the number of input arguments
error(nargchk(1,4,nargin));

if nargin < 2 || isempty(t)
    % if not given: time vector == index vector
    t = 1:length(S);
elseif length(t) ~= length(S)
    error('t and S must be of identical length!');
end
if nargin < 3
    level = 0;
end
if nargin < 4
    par.interpolation = 'linear';
end

S   = S - level;
S1  = S(1:end-1) .* S(2:end);
ind = find(S1 <= 0);

% get times of crossings
t0 = double(t(ind));

if strcmp(par.interpolation,'linear')
    % linear interpolation of crossing
    for ii=1:length(t0)
        if abs(S(ind(ii))) > eps
            % interpolate only when data point is not already zero
            t0(ii) = t0(ii) - S(ind(ii)) .* (t(ind(ii)+1) - t(ind(ii))) ./ (S(ind(ii)+1) - S(ind(ii)));
        end
    end
end

end