function varargout = tom_markFindMatchesParamHelperFcn(command, varargin)


fcn = str2func(command);

[varargout{1:nargout}] = fcn(varargin{:});



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function configdefault = getDefaultConfig(nimages)

configdefault.ncorrelate = 3;
configdefault.comparelist = zeros(0, 3);
configdefault.w = 0.075;
configdefault.imreadbinning = 0;
configdefault.filter = [];
configdefault.verbose = true;
configdefault.max_disparity = 0;
configdefault.im_shift_mode = 'none';
configdefault.im_shift = tom_mark_getAvgShift(nimages);
configdefault.min_val = 0;
configdefault.mutual = true;
configdefault.matchchain = [];
configdefault.comparelist_use = true;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reads the structure config0 and parses the interesting fields
% into config. Missing fields are set by default.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function config = parseConfig(nimages, config0, configdefault)

if (~isfield(config0, 'ncorrelate') || ~isnumeric(config0.ncorrelate) || numel(config0.ncorrelate)~=1 || round(config0.ncorrelate)~=config0.ncorrelate || config0.ncorrelate<1)
    config.ncorrelate = configdefault.ncorrelate;
else
    config.ncorrelate = config0.ncorrelate;
end;

if (~isfield(config0, 'comparelist_use') || ~(isnumeric(config0.comparelist_use)||islogical(isnumeric(config0.comparelist_use))) || numel(config0.comparelist_use)~=1)
    config.comparelist_use = configdefault.comparelist_use;
else
    config.comparelist_use = config0.comparelist_use;
end;

if (~isfield(config0, 'comparelist') || ~isnumeric(config0.comparelist) || isempty(config0.comparelist) || size(config0.comparelist,2)~=3)
    config.comparelist = configdefault.comparelist;
else
    config.comparelist = config0.comparelist;
end;

if (~isfield(config0, 'w') || ~isnumeric(config0.w) || numel(config0.w)~=1 || config0.w<=0 || (config0.w>1&&round(config0.w)~=config0.w))
    config.w = configdefault.w;
else
    config.w = config0.w;
end;

if (~isfield(config0, 'imreadbinning') || ~isnumeric(config0.imreadbinning) || numel(config0.imreadbinning)~=1 || round(config0.imreadbinning)~=config0.imreadbinning || config0.imreadbinning<0)
    config.imreadbinning = configdefault.imreadbinning;
else
    config.imreadbinning = config0.imreadbinning;
end;

if (~isfield(config0, 'filter'))
    config.filter = configdefault.filter;
else
    config.filter = config0.filter;
end;

if (~isfield(config0, 'verbose') || ~(islogical(config0.mutual)||isnumeric(config0.mutual)) || numel(config0.verbose)~=1)
    config.verbose = configdefault.verbose;
else
    config.verbose = config0.verbose;
end;


if (~isfield(config0, 'max_disparity') || ~isnumeric(config0.max_disparity) || numel(config0.max_disparity)~=1 || config0.max_disparity<0)
    config.max_disparity = configdefault.max_disparity;
else
    config.max_disparity = config0.max_disparity;
end;



if (~isfield(config0, 'im_shift_mode') || ~ischar(config0.im_shift_mode))
    config.im_shift_mode = configdefault.im_shift_mode;
else
    if (any(strcmp(config0.im_shift_mode, {'none', 'zero', 'set'})))
        config.im_shift_mode = config0.im_shift_mode;
    else
        config.im_shift_mode = configdefault.im_shift_mode;
    end;
end;

if (~isfield(config0, 'im_shift') || ~isnumeric(config0.im_shift) || ndims(config0.im_shift)~=3 || any(size(configdefault.im_shift)~=[nimages*[1 1], 2]))
    config.im_shift = configdefault.im_shift;
else
    config.im_shift = config0.im_shift;
end;

if (~isfield(config0, 'min_val') || ~isnumeric(config0.min_val) || numel(config0.min_val)~=1)
    config.min_val = configdefault.min_val;
else
    config.min_val = config0.min_val;
end;

if (~isfield(config0, 'mutual') || ~(islogical(config0.mutual)||isnumeric(config0.mutual)) || numel(config0.mutual)~=1)
    config.mutual = configdefault.mutual;
else
    config.mutual = config0.mutual;
end;


if (~isfield(config0, 'matchchain'))
    config.matchchain = configdefault.matchchain;
else
    config.matchchain = config0.matchchain;
end;

