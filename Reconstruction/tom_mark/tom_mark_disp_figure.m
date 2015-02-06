function h = tom_mark_disp_figure(name, varargin)
% 


ftag = ['tom_mark_disp_figure__' name];
h = findobj('Type', 'figure', 'Tag', ftag);

if (isempty(h))
    h = figure();
    set(h, 'Tag', ftag);
else
    if (gca ~= h)
        figure(h);
    end;
end;

if (~isempty(varargin))
    set(h, varargin{:});
end;

