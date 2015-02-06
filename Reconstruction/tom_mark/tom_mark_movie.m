function tom_mark_movie(filenames, markerset, irefmarker, w, bin, delay)


% Check input parameters and set defaults...
msize = size(markerset);
if (length(msize) < 3)
    msize(3) = 1;
end;
if (length(filenames) ~= msize(2))
    error('Number of images and number of projections in markerset differ.');    
end;
if (msize(1) > 2)
    markerset = markerset([2 3], 1:msize(2), 1:msize(3));
end;
if (~exist('irefmarker','var') || numel(irefmarker)~=1)
    irefmarker = 1;
end;
if (~exist('delay', 'var') || numel(delay) ~= 1 || delay < 0)
    delay = 0;
end;
if (~exist('bin', 'var') || numel(bin) ~= 1 || bin < 0)
    bin = 0;
end;
bin = fix(bin);
if (~exist('w', 'var') || (numel(w)~=1 && numel(w)~=2) || any(~(w>0)))
    w = [];
end;
if (~isempty(w))
    if (numel(w) == 1)
        w = [w w];
    end;
end;
w_half = fix(w/2);
w = 2*w_half + 1;



imsize = [];
for (i=1:msize(2))
    try 
        im = [];
        im = tom_reademheader(filenames{i});
    catch
        continue;
    end;
    if (isempty(im) || ~isfield(im, 'Header') || isempty(im.Header))
        continue;
    end;
    if (im.Header.Size(3) ~= 1)
        error('File is not an image');
    end;
    if (isempty(imsize))
        imsize = im.Header.Size([1 2]);
    else
        if (any(imsize ~= im.Header.Size([1 2])))
            error('dimension of images differ');
        end;
    end;
end;
if (isempty(imsize))
    error('No image found');
end;


h = figure;
haxes = gca();
set(haxes, 'YDir', 'reverse');
colormap gray;
hold on;


    
lastidx = 0;
for (i=1:msize(2))

    loaded = show_it(filenames{i}, markerset, i, irefmarker, bin, imsize, w_half);
    
    if (loaded)
        lastidx = i;
        pause(delay);
        drawnow;
    end;
end;

close(h);

% 
% if (isempty(findobj('Parent', h, 'String', 'previous')))
%     uicontrol('Parent', h, 'Style', 'pushbutton', 'Position', [  1   1 100  20], 'Callback', @callback', 'String', 'previous');
% end;
% if (isempty(findobj('Parent', h, 'String', 'next')))
%     uicontrol('Parent', h, 'Style', 'pushbutton', 'Position', [101   1 100  20], 'Callback', @callback', 'String', 'next');
% end;
% if (isempty(findobj('Parent', h, 'String', 'again')))
%     uicontrol('Parent', h, 'Style', 'pushbutton', 'Position', [  1  21 100  20], 'Callback', @callback', 'String', 'again');
% end;
% 
% 
% ss.filenames = filenames;
% ss.markerset = markerset;
% ss.lastidx = lastidx;
% ss.irefmarker = irefmarker;
% ss.bin = bin;
% ss.imsize = imsize;
% ss.w_half = w_half;
% ss.msize = msize;
% ss.delay = delay;
% 
% guidata(h, ss);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function callback(hObject, eventdata, handles)
% 
% ss = guidata(hObject);
% 
% button = get(hObject, 'String');
% switch (button)
%     case {'next','previous'}
%         if (strcmp(button, 'next'))
%             ss.lastidx = ss.lastidx + 1;
%         else
%             ss.lastidx = ss.lastidx - 1;
%         end;
%         ss.lastidx = min(ss.msize(2), max(1, ss.lastidx));
%         
%         show_it(ss.filenames{ss.lastidx}, ss.markerset, ss.lastidx, ss.irefmarker, ss.bin, ss.imsize, ss.w_half);
%         
%     otherwise
%         ss.lastidx = 0;
%         for (i=1:ss.msize(2))
% 
%             loaded = show_it(ss.filenames{i}, ss.markerset, i, ss.irefmarker, ss.bin, ss.imsize, ss.w_half);
% 
%             if (loaded)
%                 ss.lastidx = i;
%                 pause(ss.delay);
%                 drawnow;
%             end;
%         end;
%         
% end;
% 
% guidata(hObject, ss);
%         
%     
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loaded = show_it(filename, markerset, irefimage, irefmarker, bin, imsize, w_half);

loaded = false;

if (~all(markerset([1 2], irefimage, irefmarker) > 0))
    return;
end;

if (isempty(w_half))
    wrect = [[1, 1]', imsize];
    rect = wrect;
else
    wrect = [markerset([1 2], irefimage, irefmarker) - w_half', markerset([1 2], irefimage, irefmarker) + w_half'];
    rect = [ max(1, wrect(:,1)), min(imsize, wrect(:,2)) ];
end;

try 
    im = tom_emreadc(filename, 'binning', bin, 'subregion', [rect(:,1)', 1], [rect(:,2)'-rect(:,1)', 0]);
catch
    return;
end;
if (isempty(im.Value))
    return;
end;

cla;


imagesc(rect(1,:), rect(2,:), im.Value');
axis equal %off


xlim([wrect(1,1), wrect(1,2)]');
ylim([wrect(2,1), wrect(2,2)]');

idx = squeeze( and( and(markerset(1, irefimage, :) >= rect(1,1), markerset(1, irefimage, :) <= rect(1,2)), ...     
               and(markerset(2, irefimage, :) >= rect(2,1), markerset(2, irefimage, :) <= rect(2,2)) ) );
idx(irefmarker) = false;

plot(squeeze(markerset(1, irefimage, idx)), squeeze(markerset(2, irefimage, idx)), 'g+', 'MarkerSize', 8, 'LineWidth', 2)               
plot(markerset(1, irefimage, irefmarker), markerset(2, irefimage, irefmarker), 'rx', 'MarkerSize', 15, 'LineWidth', 2);



loaded = true;






