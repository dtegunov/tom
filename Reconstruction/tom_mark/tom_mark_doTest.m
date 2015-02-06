function tom_mark_doTest(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

optionlist = { ...
        'refine',           @doTest_refineMarkerSet; ...
        'getMS',            @doTest_getMarkerSetFromScrach; ...
    };

option = optionlist(1,1);
%option = '';
if (nargin >= 1 && isstr(varargin{1}))
    for (i=1:size(optionlist,1))
        if (strcmp(optionlist{i,1}, varargin{1}))
            option = varargin{1};
            varargin = varargin(2:end);
            break;
        end;
    end;
end;

function_executed = false;
for (i=1:size(optionlist,1))
    if (strcmp(option, optionlist{i,1}))
        fhandle = optionlist{i,2};
        function_executed = true;
        fhandle(varargin{:});
        break;
    end;
end;

if (~function_executed)
    warning('doTest did nothing.');
end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function doTest_getMarkerSetFromScrach(varargin)

% Relative size of the matchtemplate in relation to the original imagesize.
match1_match_w_relative = 0.125;


% Number of points in the first run.
match1_initMarker_nPoints = 10;


% The minimum imagesize handled to tom_mark_initMarker
match1_initMarker_imageSize = 256;
match1_initMarker_distancefactor = 1/2;
match1_match_imageSize = 256;
match1_match_ncorrelate = 4;


% The filenames of the images.
%filename_base = '/fs/bmsan/apps/tom_dev/Reconstruction/tom_mark/test_data/cryosection/CRY/cryosection';
filename_base = '/fs/bmsan/apps/tom_dev/Reconstruction/tom_mark/test_data/dualaxis/DUAL/pyMOD/pyMOD';



% Compute the imagefilenames
filenames = {};
i = 1;
emimage_Header = {};
while (true);
    s = [filename_base '_' num2str(i) '.em'];
    if (exist(s, 'file') && ~exist(s, 'dir'))
        TMP = tom_reademheader(s);
        emimage_Header{i} = TMP.Header;
        filenames{i} = s;
        i = i+1;
    else
        break;
    end;
end;
if (isempty(filenames))
    error('Files not found');
end;


msize = [2 length(filenames) 0];

% Read the Headers of all images.
emimage_Size = emimage_Header{1}.Size([1 2])';
emimage_Tiltangle = nan(1, msize(2));
for (i=1:msize(2))
    TMP = emimage_Header{i}.Size;
    if (TMP(3) ~= 1 || TMP(1) ~= emimage_Size(1) || TMP(2) ~= emimage_Size(2))
        error(['The ' num2str(i) 'th image ("' filenames{i} '") is either not an (2 dimensional) image or the imagesize in the sequence differs.']);
    end;
    emimage_Tiltangle(i) = emimage_Header{i}.Tiltangle;
end;



searchMarkerSize = ceil(match1_match_w_relative * min(emimage_Size));
tom_mark_initMarker_binning = max(0,fix(log2(min(emimage_Size) / match1_initMarker_imageSize)));
markerset_init = nan(2, msize(2), match1_initMarker_nPoints);

% match1_initMarker_imageIdx = fix(mean(1:msize(2)) + [-1 1]*std(1:msize(2)/3))
% match1_initMarker_imageIdx = 1:max(1,fix(msize(2)/(nPointsTot/match1_initMarker_nPoints))):msize(2)
% match1_initMarker_imageIdx = subsref(find(abs(emimage_Tiltangle) == min(abs(emimage_Tiltangle))), struct('type','()','subs',{{1}}));
match1_initMarker_imageIdx = fix(msize(2)/2);



for (i=match1_initMarker_imageIdx)

    im = tom_emreadc(filenames{i}, 'binning', tom_mark_initMarker_binning);
    
    %tic 
    ms = tom_mark_initMarker(im.Value', floor(match1_initMarker_nPoints/length(match1_initMarker_imageIdx)), ...
                             (min(emimage_Size)/sqrt(floor(match1_initMarker_nPoints/length(match1_initMarker_imageIdx)))) / 2^tom_mark_initMarker_binning * match1_initMarker_distancefactor, ...
                             searchMarkerSize*2^(-1-tom_mark_initMarker_binning)) * 2^tom_mark_initMarker_binning;

    %toc
    markerset_init([1 2], i, (msize(3)+1):msize(3)+size(ms,3)) = ms;
    msize(3) = msize(3) + size(ms, 3);
    
    if (false)
        %imagesc([1 emimage_Size(1)], [1 emimage_Size(2)], im.Value');
        im = tom_emreadc(filenames{i});
        imagesc([1 emimage_Size(1)], [1 emimage_Size(2)], im.Value');
        colormap gray;
        hold on;
        axis equal off;
        plot(squeeze(ms(1,1,:)), squeeze(ms(2,1,:)), 'gx', 'MarkerSize', 8, 'LineWidth', 2);
        disp(['initMarker read ' num2str(size(ms,3)) ' coordinates from image #' num2str(i) ]);
    end;
    
end;

%return;

if (msize(3) < size(markerset_init, 3))
    markerset_init = markerset_init([1 2], 1:msize(2), 1:msize(3));
end;



vbase = struct();
vbase.filename_base = filename_base;
vbase.match1_initMarker_imageIdx = match1_initMarker_imageIdx;
vbase.match1_initMarker_nPoints = match1_initMarker_nPoints;
vbase.match1_initMarker_imageSize = match1_initMarker_imageSize;
vbase.match1_match_w_relative = match1_match_w_relative;
vbase.match1_match_imageSize = match1_match_imageSize;
vbase.match1_match_ncorrelate = match1_match_ncorrelate;
vbase.filenames = filenames;
vbase.emimage_Header = emimage_Header;
vbase.emimage_Tiltangle = emimage_Tiltangle;
vbase.emimage_Size = emimage_Size;
vbase.markerset_init = markerset_init;
vbase.msize = msize;

match1_match_matchparam.w = match1_match_w_relative;
match1_match_matchparam.max_disparity = 0;
match1_match_matchparam.im_shift = 'default';
match1_match_matchparam.min_val = 0.;
match1_match_matchparam.mutual = true;
vbase.match1_match_matchparam = match1_match_matchparam;

assignin('base', 'doTest', vbase);


[markerset_new, match_chain, nocompare] = ...
    tom_mark_findmatches(filenames, markerset_init, match1_match_ncorrelate, [], ...
                         match1_match_matchparam, max(0, fix(log2(min(emimage_Size) / match1_match_imageSize))), []);
 

vbase.markerset_new = markerset_new;                     
vbase.match_chain = match_chain;                     
vbase.nocompare = nocompare;                     
assignin('base', 'doTest', vbase);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function doTest_refineMarkerSet(varargin)

if (nargin == 0)
    doTest = load('doTest_pyMod');
    doTest = doTest.doTest;
end;
if (nargin >= 1)
    doTest = varargin{1};
end;


match1_match_matchparam.w = doTest.match1_match_w_relative;
match1_match_matchparam.max_disparity = 100;
match1_match_matchparam.min_val = 0.;
match1_match_matchparam.mutual = true;

markerset_new = doTest.markerset_new(:,:,1:5);

global gmarkerset;
global GLOBALFILTER;
gmarkerset = markerset_new;

tic
match1_match_matchparam.w = 256;
match1_match_matchparam.max_disparity = 25;
match1_match_matchparam.min_val = 0.;
match1_match_matchparam.mutual = false;


imreadbinning = 1;


transf(1).Type = 'bandpass';
transf(end).Apply = 1;
transf(end).Value.times = 1;
transf(end).Value.low = 2;
transf(end).Value.high = 2048;
transf(end).Value.smooth = 0;
transf(end).Value.space = 'real';
transf(end).Value.mathod = 'quadr';
transf(end).Value.radius = 0;




[markerset_new] = tom_mark_refineMarkerSet(doTest.filenames, markerset_new, 2, match1_match_matchparam, imreadbinning, transf);
gmarkerset = cat(3, gmarkerset, markerset_new);
toc

return;
tic
match1_match_matchparam.w = 0.1;
match1_match_matchparam.max_disparity = 50;
match1_match_matchparam.min_val = 0.;
match1_match_matchparam.mutual = false;
[markerset_new] = tom_mark_refineMarkerSet(doTest.filenames, markerset_new, 15, match1_match_matchparam, 1, GLOBALFILTER);
gmarkerset = cat(3, gmarkerset, markerset_new);
toc

tic
match1_match_matchparam.w = 0.08;
match1_match_matchparam.max_disparity = 10;
match1_match_matchparam.min_val = 0.;
match1_match_matchparam.mutual = false;
[markerset_new] = tom_mark_refineMarkerSet(doTest.filenames, markerset_new, 15, match1_match_matchparam, 0, GLOBALFILTER);
gmarkerset = cat(3, gmarkerset, markerset_new);
toc


return;

