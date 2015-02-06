function varargout = tom_ctftool(varargin)
%TOM_AMIRA_CREATEISOSURFACE creates ...
%
%   varargout = tom_ctftool(varargin)
%
%PARAMETERS
%
%  INPUT
%   varargin            ...
%  
%  OUTPUT
%   varargout   		...
%
%EXAMPLE
%   ... = tom_amira_createisosurface(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by AK 03/01/06
%
%   Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%   Journal of Structural Biology, 149 (2005), 227-234.
%
%   Copyright (c) 2004-2007
%   TOM toolbox for Electron Tomography
%   Max-Planck-Institute of Biochemistry
%   Dept. Molecular Structural Biology
%   82152 Martinsried, Germany
%   http://www.biochem.mpg.de/tom


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_ctftool_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_ctftool_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Opening Function                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tom_ctftool_OpeningFcn(hObject, eventdata, handles, varargin)

global storage_ctf;

storage_ctf.displaycenter = 0;
storage_ctf.logscale = 1;
storage_ctf.semicircle = 0;
storage_ctf.filter.enable = 0;
storage_ctf.filter.value = 3;
storage_ctf.mean = 0;

%Specify directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(varargin,1) == 0
    pathname = uigetdir('Select a directory');

    if isequal(pathname,0) 
        error('Cancel button pressed. No data loaded.');
        return; 
    end;
    storage_ctf.path = pathname;
else
    storage_ctf.path = varargin{1};
end

[storage_ctf.dircell] = get_dircontents(storage_ctf.path,{},{});

%load good/bad file
if size(varargin,2) > 0
    s = load(varargin{2});

    lauf=1;
    dc = {};

    for i=1:length(s.align2d)
        if s.align2d(i).quality ~= 6
            dc{lauf} = s.align2d(i).filename;
            header = tom_reademheader(s.align2d(i).filename);
            header.Header.FocusIncrement = 1;
            tom_writeemheader(s.align2d(i).filename,header.Header);
            lauf=lauf+1;
        end
    end

    storage_ctf.dircell = dc;
end



storage_ctf.imagenumber = 0;
gen_ps(1);

handles.output = hObject;
guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Output Function                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = tom_ctftool_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  save image                                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_ctfsave_Callback(hObject, eventdata, handles)

global storage_ctf;

h=tom_reademheader(get(findobj('Tag','filename_text'),'String'));
h.Header.FocusIncrement = round(str2num(get(findobj('Tag','ctf_defocusvalue'),'String')).*10);
tom_writeemheader(get(findobj('Tag','filename_text'),'String'), h.Header);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  load image                                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_ctfopen_Callback(hObject, eventdata, handles)

global storage_ctf;

image = tom_emreadc;
storage_ctf.Header = image.Header;
storage_ctf.Defocus = image.Header.Defocus./10;
set(findobj('Tag','ctf_defocusvalue'),'String',num2str(storage_ctf.Defocus));

image = tom_smooth(single(image.Value),round(size(image.Value,1).*.05));% 5 % border smoothing
storage_ctf.ps = tom_ps(image);

if storage_ctf.logscale == 1
    calculate_histogram(real(log(storage_ctf.ps)));
else
    calculate_histogram(storage_ctf.ps);
end

set(findobj('Tag','filename_text'),'String',[storage_ctf.Header.Pathname storage_ctf.Header.Filename]);
render_ps();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Reset Histogram                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ctf_histo_reset_Callback(hObject, eventdata, handles)

global storage_ctf;

if storage_ctf.logscale == 1
    calculate_histogram(real(log(storage_ctf.ps)));
else
    calculate_histogram(storage_ctf.ps);
end
render_ps();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Set Histogram                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_ctf_histoset_Callback(hObject, eventdata, handles)

global storage_ctf;

waitforbuttonpress;
point1 = get(gca,'CurrentPoint');    % button down detected
finalRect = rbbox;                   % return figure units
point2 = get(gca,'CurrentPoint');    % button up detected
point1 = point1(1,1:2);              % extract x and y
point2 = point2(1,1:2);
p1 = min(point1,point2);             % calculate locations
offset = abs(point1-point2);         % and dimensions
x = [p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)];
storage_ctf.powerScale=[x(1) x(2)];
set(findobj('Tag','ctf_histo_low'),'String',storage_ctf.powerScale(1));
set(findobj('Tag','ctf_histo_high'),'String',storage_ctf.powerScale(2));
set(gca,'Xlim',[x(1) x(2)]);

render_ps();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Set Histogram Manually                                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_ctf_histo_setmanually_Callback(hObject, eventdata, handles)

global storage_ctf;

min=str2num(get(findobj('Tag','ctf_histo_low'),'String'));
max=str2num(get(findobj('Tag','ctf_histo_high'),'String'));

if max>min
    storage_ctf.powerScale = [min max];
    set(findobj('Tag','ctf_histogram'),'Xlim',[min max]);
else
    set(findobj('Tag','histogram_low'),'String',num2str(storage_ctf.powerScale(1)));
    set(findobj('Tag','histogram_high'),'String',num2str(storage_ctf.powerScale(2)));
end
render_ps();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Pick Histogram Low                                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_ctf_picklow_Callback(hObject, eventdata, handles)

global storage_ctf;

[x y] = ginput(1);
clickedimage = get(gca,'Tag');

if x > 0 & y > 0 &  x < size(storage_ctf.ps,2) & y < size(storage_ctf.ps,1) & strcmp(clickedimage,'ctf_image') == 1

    if storage_ctf.logscale == 1
        val = real(log(storage_ctf.ps(round(x),round(y))));
    else
        val = storage_ctf.ps(round(x),round(y));
    end
    
    if val >= storage_ctf.powerScale(2)
        errordlg('lower limit must be greater than upper limit!','Histrogram error');
    else
        storage_ctf.powerScale = [val storage_ctf.powerScale(2)];
        set(findobj('Tag','ctf_histo_low'),'String',num2str(val));
        render_ps();
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Pick Histogram High                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_ctf_pickhigh_Callback(hObject, eventdata, handles)

global storage_ctf;

[x y] = ginput(1);
clickedimage = get(gca,'Tag');

if x > 0 & y > 0 &  x < size(storage_ctf.ps,2) & y < size(storage_ctf.ps,1) & strcmp(clickedimage,'ctf_image') == 1

    val = storage_ctf.ps(round(x),round(y));

    if val <= storage_ctf.powerScale(1)
        errordlg('upper limit must be higher than lower limit!','Histrogram error');
    else
        storage_ctf.powerScale = [storage_ctf.powerScale(1) val];
        set(findobj('Tag','ctf_histo_high'),'String',num2str(val));
        render_ps();
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Defocus value                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ctf_defocusvalue_Callback(hObject, eventdata, handles)

global storage_ctf;

storage_ctf.Defocus = str2num(get(hObject,'String'));
render_circle();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  button CTF fit                                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_ctffit_Callback(hObject, eventdata, handles)

global storage_ctf; 

binning = str2num(get(handles.ctftool_binning,'String'));
avg = str2num(get(handles.ctftool_average,'String'));
lowcutoff = str2num(get(handles.ctftool_lowcutoff,'String'));
norings = str2num(get(handles.ctftool_norings,'String'));

if get(findobj('Tag','ctffilter_enable'),'Value') == 1
    filterval = str2num(get(findobj('Tag','ctffilter_value'),'String'));
else
    filterval = 1;
end

if get(findobj('Tag','checkbox_ctfverbose'),'Value') == 1
    [Dz, success] = tom_ctffitter(get(findobj('Tag','filename_text'),'String'), str2num(get(findobj('Tag','ctf_defocusvalue'),'String')), filterval,avg, binning, 1, norings, lowcutoff);
else
    [Dz, success] = tom_ctffitter(get(findobj('Tag','filename_text'),'String'), str2num(get(findobj('Tag','ctf_defocusvalue'),'String')), filterval, avg, binning, 0, norings, lowcutoff);
end

if success == 1
    set(findobj('Tag','statusbox'),'BackgroundColor',[0 1 0]);
    storage_ctf.Defocus = round(Dz / 10);
    set(findobj('Tag','ctf_defocusvalue'),'String',num2str(storage_ctf.Defocus));
    render_circle();
else
    set(findobj('Tag','statusbox'),'BackgroundColor',[1 0 0]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  show full circles                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_ctfcircle_Callback(hObject, eventdata, handles)

global storage_ctf;

set(findobj('Tag','button_ctf_halfcircle'),'Value',0);
storage_ctf.semicircle = 0;
render_circle();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  show half circles                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_ctf_halfcircle_Callback(hObject, eventdata, handles)

global storage_ctf;

set(findobj('Tag','button_ctfcircle'),'Value',0);
storage_ctf.semicircle = 1;
render_circle();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  center view                                                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_ctfcenter_Callback(hObject, eventdata, handles)

global storage_ctf;

set(findobj('Tag','button_viewfull'),'Value',0);
storage_ctf.displaycenter = 1;
render_ps();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  full view                                                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_viewfull_Callback(hObject, eventdata, handles)

global storage_ctf;

set(findobj('Tag','button_ctfcenter'),'Value',0);
storage_ctf.displaycenter = 0;
render_ps();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  linear scale                                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_scalelinear_Callback(hObject, eventdata, handles)

global storage_ctf;

set(findobj('Tag','button_scalelog'),'Value',0);
storage_ctf.logscale = 0;
calculate_histogram(storage_ctf.ps);
render_ps();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  log scale                                                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_scalelog_Callback(hObject, eventdata, handles)

global storage_ctf;

set(findobj('Tag','button_scalelinear'),'Value',0);
storage_ctf.logscale = 1;
calculate_histogram(real(log(storage_ctf.ps)));
render_ps();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  filter enable                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ctffilter_enable_Callback(hObject, eventdata, handles)

global storage_ctf;

storage_ctf.filter.enable = get(hObject,'Value');
storage_ctf.filter.value = round(str2num(get(findobj('Tag','ctffilter_value'),'String')));


% if storage_ctf.filter.enable == 1
%     storage_ctf.ps_orig = storage_ctf.ps;
%     storage_ctf.ps = tom_filter(storage_ctf.ps,storage_ctf.filter.value,'quadr','real');
% else
%     storage_ctf.ps = storage_ctf.ps_orig;
% end

render_ps();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  filter low                                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ctffilter_value_Callback(hObject, eventdata, handles)

global storage_ctf;

storage_ctf.filter.value = round(str2num(get(hObject,'String')));

if get(findobj('Tag','ctffilter_enable'),'Value') == 1
    storage_ctf.ps = tom_filter(storage_ctf.ps,storage_ctf.filter.value,'quadr','real');
    render_ps();
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  file slider                                                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ctf_fileslider_Callback(hObject, eventdata, handles)

global storage_ctf;

position = round(get(hObject,'Value'));
gen_ps(position);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  show image                                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_showimage_Callback(hObject, eventdata, handles)

global storage_ctf;

show_image();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  show red circles                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkbox_redcircles_Callback(hObject, eventdata, handles)

render_circle();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  show blue circles                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkbox_bluecircles_Callback(hObject, eventdata, handles)

render_circle();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                     %%
%%                                                                     %%
%%  Helper functions                                                   %%
%%                                                                     %%
%%                                                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  show original image in extra figure                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function show_image()

global storage_ctf;

if get(findobj('Tag','button_showimage'),'Value') == 1
    if isempty(findobj('Tag','ctf_origimage'))
        figure('Tag','ctf_origimage','Name','Original Image (2x binned)');axesobj = axes('Tag','ctf_origimageaxes');
    else
        axesobj = findobj('Tag','ctf_origimageaxes');
        axes(axesobj);
    end
    im = tom_emreadc(get(findobj('Tag','filename_text'),'String'),'binning',2);
    tom_imagesc(im.Value);
    set(axesobj,'Tag','ctf_origimageaxes');
else
    delete(findobj('Tag','ctf_origimage'));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  generate powerspectrum from image                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gen_ps(number)

global storage_ctf;
binning = str2num(get(findobj('Tag','ctftool_binning'),'String'));
split = str2num(get(findobj('Tag','ctftool_average'),'String'));
set(findobj('Tag','statusbox'),'BackgroundColor',[0.702 0.702 0.702]);
image = tom_emreadc([storage_ctf.path '/' storage_ctf.dircell{number}],'binning',binning);

storage_ctf.Header = image.Header;
storage_ctf.Header.Objectpixelsize = storage_ctf.Header.Objectpixelsize;
storage_ctf.Defocus = image.Header.Defocus./10;
set(findobj('Tag','ctf_defocusvalue'),'String',num2str(storage_ctf.Defocus));

image = tom_smooth(single(image.Value),size(image.Value,1)./32);
storage_ctf.ps = tom_image2ctf(double(image),0,0,split,1);%tom_ps(image);

if storage_ctf.logscale == 1
    calculate_histogram(real(log(storage_ctf.ps)));
else
    calculate_histogram(storage_ctf.ps);
end

set(findobj('Tag','filename_text'),'String',[storage_ctf.Header.Pathname storage_ctf.Header.Filename]);
render_ps();

storage_ctf.imagenumber = number;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display powerspectrum                                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function render_ps()

%warning('off');

global storage_ctf;

if storage_ctf.displaycenter == 1
    %cut out central part of power spectrum
    size_ps = size(storage_ctf.ps);
    center = floor(size_ps ./ 2 + 1);
    width = floor(center ./ 2 - 1);
    ps_out = storage_ctf.ps(center(1)-width(1):center(1)+width(1)-1,center(2)-width(2):center(2)+width(2)-1);
else
    ps_out = storage_ctf.ps;
end

if get(findobj('Tag','ctffilter_enable'),'Value') == 1
    ps_out = tom_filter(ps_out,storage_ctf.filter.value,'quadr','real');
end

%display power spectrum logarithmic
if storage_ctf.logscale == 1
    ps_out = real(log(ps_out));
end

tmpobj = findobj('Tag','ctf_image');
axes(tmpobj);
imagesc(ps_out',storage_ctf.powerScale);colormap gray;axis image;
set(tmpobj, 'Tag', 'ctf_image');
render_circle();
show_image();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  calculate histogram                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function calculate_histogram(image)

global storage_ctf;

[mean max min std] = tom_dev(image,'noinfo');
min = mean-3.*std;
max = mean+3.*std;
storage_ctf.mean = mean;

[h,n] = tom_hist3d(image);
h = 200 .* h ./ (100.*size(image,1) .* size(image,2));
axesobj = findobj('Tag','ctf_histogram');
axes(axesobj); 
bar(n,h);
axis auto;
set(axesobj,'Tag','ctf_histogram');
set(findobj('Tag','ctf_histo_low'),'String',num2str(min));
set(findobj('Tag','ctf_histo_high'),'String',num2str(max));
storage_ctf.powerScale = [min, max];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  plot ctf                                                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_ctf(image)

global storage_ctf;

if storage_ctf.displaycenter == 1
    %cut out central part of power spectrum
    size_ps = size(storage_ctf.ps);
    center = floor(size_ps ./ 2 + 1);
    width = floor(center ./ 2 - 1);
    ps_out = storage_ctf.ps(center(1)-width(1):center(1)+width(1)-1,center(2)-width(2):center(2)+width(2)-1);
else
    size_ps = size(storage_ctf.ps);
    center = floor(size_ps ./ 2 + 1);
    width = floor(center - 1);
    ps_out = storage_ctf.ps;
end

%display power spectrum logarithmic
if storage_ctf.logscale == 1
    ps_out = real(log(ps_out));
end

if get(findobj('Tag','ctffilter_enable'),'Value') == 1
    ps_out = tom_filter(ps_out,storage_ctf.filter.value,'quadr','real');
end

%circular integration over ps
integ_ps = tom_cart2polar(ps_out);
integ_ps = sum(integ_ps,2);
tmpobj = findobj('Tag','ctf_graph');
axes(tmpobj);
cla;

%plot theoretical ctf
Dz = storage_ctf.Defocus./1000;
pix_size = storage_ctf.Header.Objectpixelsize./10;
voltage = storage_ctf.Header.Voltage./1000;

if storage_ctf.displaycenter == 1
    pixs = size(ps_out,1).*2;
else
    pixs = size(ps_out,1);
end

Cs = storage_ctf.Header.Cs;
alpha = 0.02;
Cc = 2.2;
deltaE = 0.8;
ctf_out = tom_myctf(Dz, pix_size, voltage, pixs, Cs, alpha, Cc, deltaE, 0);
hold on;
plot(tom_norm(integ_ps,3),'g','LineWidth',1.5);
hold off;
set(tmpobj,'Tag','ctf_graph','XLim',[1 width(1)]);

storage_ctf.ctf_out = ctf_out;

%plot ctf from header
Dz = storage_ctf.Header.Defocus./10000;
ctf_out = tom_myctf(Dz, pix_size, voltage, pixs, Cs, alpha, Cc, deltaE, 0, 0);
hold on;
plot(ctf_out,'--b','LineWidth',1.5);
hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  draw circle                                                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function render_circle()

global storage_ctf;

defocus = storage_ctf.Defocus./1000;
pix_size = storage_ctf.Header.Objectpixelsize./10;
voltage = storage_ctf.Header.Voltage./1000;
pixs = size(storage_ctf.ps,1);
Cs = storage_ctf.Header.Cs;



tmpobj = findobj('Tag','ctf_image');
axes(tmpobj);

hold on;

if storage_ctf.displaycenter == 1
    Center= ((size(storage_ctf.ps,1)./4)+1) + ((size(storage_ctf.ps,2)./4)+1)*sqrt(-1);
else
    Center= ((size(storage_ctf.ps,1)./2)+1) + ((size(storage_ctf.ps,2)./2)+1)*sqrt(-1);
end

Gridpt = 200;

delete(findobj('Tag','ctfcircle'));

%draw red circles
if get(findobj('Tag','checkbox_redcircles'),'Value') == 1
    out=tom_ctfzero(defocus,pix_size, voltage,pixs,Cs);

    for i=1:size(out,2)
        Radius = out(i);
        if storage_ctf.semicircle == 1
            [u,v]=circle(Center,Radius,Gridpt,'semi-circle');
        else
            [u,v]=circle(Center,Radius,Gridpt,'circle');
        end
        line(u,v,'LineWidth',1,'LineStyle','-.','Color',[1 0 0],'Tag','ctfcircle');%red dark
    end

end

%draw blue circles
if get(findobj('Tag','checkbox_bluecircles'),'Value') == 1
    out_blue=tom_ctfzero(storage_ctf.Header.Defocus./10000,pix_size, voltage,pixs,Cs);
    for i=1:size(out_blue,2)
        Radius = out_blue(i);
        if storage_ctf.semicircle == 1
            [u,v]=circle(Center,Radius,Gridpt,'semi-circle2');
        else
            [u,v]=circle(Center,Radius,Gridpt,'circle');
        end
        line(u,v,'LineWidth',1,'LineStyle','-.','Color',[0 0 1],'Tag','ctfcircle');%blue
    end
end;

hold off;
plot_ctf();


function [X, Y] = circle(w,r,n,What)
%   This function is used to calculate a matrix of a 1/2 circle (between pi/2 and -pi/2) or a circle.
%   Syntax: [X,Y]=circle(w,r,n,What)
%       Input:
%           w: Center of the 1/2 circle. W must be a complex number as w=x + yi 
%              (x and y are the coordinate)
%           r: Radius of the circle
%           n: nb of point. 
%           What: case 'semi-circle': matrix of a 1/2 circle (between pi/2 and -pi/2)
%                 case 'circle':matrix of a circle
%       Output:
%           X: it is a matrix of coordinate to draw the circle
%           Y: it is a matrix of coordinate to draw the circle
%
w1 = real(w);
w2 = imag(w);
for k = 1:n
    switch What
        case 'semi-circle'
            t = -k*(pi)/n;
            X(k) = w1 + r*sin(t);
            Y(k) = w2 + r*cos(t);
       case 'semi-circle2'
            t = -k*(pi)/n;
            X(k) = w1 - r*sin(t);
            Y(k) = w2 - r*cos(t);
        case 'circle'             
            t = k*pi/n;
            X(k) = w1 + r*cos(t);
            Y(k) = w2 + r*sin(t);
            X(n+k) = w1 - r*cos(t);
            Y(n+k) = w2 - r*sin(t);
            X(2*n+1) = X(1);
            Y(2*n+1) = Y(1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  ctf calculation                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y= tom_myctf(Dz, pix_size, voltage, pixs, Cs, alpha, Cc, deltaE, infoflag, plotflag)
%TOM_CTF calculates and plots CTF (pure phase contrast)
%
%   ctf_out = tom_ctf(Dz, pix_size, voltage, pixs, Cs, alpha, Cc, deltaE)
%
%   The one-dimensional CTF is calculated and plotted. In case no figure is
%   displayed open one (see also example). The function can be used to get
%   a quick overview of the CTF for cterian imaging conditions.
%
%PARAMETERS
%   INPUT
%   Dz       : Defocus (<0 underfocus, >0 overfocus) (in \mu m);
%   pix_size : pixel size (in nm) (default: 0.72 nm)
%   voltage  : accelerating Voltage (in keV) (default: 300 kV)
%   pixs     : Number of pixels used (default: 2048)
%   Cs       : sperical aberration in mm(default: 2 mm)
%   alpha    : illumination aperture in mrad (default 0.02 - reasonable for FEG)
%   Cc       : chrmatic aberration in mm -(default 2.2 mm)
%   deltaE   : energy width in eV (default 0.8 eV)
%   infoflag : display axes information (default 1)
%
%   OUTPUT
%   ctf_out  : vector of dim pixs/2 containig the plotted ctf (if pixs is 
%                                    not even: dim = ceil(pixs/2))
%              the ouput is beeing plotted where the absyssis are inverse
%              pixels.
%EXAMPLE
%   figure;
%   ctf = tom_ctf(-4,1, 120,1024,2);
%
% SEE ALSO
%    TOM_CREATE_CTF TOM_FOURIER TOM_IFOURIER 
%
%    Copyright (c) 2004
%    TOM toolbox for Electron Tomography
%    Max-Planck-Institute for Biochemistry
%    Dept. Molecular Structural Biology
%    82152 Martinsried, Germany
%    http://www.biochem.mpg.de/tom
%
%
%   10/09/02 FF

global xdata;

error(nargchk(1,10,nargin));

if nargin < 10
    plotflag = 1;
end
if nargin<9
    infoflag = 1;
end

if nargin<8
  deltaE=0.8; 
end;
if nargin<7
  Cc=2.2/1000; 
else
   Cc=Cc/1000; 
end;
if nargin<6
  alpha=0.02/1000; 
else
   alpha=alpha/1000; 
end;
if nargin<5
  Cs=2*10^(-3); 
else
   Cs=Cs*10^(-3); 
end;
if nargin<4
    pixs=2048;
end;
if nargin<3
    voltage=300000; 
else
    voltage = voltage * 1000;
end;
if nargin <2
    pix_size=0.72*10^(-9);
else
    pix_size=pix_size*10^(-9);
end;
Dz=Dz*10^(-6); %from micrometer to meter
Dzn=Dz*1000000; %for display
Csn=Cs*1000;%for display
voltagen=voltage/1000;%for display
voltagest=voltage*(1+voltage/1022000); %for relativistic calc
lambda=sqrt(150.4/voltagest)*10^-10;
q=0:1/(pixs*pix_size):1/(2*pix_size);% von, Increment, Nyqvist
xdata = q;
nyqvist = 2*pix_size*10^9;;
y = sin( pi/2* (Cs*lambda.^3.*q.^4 - 2*Dz*lambda*q.^2) );
% new function: include also envelope function:
% 1) spatial coherence
%  alpha: illumination aperture - assume 0.02mrad
Ks = exp(- ((pi*(Cs.*lambda.^2.*q.^3 - Dz.*q)*alpha).^2)/log(2));
% 2) temporal coherence
delta = Cc*deltaE/voltage;
%Kt = exp(- (pi*lambda*delta*q.^2/2));% alter according to Reimer
Kt = exp(- (pi*lambda*delta*q.^2/(4*log(2))).^2);
K = Kt.*Ks;
%1st zero (approximately for high defocus)
ctf_zero = sqrt(lambda*abs(Dz)*10^18);


if plotflag == 1
    %rings
    plot(0:floor(pixs/2),K.*y,'Color',[1 0 0],'LineWidth',1.5);

    hold on; 

    %envelope function
    plot(0:floor(pixs/2),K,'Color',[0 0 1],'LineWidth',1.5);hold off;

    grid on;
    axis([0 floor(pixs/2) -1.1 1.1] );

    if infoflag == 1
        title(['CTF for: Defocus:',num2str(Dzn),'\mum, Voltage:',num2str(voltagen),'kV, C_s:',num2str(Csn),...
                'mm, Nyqvist:', num2str(nyqvist),'nm, 1st zero:', num2str(ctf_zero,2),'nm']);
        ylabel('CTF');
        xlabel('Frequency');
    end
end
y = K.*y;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  get all the em-files in a directory                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dircell, flagcell_out] = get_dircontents(directory, dircell, flagcell)

%Get list of EM-Files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = waitbar(0,'Getting directory list...');

dirlist = dir(directory);
files = size(dirlist,1);
j = size(dircell,2);
if j == 0
    j = 1;
end

%sort out all the em-files
for i = 1:files
    if dirlist(i).isdir == 0
        if isempty(strmatch(dirlist(i).name,dircell))
            if tom_isemfile([directory '/' dirlist(i).name]) == 1
                dircell{j} = dirlist(i).name;
                try
                    flagcell_out{j} = flagcell{j};
                catch
                    flagcell_out{j} = 1;
                end
                j = size(dircell,2) + 1;
            end
        end
    end
    waitbar(i./files,h,[num2str(i), ' of ', num2str(files), ' files scanned.']);
end

close(h);

if size(dircell,2) == 0
    errordlg('No EM files could be found in this directory.');
    return;
end

set(findobj('Tag','ctf_fileslider'),'Max',size(dircell,2),'SliderStep',[1./size(dircell,2) 1./size(dircell,2)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Unused Callbacks/Create functions                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ctf_defocusvalue_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ctf_histo_low_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ctf_histo_high_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ctf_histo_low_Callback(hObject, eventdata, handles)
function ctf_histo_high_Callback(hObject, eventdata, handles)
function ctffilter_value_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ctfbandpass_high_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ctf_fileslider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function checkbox_ctfverbose_Callback(hObject, eventdata, handles)
function ctftool_binning_Callback(hObject, eventdata, handles)
function ctftool_binning_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ctftool_average_Callback(hObject, eventdata, handles)
function ctftool_average_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ctftool_lowcutoff_Callback(hObject, eventdata, handles)
function ctftool_lowcutoff_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ctftool_norings_Callback(hObject, eventdata, handles)
function ctftool_norings_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


