function varargout = tom_ctftool2(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_ctftool2_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_ctftool2_OutputFcn, ...
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


% -------------------------------------------------------------------------
% Opening function
% -------------------------------------------------------------------------
function tom_ctftool2_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

set(handles.tom_ctftool2,'Currentaxes',handles.axes_psd);
axis off;

set(handles.viewpanel,'SelectionChangeFcn',@helper_selectionchange);

handles.clim_micrograph = [];
handles.clim_psd = [];
handles.clim_enhancedpsd = [];
handles.clim_modelhalf = [];
handles.clim_modelquadrant = [];
handles.micrograph = [];
handles.psd = [];
handles.enhancedpsd = [];
handles.rendermode = 'micrograph';
handles.st.DeltafU = 0;
handles.st.DeltafV = 0;
handles.st.AzimuthalAngle = 0;
handles.st.CTFmodelhalf = zeros(256,256);
handles.st.CTFmodelquadrant = zeros(256,256);
handles.st.zeros = [];
handles.line = zeros(7,1);
handles.astigline = 0;
handles.xlim_val = [0 0];
handles.ylim_val = [0 0];
handles.path = '';
handles.dircell = {};
handles.output = [];
handles.resflag = false;
handles.helpercircle = [0 0];


if nargin > 3
    handles.micrograph = varargin{1};
    set(handles.edit_dv,'String',handles.micrograph.Header.Defocus./10);
    set(handles.edit_du,'String',handles.micrograph.Header.Defocus./10);
    set(handles.edit_astigangle,'String',handles.micrograph.Header.AstigmatismAngle);
    handles = helper_rendering(handles);
end

guidata(hObject, handles);

if nargin > 3
    set(handles.button_exit,'Visible','on');
    set(handles.edit_filename,'Enable','off');
    set(handles.button_browse,'Enable','off');
    set(handles.slider_file,'Enable','off');
    set(handles.text_dir,'Enable','off');
    set(handles.button_save,'Enable','off');
    handles.resflag = true;
    uiwait(handles.tom_ctftool2);
end


% -------------------------------------------------------------------------
% Output function
% -------------------------------------------------------------------------
function varargout = tom_ctftool2_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;
if handles.resflag == true
    delete(handles.tom_ctftool2);
end


% -------------------------------------------------------------------------
% edit filename
% -------------------------------------------------------------------------
function edit_filename_Callback(hObject, eventdata, handles)

handles = helper_getdircontents(handles);
handles = helper_resetcache(handles);
handles = helper_loadfile(handles);
handles = helper_rendering(handles);

guidata(hObject, handles);


% -------------------------------------------------------------------------
% button browse
% -------------------------------------------------------------------------
function button_browse_Callback(hObject, eventdata, handles)

[FileName,PathName] = uigetfile({'*.em'},'Select file to open');

if FileName == 0
    return;
end

set(handles.edit_filename,'String',[PathName FileName]);
handles = helper_getdircontents(handles);

handles = helper_resetcache(handles);
handles = helper_loadfile(handles);
handles = helper_rendering(handles);

guidata(hObject, handles);


% -------------------------------------------------------------------------
% file slider
% -------------------------------------------------------------------------
function slider_file_Callback(hObject, eventdata, handles)

path = fileparts(get(handles.edit_filename,'String'));

set(handles.edit_filename,'String',[path filesep handles.dircell{get(handles.slider_file,'Value')}]);

handles = helper_resetcache(handles);
handles = helper_loadfile(handles);
handles = helper_rendering(handles);

guidata(hObject, handles);


% -------------------------------------------------------------------------
% button save
% -------------------------------------------------------------------------
function button_save_Callback(hObject, eventdata, handles)

handles.micrograph.Header.Defocus = str2double(get(handles.edit_dv,'String')).*10;
handles.micrograph.Header.FocusIncrement = str2double(get(handles.edit_du,'String')).*10;
handles.micrograph.Header.AstigmatismAngle = str2double(get(handles.astigangle,'String'));
tom_writeemheader(handles.micrograph.Header.Filename,handles.micrograph.Header);

guidata(hObject, handles);


% -------------------------------------------------------------------------
% button exit
% -------------------------------------------------------------------------
function button_exit_Callback(hObject, eventdata, handles)

handles.output = [str2double(get(handles.edit_dv,'String')) str2double(get(handles.edit_du,'String')) str2double(get(handles.edit_astigangle,'String'))];
guidata(hObject, handles);
uiresume(handles.tom_ctftool2);


% -------------------------------------------------------------------------
% psd size
% -------------------------------------------------------------------------
function edit_psdsize_Callback(hObject, eventdata, handles)

handles = helper_calcpsd(handles);

guidata(hObject, handles);


% -------------------------------------------------------------------------
% show center
% -------------------------------------------------------------------------
function checkbox_center_Callback(hObject, eventdata, handles)

handles = helper_center_rendering(handles);

guidata(hObject, handles);


% -------------------------------------------------------------------------
% show circles
% -------------------------------------------------------------------------
function checkbox_circles_Callback(hObject, eventdata, handles)

handles = helper_renderellipses(handles);

guidata(hObject, handles);


% -------------------------------------------------------------------------
% show contrast tool
% -------------------------------------------------------------------------
function checkbox_contrast_Callback(hObject, eventdata, handles)

if get(handles.checkbox_contrast,'Value') == 1
    handles.contrasttool = imcontrast(handles.axes_psd);
else
    delete(handles.contrasttool);
end

guidata(hObject, handles);


% -------------------------------------------------------------------------
% defocus V
% -------------------------------------------------------------------------
function edit_dv_Callback(hObject, eventdata, handles)

handles = helper_calczeros(handles);
handles = helper_renderellipses(handles);

guidata(hObject, handles);


% -------------------------------------------------------------------------
% defocus U
% -------------------------------------------------------------------------
function edit_du_Callback(hObject, eventdata, handles)

handles = helper_calczeros(handles);
handles = helper_renderellipses(handles);

guidata(hObject, handles);


% -------------------------------------------------------------------------
% astigmatism angle
% -------------------------------------------------------------------------
function edit_astigangle_Callback(hObject, eventdata, handles)

handles = helper_calczeros(handles);
handles = helper_renderellipses(handles);

guidata(hObject, handles);


% -------------------------------------------------------------------------
% min freq
% -------------------------------------------------------------------------
function edit_minfreq_Callback(hObject, eventdata, handles)

if strcmp(handles.rendermode,'enhanced PSD')
    handles.enhancedpsd = [];
    handles = helper_rendering(handles);
end

handles = helper_renderellipses(handles);

guidata(hObject, handles);


% -------------------------------------------------------------------------
% max freq
% -------------------------------------------------------------------------
function edit_maxfreq_Callback(hObject, eventdata, handles)

if strcmp(handles.rendermode,'enhanced PSD')
    handles.enhancedpsd = [];
    handles = helper_rendering(handles);
end

handles = helper_renderellipses(handles);

guidata(hObject, handles);


% -------------------------------------------------------------------------
% enhance min
% -------------------------------------------------------------------------
function edit_enhancemin_Callback(hObject, eventdata, handles)

if strcmp(handles.rendermode,'enhanced PSD')
    handles.enhancedpsd = [];
    handles = helper_rendering(handles);
end

guidata(hObject, handles);


% -------------------------------------------------------------------------
% enhance max
% -------------------------------------------------------------------------
function edit_enhancemax_Callback(hObject, eventdata, handles)

if strcmp(handles.rendermode,'enhanced PSD')
    handles.enhancedpsd = [];
    handles = helper_rendering(handles);
end

guidata(hObject, handles);


% -------------------------------------------------------------------------
% enhance weight
% -------------------------------------------------------------------------
function edit_enhanceweight_Callback(hObject, eventdata, handles)


guidata(hObject, handles);


% -------------------------------------------------------------------------
% Button Fit
% -------------------------------------------------------------------------
function button_fit_Callback(hObject, eventdata, handles)

if isempty(handles.psd)
    handles = helper_calcpsd(handles);
end

Dz = str2double(get(handles.edit_dv,'String')).*10;
voltage = handles.micrograph.Header.Voltage./1000;
objectpixelsize = handles.micrograph.Header.Objectpixelsize;
Cs = handles.micrograph.Header.Cs;
min_freq = str2double(get(handles.edit_minfreq,'String'));
max_freq = str2double(get(handles.edit_maxfreq,'String'));
Ca = 2;
enhance_filter_min = str2double(get(handles.edit_enhancemin,'String'));
enhance_filter_max = str2double(get(handles.edit_enhancemax,'String'));
enhance_weight = str2double(get(handles.edit_enhanceweight,'String'));

ctfmodelsize = str2double(get(handles.edit_psdsize,'String'));

handles.st = tom_xmipp_adjust_ctf(handles.psd,Dz,voltage,objectpixelsize,ctfmodelsize,Cs,min_freq,max_freq,Ca,enhance_filter_min,enhance_filter_max,enhance_weight);

set(handles.edit_dv,'String',num2str(round(handles.st.DeltafV./10)));
set(handles.edit_du,'String',num2str(round(handles.st.DeltafU./10)));
set(handles.edit_astigangle,'String',num2str(round(handles.st.AzimuthalAngle)));

handles = helper_rendering(handles);
guidata(hObject, handles);


% -------------------------------------------------------------------------
% Button Correct
% -------------------------------------------------------------------------
function button_correct_Callback(hObject, eventdata, handles)

im_corr = tom_xmipp_ctf_correct_phase(handles.micrograph.Value,handles.st,'leave',0);
[PATHSTR,NAME,EXT] = fileparts(handles.micrograph.Header.Filename);

tom_emwrite([handles.path filesep NAME '_corrected.em'],im_corr);

guidata(hObject, handles);


% -------------------------------------------------------------------------
% helper: load file
% -------------------------------------------------------------------------
function handles = helper_loadfile(handles)

fn = get(handles.edit_filename,'String');

if tom_isemfile(fn)
    handles.micrograph = tom_emreadc(fn);
elseif tom_isspiderfile(fn)
    handles.micrograph = tom_spiderread(fn);
else
    error('File format not recognized.');
end
set(handles.edit_dv,'String',handles.micrograph.Header.Defocus./10);
set(handles.edit_du,'String',handles.micrograph.Header.Defocus./10);
set(handles.edit_astigangle,'String',handles.micrograph.Header.AstigmatismAngle);


% -------------------------------------------------------------------------
% helper: calculate psd
% -------------------------------------------------------------------------
function handles = helper_calcpsd(handles)

handles.micrograph.Value = tom_smooth(handles.micrograph.Value, size(handles.micrograph.Value,1)/32);
handles.psd = tom_calc_periodogram(handles.micrograph.Value,str2double(get(handles.edit_psdsize,'String')));


% -------------------------------------------------------------------------
% helper: calculate enhanced PSD
% -------------------------------------------------------------------------
function handles = helper_calcenhancedpsd(handles)

if isempty(handles.psd)
    handles.psd = tom_calc_periodogram(handles.micrograph.Value,str2double(get(handles.edit_psdsize,'String')));
end
handles.enhancedpsd = tom_xmipp_psd_enhance(handles.psd,true,true,str2double(get(handles.edit_enhancemin,'String')),str2double(get(handles.edit_enhancemax,'String')),0.02,0.025,0.5);


% -------------------------------------------------------------------------
% helper: rendering
% -------------------------------------------------------------------------
function handles = helper_rendering(handles)


set(handles.tom_ctftool2,'Currentaxes',handles.axes_psd);
axis auto;
handles.xlim_val = [0 0];
handles.ylim_val = [0 0];

switch handles.rendermode

    case 'micrograph'

        imagesc(tom_bin(handles.micrograph.Value,2)');

        if ~isempty(handles.clim_micrograph)
            set(handles.axes_psd,'Clim',handles.clim_micrograph);
        end

        
    case 'PSD'
        if isempty(handles.psd)
            handles = helper_calcpsd(handles);
        end
        imagesc(log10(fftshift(handles.psd))');

        if ~isempty(handles.clim_psd)
            set(handles.axes_psd,'Clim',handles.clim_psd);
        end
        handles = helper_renderellipses(handles);

        
    case 'enhanced PSD'
        if isempty(handles.enhancedpsd)
            handles = helper_calcenhancedpsd(handles);
        end
        imagesc(handles.enhancedpsd');

        if ~isempty(handles.clim_enhancedpsd)
            set(handles.axes_psd,'Clim',handles.clim_enhancedpsd);
        end
        handles = helper_renderellipses(handles);


    case 'ctfmodel half'
        imagesc(handles.st.CTFmodelhalf');
        if ~isempty(handles.clim_modelhalf)
            set(handles.axes_psd,'Clim',handles.clim_modelhalf);
        end
        handles = helper_renderellipses(handles);

    case 'ctfmodel quadrant'
        imagesc(handles.st.CTFmodelquadrant');
        if ~isempty(handles.clim_modelquadrant)
            set(handles.axes_psd,'Clim',handles.clim_modelquadrant);
        end
        handles = helper_renderellipses(handles);
    
end

axis off;axis ij;colormap gray;

handles = helper_center_rendering(handles);

if get(handles.checkbox_contrast,'Value') == 1
    handles.contrasttool = imcontrast(handles.axes_psd);
end

handles = helper_plotctf(handles);


% -------------------------------------------------------------------------
% helper: View mode change callback function
% -------------------------------------------------------------------------
function handles = helper_selectionchange(hObject,eventdata)

handles = guidata(hObject);

handles.rendermode = get(eventdata.NewValue,'String');

switch get(eventdata.OldValue,'String')
    
    case 'micrograph'
        handles.clim_micrograph = get(handles.axes_psd,'Clim');
    case 'PSD'
        handles.clim_psd = get(handles.axes_psd,'Clim');
    case 'enhanced PSD'
        handles.clim_enhancedpsd = get(handles.axes_psd,'Clim');
    case 'ctfmodel half'
        handles.clim_modelhalf = get(handles.axes_psd,'Clim');
    case 'ctfmodel quadrant'
        handles.clim_modelquadrant = get(handles.axes_psd,'Clim');
end

handles = helper_rendering(handles);

guidata(hObject, handles);


% -------------------------------------------------------------------------
% helper: file changed, reset cache
% -------------------------------------------------------------------------
function handles = helper_resetcache(handles)

sz = str2double(get(handles.edit_psdsize,'String'));
handles.psd = [];
handles.enhancedpsd = [];
handles.st.DeltafU = 0;
handles.st.DeltafV = 0;
handles.st.AzimuthalAngle = 0;
handles.st.CTFmodelhalf = zeros(sz,sz);
handles.st.CTFmodelquadrant = zeros(sz,sz);
handles.st.zeros = [];


% -------------------------------------------------------------------------
% helper: center the rendering
% -------------------------------------------------------------------------
function handles = helper_center_rendering(handles)

if get(handles.checkbox_center,'Value') == 1
    xlim_val = get(handles.axes_psd,'Xlim');
    ylim_val = get(handles.axes_psd,'Ylim');
    
    handles.xlim_val = xlim_val;
    handles.ylim_val = ylim_val;
    
    center = xlim_val(2)./4+1;
    xlim_val(1) = xlim_val(1)+center;
    xlim_val(2) = xlim_val(2)-center;    
    
    center = ylim_val(2)./4+1;
    ylim_val(1) = ylim_val(1)+center;
    ylim_val(2) = ylim_val(2)-center;    
    
    set(handles.axes_psd,'Xlim',xlim_val);
    set(handles.axes_psd,'Ylim',ylim_val);
    
else
    if sum(handles.xlim_val) > 0
        set(handles.axes_psd,'Xlim',handles.xlim_val);
        set(handles.axes_psd,'Ylim',handles.ylim_val);
    end
end


% -------------------------------------------------------------------------
% helper: render ellipses
% -------------------------------------------------------------------------
function handles = helper_renderellipses(handles)

if get(handles.checkbox_circles,'Value') == 0 || strcmp(handles.rendermode,'micrograph') == true
    for i=1:7
        try
            delete(handles.line(i));
        end
    end
    try
        delete(handles.astigline);
    end
    try
        delete(handles.helpercircle(1));
        delete(handles.helpercircle(2));
    end

    return;
end

sz = str2double(get(handles.edit_psdsize,'String'));

set(handles.tom_ctftool2,'Currentaxes',handles.axes_psd);

if ~isempty(handles.st.zeros)
    try
        delete(handles.astigline);
    end

    for i=1:7
        try
            delete(handles.line(i));
            handles.line(i) = 0;
        end
        handles.line(i) = line(handles.st.zeros(1,:,i),handles.st.zeros(2,:,i),'Color',[1 0 0]);

    end
end

try
    delete(handles.helpercircle(1));
    delete(handles.helpercircle(2));
end
    
min_freq = str2double(get(handles.edit_minfreq,'String'));
max_freq = str2double(get(handles.edit_maxfreq,'String'));
handles.helpercircle(1) = tom_ellipse('Axes',[min_freq.*sz min_freq.*sz],'Position',[(sz./2+1) (sz./2+1)],'Major','off','Minor','off','Angle',0,'Type','patch');
handles.helpercircle(2) = tom_ellipse('Axes',[max_freq.*sz max_freq.*sz],'Position',[(sz./2+1) (sz./2+1)],'Major','off','Minor','off','Angle',0,'Type','patch');
set(handles.helpercircle(1),'EdgeColor',[0 1 0]);
set(handles.helpercircle(2),'EdgeColor',[0 1 0]);


% -------------------------------------------------------------------------
% helper: get directory contents
% -------------------------------------------------------------------------
function handles = helper_getdircontents(handles)

filename = get(handles.edit_filename,'String');
[path,fn,ext] = fileparts(filename);
if strcmp(handles.path,path) ~= true
    handles.dircell = tom_HT_getdircontents(path,{'em'},1);
    handles.path = path;
end
numfiles = length(handles.dircell);
pos = strfind(handles.dircell, [fn ext]);
[i,j] = find(~cellfun(@isempty, pos));
set(handles.slider_file,'Min',1,'Max',numfiles,'SliderStep',[1./(numfiles-1) 1./(numfiles-1)],'Value',j(j==1));


% -------------------------------------------------------------------------
% helper: calculate ctf zeros
% -------------------------------------------------------------------------
function handles = helper_calczeros(handles)

sz = str2double(get(handles.edit_psdsize,'String'));
DeltafU = str2double(get(handles.edit_du,'String'));
DeltafV = str2double(get(handles.edit_dv,'String'));
AzimuthalAngle = str2double(get(handles.edit_astigangle,'String'));
handles.st.zeros = tom_xmipp_ctf_zeros(DeltafU*10,DeltafV*10,AzimuthalAngle,handles.micrograph.Header.Voltage/1000,handles.micrograph.Header.Objectpixelsize,sz,handles.micrograph.Header.Cs);


% -------------------------------------------------------------------------
% helper: calc ctf
% -------------------------------------------------------------------------
function handles = helper_plotctf(handles)
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


plotflag = 1;
infoflag = 0;
deltaE=0.8; 
Cc=0.0022; 

alpha=2e-5; 
Cs=handles.micrograph.Header.Cs.*10^(-3); 

pixs=handles.micrograph.Header.Size(1)./2;
voltage = handles.micrograph.Header.Voltage;
pix_size=(handles.micrograph.Header.Objectpixelsize)*10^(-10);
Du = str2double(get(handles.edit_du,'String'));
Dz = str2double(get(handles.edit_dv,'String'));

if Du > 0
    Dz = (Dz + Du) ./2;
end

if Dz == 0
    Dz = handles.micrograph.Header.Defocus./10;
end
Dz = Dz ./ 1000;
Dz = Dz *10^(-6); %from micrometer to meter

Dzn=Dz*1000000; %for display
Csn=Cs*1000;%for display
voltagen=voltage/1000;%for display
voltagest=voltage*(1+voltage/1022000); %for relativistic calc
lambda=sqrt(150.4/voltagest)*10^-10;
q=0:1/(pixs*pix_size):1/(2*pix_size);% von, Increment, Nyqvist
%xdata = q;
nyqvist = 2*pix_size*10^9;
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
    set(0,'CurrentFigure',handles.tom_ctftool2);
    set(handles.tom_ctftool2,'Currentaxes',handles.axes_ctf);
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

    if strcmp(handles.rendermode,'PSD')
    elseif strcmp(handles.rendermode,'enhanced PSD') || strcmp(handles.rendermode,'ctfmodel half') || strcmp(handles.rendermode,'ctfmodel quadrant')
%        plot(sum(tom_cart2polar(handles.enhanced_psd),2),'--b','LineWidth',1.5);
    end
    
end
%y = K.*y;




% -------------------------------------------------------------------------
% Create functions
% -------------------------------------------------------------------------
function edit_filename_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_enhanceweight_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_enhancemax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_enhancemin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_maxfreq_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_minfreq_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_astigangle_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_du_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_dv_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_psdsize_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function slider_file_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


