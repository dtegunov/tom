function varargout = tom_av2_alignment(varargin)
%TOM_AV2_ALIGNMENT is a GUI for 2D alignment of particles
%
%   varargout = tom_av2_alignment(varargin)
%
%PARAMETERS
%
%  INPUT
%   varargin            ...
%  
%  OUTPUT
%   varargout           ...
%
%EXAMPLE
%   ... = tom_av2_alignment(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by AK 01/25/06
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
                   'gui_OpeningFcn', @tom_av2_alignment_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_av2_alignment_OutputFcn, ...
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
function tom_av2_alignment_OpeningFcn(hObject, eventdata, handles, varargin)


global storage_align2dgui;


%initialise storage structure
storage_align2dgui = struct();
storage_align2dgui.parstruct = struct();

storage_align2dgui.input.stack = '';
storage_align2dgui.input.alignfile = '';
storage_align2dgui.input.refstack = '';

storage_align2dgui.output.stack = '';
storage_align2dgui.output.alignfile = '';
storage_align2dgui.output.refstack = '';

storage_align2dgui.iterations.alignment = 5;
storage_align2dgui.iterations.refinement = 5;
storage_align2dgui.demomode.alignment = 0;
storage_align2dgui.demomode.multiref = 0;

storage_align2dgui.filter.align.off = 0;
storage_align2dgui.filter.align.defaults = 0;
storage_align2dgui.filter.align.times = '';
storage_align2dgui.filter.align.bandpass.enable = 0;
storage_align2dgui.filter.align.bandpass.low = '';
storage_align2dgui.filter.align.bandpass.high = '';
storage_align2dgui.filter.align.bandpass.smooth = '';
storage_align2dgui.filter.align.kernel.enable = 0;
storage_align2dgui.filter.align.kernel.real = 0;
storage_align2dgui.filter.align.kernel.fourier = 0;
storage_align2dgui.filter.align.kernel.quadr = 0;
storage_align2dgui.filter.align.kernel.circ = 0;
storage_align2dgui.filter.align.kernel.radius = '';

storage_align2dgui.filter.classify.off = 0;
storage_align2dgui.filter.classify.defaults = 0;
storage_align2dgui.filter.classify.times = '';
storage_align2dgui.filter.classify.bandpass.enable = 0;
storage_align2dgui.filter.classify.bandpass.low = '';
storage_align2dgui.filter.classify.bandpass.high = '';
storage_align2dgui.filter.classify.bandpass.smooth = '';
storage_align2dgui.filter.classify.kernel.enable = 0;
storage_align2dgui.filter.classify.kernel.real = 0;
storage_align2dgui.filter.classify.kernel.fourier = 0;
storage_align2dgui.filter.classify.kernel.quadr = 0;
storage_align2dgui.filter.classify.kernel.circ = 0;
storage_align2dgui.filter.classify.kernel.radius = '';

storage_align2dgui.mask.align.off = 0;
storage_align2dgui.mask.align.defaults = 0;
storage_align2dgui.mask.align.sphere.enable = 0;
storage_align2dgui.mask.align.sphere.radius = '';
storage_align2dgui.mask.align.sphere.sigma = '';
storage_align2dgui.mask.align.sphere.center.x = '';
storage_align2dgui.mask.align.sphere.center.y = '';
storage_align2dgui.mask.align.rect.enable = 0;
storage_align2dgui.mask.align.rect.radius.x = '';
storage_align2dgui.mask.align.rect.radius.y = '';
storage_align2dgui.mask.align.rect.sigma = '';
storage_align2dgui.mask.align.rect.center.x = '';
storage_align2dgui.mask.align.rect.center.y = '';

storage_align2dgui.mask.classify1.off = 0;
storage_align2dgui.mask.classify1.defaults = 0;
storage_align2dgui.mask.classify1.sphere.enable = 0;
storage_align2dgui.mask.classify1.sphere.radius = '';
storage_align2dgui.mask.classify1.sphere.sigma = '';
storage_align2dgui.mask.classify1.sphere.center.x = '';
storage_align2dgui.mask.classify1.sphere.center.y = '';
storage_align2dgui.mask.classify1.rect.enable = 0;
storage_align2dgui.mask.classify1.rect.radius.x = '';
storage_align2dgui.mask.classify1.rect.radius.y = '';
storage_align2dgui.mask.classify1.rect.sigma = '';
storage_align2dgui.mask.classify1.rect.center.x = '';
storage_align2dgui.mask.classify1.rect.center.y = '';

storage_align2dgui.mask.classify2.off = 0;
storage_align2dgui.mask.classify2.defaults = 0;
storage_align2dgui.mask.classify2.sphere.enable = 0;
storage_align2dgui.mask.classify2.sphere.radius = '';
storage_align2dgui.mask.classify2.sphere.sigma = '';
storage_align2dgui.mask.classify2.sphere.center.x = '';
storage_align2dgui.mask.classify2.sphere.center.y = '';
storage_align2dgui.mask.classify2.rect.enable = 0;
storage_align2dgui.mask.classify2.rect.radius.x = '';
storage_align2dgui.mask.classify2.rect.radius.y = '';
storage_align2dgui.mask.classify2.rect.sigma = '';
storage_align2dgui.mask.classify2.rect.center.x = '';
storage_align2dgui.mask.classify2.rect.center.y = '';

storage_align2dgui.mask.cc_rot.off = 0;
storage_align2dgui.mask.cc_rot.defaults = 0;
storage_align2dgui.mask.cc_rot.sphere.enable = 0;
storage_align2dgui.mask.cc_rot.sphere.radius = '';
storage_align2dgui.mask.cc_rot.sphere.sigma = '';
storage_align2dgui.mask.cc_rot.sphere.center.x = '';
storage_align2dgui.mask.cc_rot.sphere.center.y = '';
storage_align2dgui.mask.cc_rot.rect.enable = 0;
storage_align2dgui.mask.cc_rot.rect.radius.x = '';
storage_align2dgui.mask.cc_rot.rect.radius.y = '';
storage_align2dgui.mask.cc_rot.rect.sigma = '';
storage_align2dgui.mask.cc_rot.rect.center.x = '';
storage_align2dgui.mask.cc_rot.rect.center.y = '';

storage_align2dgui.mask.cc_trans.off = 0;
storage_align2dgui.mask.cc_trans.defaults = 0;
storage_align2dgui.mask.cc_trans.sphere.enable = 0;
storage_align2dgui.mask.cc_trans.sphere.radius = '';
storage_align2dgui.mask.cc_trans.sphere.sigma = '';
storage_align2dgui.mask.cc_trans.sphere.center.x = '';
storage_align2dgui.mask.cc_trans.sphere.center.y = '';
storage_align2dgui.mask.cc_trans.rect.enable = 0;
storage_align2dgui.mask.cc_trans.rect.radius.x = '';
storage_align2dgui.mask.cc_trans.rect.radius.y = '';
storage_align2dgui.mask.cc_trans.rect.sigma = '';
storage_align2dgui.mask.cc_trans.rect.center.x = '';
storage_align2dgui.mask.cc_trans.rect.center.y = '';

storage_align2dgui.mask.off = 0;
storage_align2dgui.mask.defaults = 0;

storage_align2dgui.parallelmode = 0;
storage_align2dgui.muratflag = 0;
if nargin > 3
    storage_align2dgui.muratflag = 1;
    storage_align2dgui.particlesize = varargin{1};
end


% Choose default command line output for tom_av2_alignment
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes tom_av2_alignment wait for user response (see UIRESUME)
% uiwait(handles.figure1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Output Function                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = tom_av2_alignment_OutputFcn(hObject, eventdata, handles) 
% Get default command line output from handles structure
varargout{1} = handles.output;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Call the parallel options gui                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_align2d_changeparallel_Callback(hObject, eventdata, handles)

global storage_align2dgui;

if ~isfield(storage_align2dgui.parstruct,'jobmanager')
    storage_align2dgui.parstruct = tom_parallelsettings();
else
    storage_align2dgui.parstruct = tom_parallelsettings(storage_align2dgui.parstruct);
end

infostring = ['host: ' storage_align2dgui.parstruct.jobmanager];
infostring = strvcat(infostring,['number of tasks: ' num2str(storage_align2dgui.parstruct.number_of_tasks)]);
infostring = strvcat(infostring,['workers: ',num2str(storage_align2dgui.parstruct.workers.min) ' - ' num2str(storage_align2dgui.parstruct.workers.max)]);
infostring = strvcat(infostring,['timeout: ',num2str(storage_align2dgui.parstruct.timeout)]);

set(findobj('Tag','align2d_parallelsettings'),'String',infostring);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Call the parallel options gui                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkbox_align2d_parallelmode_Callback(hObject, eventdata, handles)

global storage_align2dgui;

storage_align2dgui.parallelmode = get(hObject,'Value');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                       Input Callbacks                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  directly edit input stack                                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_align2d_instack_Callback(hObject, eventdata, handles)

global storage_align2dgui;

filename = get(hObject, 'String');

if ischar(filename)
    if tom_isemfile(filename) ~= 1
        errordlg('File is not an EM file.','File Error');
    else
        storage_align2dgui.input.stack = filename;
        loadstack();
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  browse for input stack                                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_align2d_instackbrowse_Callback(hObject, eventdata, handles)

global storage_align2dgui;

[filename, pathname] = uigetfile({'*.em';'*.vol'}, 'Pick an EM-file');
if ischar(filename)
    if tom_isemfile([pathname '/' filename]) ~= 1
        errordlg('File is not an EM file.','File Error');
    else
        storage_align2dgui.input.stack = [pathname, filename];
        set(findobj('Tag','input_align2d_instack'),'String',storage_align2dgui.input.stack);
        loadstack();
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  directly edit input alignment file                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_align2d_inalign_Callback(hObject, eventdata, handles)

global storage_align2dgui;

filename = get(hObject, 'String');

if ischar(filename)
    storage_align2dgui.input.alignfile = filename;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  browse for input alignment file                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_align2d_inalignbrowse_Callback(hObject, eventdata, handles)

global storage_align2dgui;

[filename, pathname] = uigetfile({'*.mat'}, 'Pick a MAT-file');
if ischar(filename)
    storage_align2dgui.input.alignfile = [pathname, filename];
    set(findobj('Tag','input_align2d_inalign'),'String',storage_align2dgui.input.alignfile);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  directly edit input reference stack                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_align2d_inref_Callback(hObject, eventdata, handles)

global storage_align2dgui;

filename = get(hObject, 'String');

if ischar(filename)
    if tom_isemfile(filename) ~= 1
        errordlg('File is not an EM file.','File Error');
    else
        storage_align2dgui.input.refstack = filename;
        loadref();
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  browse for input reference stack                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_align2d_inrefbrowse_Callback(hObject, eventdata, handles)

global storage_align2dgui;

[filename, pathname] = uigetfile({'*.em';'*.vol'}, 'Pick an EM-file');
if ischar(filename)
    if tom_isemfile([pathname '/' filename]) ~= 1
        errordlg('File is not an EM file.','File Error');
    else
        storage_align2dgui.input.refstack = [pathname, filename];
        set(findobj('Tag','input_align2d_inref'),'String',storage_align2dgui.input.refstack);
        loadref();
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                      Output Callbacks                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  directly edit output stack                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_align2d_outstack_Callback(hObject, eventdata, handles)

global storage_align2dgui;

filename = get(hObject, 'String');

if ischar(filename)
    storage_align2dgui.output.stack = filename;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  browse for output stack                                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_align2d_outstackbrowse_Callback(hObject, eventdata, handles)

global storage_align2dgui;

[filename, pathname] = uiputfile({'*.em';'*.vol'}, 'Pick an EM-file');
if ischar(filename)
    storage_align2dgui.output.stack = [pathname, filename];
    set(findobj('Tag','input_align2d_outstack'),'String',storage_align2dgui.output.stack);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  directly edit output alignment file                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_align2d_outalign_Callback(hObject, eventdata, handles)

global storage_align2dgui;

filename = get(hObject, 'String');

if ischar(filename)
    storage_align2dgui.output.alignfile = filename;
end  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  browse for output alignment file                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_align2d_outalignbrowse_Callback(hObject, eventdata, handles)

global storage_align2dgui;

[filename, pathname] = uiputfile({'*.mat'}, 'Pick an MAT-file');
if ischar(filename)
    storage_align2dgui.output.alignfile = [pathname, filename];
    set(findobj('Tag','input_align2d_outalign'),'String',storage_align2dgui.output.alignfile);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  directly edit output reference stack                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_align2d_outref_Callback(hObject, eventdata, handles)

global storage_align2dgui;

filename = get(hObject, 'String');

if ischar(filename)
    storage_align2dgui.output.refstack = filename;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  browse for output reference stack                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_align2d_outrefbrowse_Callback(hObject, eventdata, handles)

global storage_align2dgui;

[filename, pathname] = uiputfile({'*.em';'*.vol'}, 'Pick an EM-file');
if ischar(filename)
    storage_align2dgui.output.refstack = [pathname, filename];
    set(findobj('Tag','input_align2d_outref'),'String',storage_align2dgui.output.refstack);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                      Filter Callbacks                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  filter alignment                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radio_align2d_filter_alignment_Callback(hObject, eventdata, handles)

if get(hObject,'Value') == 1
    load_filtervalues(get_filtertype());
else
    set(hObject,'Value',1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  filter classification                                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radio_align2d_filter_classification_Callback(hObject, eventdata, handles)

if get(hObject,'Value') == 1
    load_filtervalues(get_filtertype());
else
    set(hObject,'Value',1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  filter defaults checkbox                                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkbox_align2d_filter_defaults_Callback(hObject, eventdata, handles)

global storage_align2dgui;

type = get_filtertype();

if get(hObject,'Value') == 1
    storage_align2dgui.filter.(type).defaults = 1;
    set(findobj('Tag','checkbox_align2d_filter_off'),'Enable','off');
    toggle_filtersettings('off');
else
    storage_align2dgui.filter.(type).defaults = 0;
    set(findobj('Tag','checkbox_align2d_filter_off'),'Enable','on');
    if storage_align2dgui.filter.(type).off == 0
        toggle_filtersettings('on');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  filter off checkbox                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkbox_align2d_filter_off_Callback(hObject, eventdata, handles)

global storage_align2dgui;

type = get_filtertype();

if get(hObject,'Value') == 1
    storage_align2dgui.filter.(type).off = 1;
    set(findobj('Tag','checkbox_align2d_filter_defaults'),'Enable','off');
    toggle_filtersettings('off');
else
    storage_align2dgui.filter.(type).off = 0;
    set(findobj('Tag','checkbox_align2d_filter_defaults'),'Enable','on');
    if storage_align2dgui.filter.(type).defaults == 0
        toggle_filtersettings('on');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  filter times                                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_align2d_filter_times_Callback(hObject, eventdata, handles)

global storage_align2dgui;

type = get_filtertype();

val = abs(round(str2num(get(hObject,'String'))));
set(hObject,'String',num2str(val));
storage_align2dgui.filter.(type).times = val;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  filter bandpass radio button                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radio_align2d_filter_bandpass_Callback(hObject, eventdata, handles)

global storage_align2dgui;

type = get_filtertype();

if get(hObject,'Value') == 1
    set(findobj('Tag','radio_align2d_filter_kernel'),'Value', 0);
    storage_align2dgui.filter.(type).bandpass.enable = 1;
    storage_align2dgui.filter.(type).kernel.enable = 0;
else
    set(hObject,'Value', 1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  filter bandpass low                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_align2d_filter_bandpass_low_Callback(hObject, eventdata, handles)

global storage_align2dgui;

type = get_filtertype();

val = abs(round(str2num(get(hObject,'String'))));
set(hObject,'String',num2str(val));
storage_align2dgui.filter.(type).bandpass.low = val;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  filter bandpass high                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_align2d_filter_bandpass_high_Callback(hObject, eventdata, handles)

global storage_align2dgui;

type = get_filtertype();

val = abs(round(str2num(get(hObject,'String'))));
set(hObject,'String',num2str(val));
storage_align2dgui.filter.(type).bandpass.high = val;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  filter bandpass smooth                                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_align2d_filter_bandpass_smooth_Callback(hObject, eventdata, handles)

global storage_align2dgui;

type = get_filtertype();

val = abs(round(str2num(get(hObject,'String'))));
set(hObject,'String',num2str(val));
storage_align2dgui.filter.(type).bandpass.smooth = val;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  filter kernel radio button                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radio_align2d_filter_kernel_Callback(hObject, eventdata, handles)

global storage_align2dgui;

type = get_filtertype();

if get(hObject,'Value') == 1
    set(findobj('Tag','radio_align2d_filter_bandpass'),'Value', 0);
    storage_align2dgui.filter.(type).bandpass.enable = 0;
    storage_align2dgui.filter.(type).kernel.enable = 1;
else
    set(hObject,'Value', 1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  filter kernel real space                                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radio_align2d_filter_kernel_real_Callback(hObject, eventdata, handles)

global storage_align2dgui;

type = get_filtertype();

if get(hObject,'Value') == 1
    set(findobj('Tag','radio_align2d_filter_kernel_fourier'),'Value', 0);
    storage_align2dgui.filter.(type).kernel.real = 1;
    storage_align2dgui.filter.(type).kernel.fourier = 0;
else
    set(hObject,'Value', 1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  filter kernel fourier space                                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radio_align2d_filter_kernel_fourier_Callback(hObject, eventdata, handles)

global storage_align2dgui;

type = get_filtertype();

if get(hObject,'Value') == 1
    set(findobj('Tag','radio_align2d_filter_kernel_real'),'Value', 0);
    storage_align2dgui.filter.(type).kernel.real = 0;
    storage_align2dgui.filter.(type).kernel.fourier = 1;
else
    set(hObject,'Value', 1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  filter kernel circular                                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radio_align2d_filter_kernel_circular_Callback(hObject, eventdata, handles)

global storage_align2dgui;

type = get_filtertype();

if get(hObject,'Value') == 1
    set(findobj('Tag','radio_align2d_filter_kernel_quadratic'),'Value', 0);
    storage_align2dgui.filter.(type).kernel.quadr = 0;
    storage_align2dgui.filter.(type).kernel.circ = 1;
else
    set(hObject,'Value', 1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  filter kernel quadratic                                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radio_align2d_filter_kernel_quadratic_Callback(hObject, eventdata, handles)

global storage_align2dgui;

type = get_filtertype();

if get(hObject,'Value') == 1
    set(findobj('Tag','radio_align2d_filter_kernel_circular'),'Value', 0);
    storage_align2dgui.filter.(type).kernel.quadr = 1;
    storage_align2dgui.filter.(type).kernel.circ = 0;
else
    set(hObject,'Value', 1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  filter kernel radius                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_align2d_filter_kernel_radius_Callback(hObject, eventdata, handles)

global storage_align2dgui;

type = get_filtertype();

val = abs(round(str2num(get(hObject,'String'))));
set(hObject,'String',num2str(val));
storage_align2dgui.filter.(type).kernel.radius = val;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                        Mask Callbacks                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  mask alignment                                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radio_align2d_mask_alignment_Callback(hObject, eventdata, handles)

if get(hObject,'Value') == 1
    load_maskvalues(get_masktype());
else
    set(hObject,'Value',1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  mask classification 1                                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radio_align2d_mask_classification1_Callback(hObject, eventdata, handles)

if get(hObject,'Value') == 1
    load_maskvalues(get_masktype());
else
    set(hObject,'Value',1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  mask classification 2                                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radio_align2d_mask_classification2_Callback(hObject, eventdata, handles)

if get(hObject,'Value') == 1
    load_maskvalues(get_masktype());
else
    set(hObject,'Value',1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  mask cc rot                                                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radio_align2d_mask_cc_rot_Callback(hObject, eventdata, handles)

if get(hObject,'Value') == 1
    load_maskvalues(get_masktype());
else
    set(hObject,'Value',1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  mask cc trans                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radio_align2d_mask_cc_trans_Callback(hObject, eventdata, handles)

if get(hObject,'Value') == 1
    load_maskvalues(get_masktype());
else
    set(hObject,'Value',1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  mask defaults                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkbox_align2d_mask_defaults_Callback(hObject, eventdata, handles)

global storage_align2dgui;

type = get_masktype();

if get(hObject,'Value') == 1
    storage_align2dgui.mask.(type).defaults = 1;
    set(findobj('Tag','checkbox_align2d_mask_off'),'Enable','off');
    toggle_masksettings('off');
else
    storage_align2dgui.mask.(type).defaults = 0;
    set(findobj('Tag','checkbox_align2d_mask_off'),'Enable','on');
    if storage_align2dgui.mask.(type).off == 0
        toggle_masksettings('on');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  mask off                                                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkbox_align2d_mask_off_Callback(hObject, eventdata, handles)

global storage_align2dgui;

type = get_masktype();

if get(hObject,'Value') == 1
    storage_align2dgui.mask.(type).off = 1;
    set(findobj('Tag','checkbox_align2d_mask_defaults'),'Enable','off');
    toggle_masksettings('off');
else
    storage_align2dgui.mask.(type).off = 0;
    set(findobj('Tag','checkbox_align2d_mask_defaults'),'Enable','on');
    if storage_align2dgui.mask.(type).defaults == 0
        toggle_masksettings('on');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  spheremask                                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radio_align2d_mask_sphere_Callback(hObject, eventdata, handles)

global storage_align2dgui;

type = get_masktype();

if get(hObject,'Value') == 1
    set(findobj('Tag','radio_align2d_mask_rect'),'Value', 0);
    storage_align2dgui.mask.(type).rect.enable = 0;
    storage_align2dgui.mask.(type).sphere.enable = 1;
else
    set(hObject,'Value', 1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  spheremask radius                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_align2d_mask_sphere_radius_Callback(hObject, eventdata, handles)

global storage_align2dgui;

type = get_masktype();

val = abs(round(str2num(get(hObject,'String'))));
set(hObject,'String',num2str(val));
storage_align2dgui.mask.(type).sphere.radius = val;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  spheremask sigma                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_align2d_mask_sphere_sigma_Callback(hObject, eventdata, handles)

global storage_align2dgui;

type = get_masktype();

val = abs(round(str2num(get(hObject,'String'))));
set(hObject,'String',num2str(val));
storage_align2dgui.mask.(type).sphere.sigma = val;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  spheremask center x                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_align2d_mask_sphere_centerx_Callback(hObject, eventdata, handles)

global storage_align2dgui;

type = get_masktype();

val = abs(round(str2num(get(hObject,'String'))));
set(hObject,'String',num2str(val));
storage_align2dgui.mask.(type).sphere.center.x = val;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  spheremask center y                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_align2d_mask_sphere_centery_Callback(hObject, eventdata, handles)

global storage_align2dgui;

type = get_masktype();

val = abs(round(str2num(get(hObject,'String'))));
set(hObject,'String',num2str(val));
storage_align2dgui.mask.(type).sphere.center.y = val;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  rectmask                                                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radio_align2d_mask_rect_Callback(hObject, eventdata, handles)

global storage_align2dgui;

type = get_masktype();

if get(hObject,'Value') == 1
    set(findobj('Tag','radio_align2d_mask_sphere'),'Value', 0);
    storage_align2dgui.mask.(type).rect.enable = 1;
    storage_align2dgui.mask.(type).sphere.enable = 0;
else
    set(hObject,'Value', 1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  rectmask sigma                                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_align2d_mask_rect_sigma_Callback(hObject, eventdata, handles)

global storage_align2dgui;

type = get_masktype();

val = abs(round(str2num(get(hObject,'String'))));
set(hObject,'String',num2str(val));
storage_align2dgui.mask.(type).rect.sigma = val;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  rectmask center x                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_align2d_mask_rect_centerx_Callback(hObject, eventdata, handles)

global storage_align2dgui;

type = get_masktype();

val = abs(round(str2num(get(hObject,'String'))));
set(hObject,'String',num2str(val));
storage_align2dgui.mask.(type).rect.center.x = val;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  rectmask center y                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_align2d_mask_rect_centery_Callback(hObject, eventdata, handles)

global storage_align2dgui;

type = get_masktype();

val = abs(round(str2num(get(hObject,'String'))));
set(hObject,'String',num2str(val));
storage_align2dgui.mask.(type).rect.center.y = val;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  rectmask radius x                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_align2d_mask_rect_radiusx_Callback(hObject, eventdata, handles)

global storage_align2dgui;

type = get_masktype();

val = abs(round(str2num(get(hObject,'String'))));
set(hObject,'String',num2str(val));
storage_align2dgui.mask.(type).rect.radius.x = val;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  rectmask radius y                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_align2d_mask_rect_radiusy_Callback(hObject, eventdata, handles)

global storage_align2dgui;

type = get_masktype();

val = abs(round(str2num(get(hObject,'String'))));
set(hObject,'String',num2str(val));
storage_align2dgui.mask.(type).rect.radius.y = val;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                        MISC Callbacks                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  apply filter stack                                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_align2d_applyfilter_Callback(hObject, eventdata, handles)

global storage_align2dgui;

if isfield(storage_align2dgui,'sampleimage')
     if get(findobj('Tag','checkbox_align2d_filter_defaults'),'Value') == 1
        filter_st.Apply = 2;
     elseif get(findobj('Tag','checkbox_align2d_filter_off'),'Value') == 1   
         filter_st.Apply = 0;
     else    
        filter_st.Apply = 1;
        filter_st.Times = str2num(get(findobj('Tag','input_align2d_filter_times'),'String'));                                           
        %bandpass filter
        if get(findobj('Tag','radio_align2d_filter_bandpass'),'Value') == 1
            filter_st.Method='bandpass';
            filter_st.Value = [str2num(get(findobj('Tag','input_align2d_filter_bandpass_low'),'String')) str2num(get(findobj('Tag','input_align2d_filter_bandpass_high'),'String')) str2num(get(findobj('Tag','input_align2d_filter_bandpass_smooth'),'String'))];
        %kernel filter
        elseif get(findobj('Tag','radio_align2d_filter_kernel'),'Value') == 1
            filter_st.Value = [str2num(get(findobj('Tag','input_align2d_filter_kernel_radius'),'String'))];
            if get(findobj('Tag','radio_align2d_filter_kernel_real'),'Value') == 1
                filter_st.Space = 'real';
            elseif get(findobj('Tag','radio_align2d_filter_kernel_fourier'),'Value') == 1
                filter_st.Space = 'fourier';
            else
                errordlg('Specify filter first.');
                return;
            end
            
            if get(findobj('Tag','radio_align2d_filter_kernel_quadratic'),'Value') == 1
                filter_st.Method = 'quadr';
            elseif get(findobj('Tag','radio_align2d_filter_kernel_circular'),'Value') == 1
                filter_st.Method = 'circ';
            else
                errordlg('Specify filter first.');
                return;
            end
        else
            errordlg('Specify filter first.');
            return;
        end
        
     end
     
     im=tom_apply_filter(storage_align2dgui.sampleimage.Value,filter_st); 
     tmpobj = findobj('Tag','im_align2d_sample');
     axes(tmpobj);
     tom_imagesc(im,'noinfo');
     set(tmpobj,'Tag','im_align2d_sample');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  apply filter ref                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_align2d_applyfilter_ref_Callback(hObject, eventdata, handles)

global storage_align2dgui;

if isfield(storage_align2dgui,'samplerefimage')
     if get(findobj('Tag','checkbox_align2d_filter_defaults'),'Value') == 1
        filter_st.Apply = 2;
     elseif get(findobj('Tag','checkbox_align2d_filter_off'),'Value') == 1   
         filter_st.Apply = 0;
     else    
        filter_st.Apply = 1;
        filter_st.Times = str2num(get(findobj('Tag','input_align2d_filter_times'),'String'));                                           
        %bandpass filter
        if get(findobj('Tag','radio_align2d_filter_bandpass'),'Value') == 1
            filter_st.Method='bandpass';
            filter_st.Value = [str2num(get(findobj('Tag','input_align2d_filter_bandpass_low'),'String')) str2num(get(findobj('Tag','input_align2d_filter_bandpass_high'),'String')) str2num(get(findobj('Tag','input_align2d_filter_bandpass_smooth'),'String'))];
        %kernel filter
        elseif get(findobj('Tag','radio_align2d_filter_kernel'),'Value') == 1
            filter_st.Value = [str2num(get(findobj('Tag','input_align2d_filter_kernel_radius'),'String'))];
            if get(findobj('Tag','radio_align2d_filter_kernel_real'),'Value') == 1
                filter_st.Space = 'real';
            elseif get(findobj('Tag','radio_align2d_filter_kernel_fourier'),'Value') == 1
                filter_st.Space = 'fourier';
            else
                errordlg('Specify filter first.');
                return;
            end
            
            if get(findobj('Tag','radio_align2d_filter_kernel_quadratic'),'Value') == 1
                filter_st.Method = 'quadr';
            elseif get(findobj('Tag','radio_align2d_filter_kernel_circular'),'Value') == 1
                filter_st.Method = 'circ';
            else
                errordlg('Specify filter first.');
                return;
            end
        else
            errordlg('Specify filter first.');
            return;
        end
        
     end
     
     im=tom_apply_filter(storage_align2dgui.samplerefimage.Value,filter_st); 
     tmpobj = findobj('Tag','im_align2d_refsample');
     axes(tmpobj);
     tom_imagesc(im,'noinfo');
     set(tmpobj,'Tag','im_align2d_refsample');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  apply mask stack                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_align2d_applymask_Callback(hObject, eventdata, handles)

global storage_align2dgui;

if isfield(storage_align2dgui,'sampleimage')
    if get(findobj('Tag','checkbox_align2d_mask_defaults'),'Value') == 1
        mask_st.Apply = 2;
        mask_st.Value = [storage_align2dgui.particlesize(1) storage_align2dgui.particlesize(2)];
    elseif get(findobj('Tag','checkbox_align2d_mask_off'),'Value') == 1
        mask_st.Apply = 0;
        mask_st.Value = [storage_align2dgui.particlesize(1) storage_align2dgui.particlesize(2)];
    else
        mask_st.Apply = 1;
        if get(findobj('Tag','radio_align2d_mask_sphere'),'Value') == 1
            mask_st.Method = 'sphere';
            mask_st.Value = [storage_align2dgui.particlesize(1) storage_align2dgui.particlesize(2) str2num(get(findobj('Tag','input_align2d_mask_sphere_radius'),'String')) str2num(get(findobj('Tag','input_align2d_mask_sphere_sigma'),'String')) str2num(get(findobj('Tag','input_align2d_mask_sphere_centerx'),'String')) str2num(get(findobj('Tag','input_align2d_mask_sphere_centery'),'String'))];
        elseif get(findobj('Tag','radio_align2d_mask_rectangle'),'Value') == 1
            mask_st.Method = 'rectangle';
            mask_st.Value = [storage_align2dgui.particlesize(1) storage_align2dgui.particlesize(2) str2num(get(findobj('Tag','input_align2d_mask_rect_radiusx'),'String')) str2num(get(findobj('Tag','input_align2d_mask_rect_radiusy'),'String')) str2num(get(findobj('Tag','input_align2d_mask_rect_sigma'),'String')) str2num(get(findobj('Tag','input_align2d_mask_rect_centerx'),'String')) str2num(get(findobj('Tag','input_align2d_mask_rect_centery'),'String'))];
        else
            errordlg('Specify mask first.');
            return;
        end
        
    end
    
    mask=tom_create_mask(mask_st);
    tmpobj = findobj('Tag','im_align2d_sample');
    axes(tmpobj);
    tom_imagesc(storage_align2dgui.sampleimage.Value .*mask,'noinfo');
    set(tmpobj,'Tag','im_align2d_sample');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  apply mask ref                                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_align2d_applymask_ref_Callback(hObject, eventdata, handles)

global storage_align2dgui;

if isfield(storage_align2dgui,'samplerefimage')
    if get(findobj('Tag','checkbox_align2d_mask_defaults'),'Value') == 1
        mask_st.Apply = 2;
        mask_st.Value = [storage_align2dgui.refparticlesize(1) storage_align2dgui.refparticlesize(2)];
    elseif get(findobj('Tag','checkbox_align2d_mask_off'),'Value') == 1
        mask_st.Apply = 0;
        mask_st.Value = [storage_align2dgui.refparticlesize(1) storage_align2dgui.refparticlesize(2)];
    else
        mask_st.Apply = 1;
        if get(findobj('Tag','radio_align2d_mask_sphere'),'Value') == 1
            mask_st.Method = 'sphere';
            mask_st.Value = [storage_align2dgui.refparticlesize(1) storage_align2dgui.refparticlesize(2) str2num(get(findobj('Tag','input_align2d_mask_sphere_radius'),'String')) str2num(get(findobj('Tag','input_align2d_mask_sphere_sigma'),'String')) str2num(get(findobj('Tag','input_align2d_mask_sphere_centerx'),'String')) str2num(get(findobj('Tag','input_align2d_mask_sphere_centery'),'String'))];
        elseif get(findobj('Tag','radio_align2d_mask_rectangle'),'Value') == 1
            mask_st.Method = 'rectangle';
            mask_st.Value = [storage_align2dgui.refparticlesize(1) storage_align2dgui.refparticlesize(2) str2num(get(findobj('Tag','input_align2d_mask_rect_radiusx'),'String')) str2num(get(findobj('Tag','input_align2d_mask_rect_radiusy'),'String')) str2num(get(findobj('Tag','input_align2d_mask_rect_sigma'),'String')) str2num(get(findobj('Tag','input_align2d_mask_rect_centerx'),'String')) str2num(get(findobj('Tag','input_align2d_mask_rect_centery'),'String'))];
        else
            errordlg('Specify mask first.');
            return;
        end
        
    end
    
    mask=tom_create_mask(mask_st);
    tmpobj = findobj('Tag','im_align2d_refsample');
    axes(tmpobj);
    tom_imagesc(storage_align2dgui.samplerefimage.Value .*mask,'noinfo');
    set(tmpobj,'Tag','im_align2d_refsample');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  iterations refinement                                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_align2d_iterations_refinement_Callback(hObject, eventdata, handles)

global storage_align2dgui;

storage_align2dgui.iterations.refinement = abs(round(str2num(get(hObject,'String'))));
set(hObject,'String',num2str(storage_align2dgui.iterations.refinement));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  iterations alignment                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_align2d_iterations_alignment_Callback(hObject, eventdata, handles)

global storage_align2dgui;

storage_align2dgui.iterations.alignment = abs(round(str2num(get(hObject,'String'))));
set(hObject,'String',num2str(storage_align2dgui.iterations.alignment));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  demomode alignment                                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkbox_align2d_demo_alignment_Callback(hObject, eventdata, handles)

global storage_align2dgui;

if get(hObject,'Value') == 1
    storage_align2dgui.demomode.alignment = 1;
else
    storage_align2dgui.demomode.alignment = 0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  demomode multiref                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkbox_align2d_demo_multiref_Callback(hObject, eventdata, handles)

global storage_align2dgui;

if get(hObject,'Value') == 1
    storage_align2dgui.demomode.multiref = 1;
else
    storage_align2dgui.demomode.multiref = 0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  load settings                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_align2d_settings_load_Callback(hObject, eventdata, handles)

global storage_align2dgui;

[filename, pathname] = uigetfile('*.mat', 'Load Current Settings');
if ischar(filename)
    s = load([pathname '/' filename]);

    if isfield(s,'storagesave_alignment2d') ~= 1
        errordlg('This is not a settings file.','File Error');
    else
        
        storage_align2dgui = s.storagesave_alignment2d;
        
        %input / output fields
        set(findobj('Tag','input_align2d_instack'),'String',storage_align2dgui.input.stack);
        set(findobj('Tag','input_align2d_inalign'),'String',storage_align2dgui.input.alignfile);
        set(findobj('Tag','input_align2d_inref'),'String',storage_align2dgui.input.refstack);
        set(findobj('Tag','input_align2d_outstack'),'String',storage_align2dgui.output.stack);
        set(findobj('Tag','input_align2d_outalign'),'String',storage_align2dgui.output.alignfile);
        set(findobj('Tag','input_align2d_outref'),'String',storage_align2dgui.output.refstack);
        loadstack();
        loadref();
        
        %parallel info box
        infostring = ['Host:' storage_align2dgui.parstruct.jobmanager];
        infostring = strvcat(infostring,['number of tasks: ' num2str(storage_align2dgui.parstruct.number_of_tasks)]);
        infostring = strvcat(infostring,['workers: ',num2str(storage_align2dgui.parstruct.workers.min) ' - ' num2str(storage_align2dgui.parstruct.workers.max)]);
        infostring = strvcat(infostring,['timeout: ',num2str(storage_align2dgui.parstruct.timeout)]);
        set(findobj('Tag','align2d_parallelsettings'),'String',infostring);

        if storage_align2dgui.parallelmode == 1
            set(findobj('Tag','checkbox_align2d_parallelmode'),'Value',1);
        else
            set(findobj('Tag','checkbox_align2d_parallelmode'),'Value',0);
        end

        %demo mode
        if storage_align2dgui.demomode.alignment == 1
            set(findobj('Tag','checkbox_align2d_demo_alignment'),'Value',1);
        else
            set(findobj('Tag','checkbox_align2d_demo_alignment'),'Value',0);
        end
        if storage_align2dgui.demomode.multiref == 1
            set(findobj('Tag','checkbox_align2d_demo_multiref'),'Value',1);
        else
            set(findobj('Tag','checkbox_align2d_demo_multiref'),'Value',0);
        end
        
        %iterations
        set(findobj('Tag','input_align2d_iterations_refinement'),'String',num2str(storage_align2dgui.iterations.refinement));
        set(findobj('Tag','input_align2d_iterations_alignment'),'String',num2str(storage_align2dgui.iterations.alignment));
       
        %filter values
        set(findobj('Tag','radio_align2d_filter_alignment'),'Value',1);
        set(findobj('Tag','radio_align2d_filter_classification'),'Value',0);    
        load_filtervalues('align');
        
        %mask values
        set(findobj('Tag','radio_align2d_mask_alignment'),'Value',1);
        set(findobj('Tag','radio_align2d_mask_classification1'),'Value',0);
        set(findobj('Tag','radio_align2d_mask_classification2'),'Value',0);
        set(findobj('Tag','radio_align2d_mask_cc_rot'),'Value',0);
        set(findobj('Tag','radio_align2d_mask_cc_trans'),'Value',0);
        load_maskvalues('align');
        
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  save settings                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_align2d_settings_save_Callback(hObject, eventdata, handles)

global storage_align2dgui;

[filename, pathname] = uiputfile('*.mat', 'Save Current Settings as');
if ischar(filename)
    storagesave_alignment2d = storage_align2dgui;
    save([pathname '/' filename], 'storagesave_alignment2d');
    disp('Settings Saved');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  run alignment                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_align2d_runalignment_Callback(hObject, eventdata, handles)

global storage_align2dgui;

storage_align2dgui.outstruct = struct();

%masks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%align
%%%%%%%%%%%%
%ones
if storage_align2dgui.mask.align.off == 1
    storage_align2dgui.outstruct.mask.align.Apply = 0;
    storage_align2dgui.outstruct.mask.align.Type = 'sphere2d';
    storage_align2dgui.outstruct.mask.align.Value.size_x = storage_align2dgui.particlesize(1);
    storage_align2dgui.outstruct.mask.align.Value.size_y = storage_align2dgui.particlesize(2);
%default values
elseif storage_align2dgui.mask.align.defaults == 1
    storage_align2dgui.outstruct.mask.align.Apply = 2;
    storage_align2dgui.outstruct.mask.align.Type = 'sphere2d';
    storage_align2dgui.outstruct.mask.align.Value.size_x = storage_align2dgui.particlesize(1);
    storage_align2dgui.outstruct.mask.align.Value.size_y = storage_align2dgui.particlesize(2);
%real mask
else
    storage_align2dgui.outstruct.mask.align.Apply = 1;
    %sphere mask
    if storage_align2dgui.mask.align.sphere.enable == 1
        storage_align2dgui.outstruct.mask.align.Type = 'sphere2d';
        storage_align2dgui.outstruct.mask.align.Value.size_x = storage_align2dgui.particlesize(1);
        storage_align2dgui.outstruct.mask.align.Value.size_y = storage_align2dgui.particlesize(2);
        storage_align2dgui.outstruct.mask.align.Value.radius = storage_align2dgui.mask.align.sphere.radius;
        storage_align2dgui.outstruct.mask.align.Value.sigma = storage_align2dgui.mask.align.sphere.sigma;
        storage_align2dgui.outstruct.mask.align.Value.center_x = storage_align2dgui.mask.align.sphere.center.x;
        storage_align2dgui.outstruct.mask.align.Value.center_y = storage_align2dgui.mask.align.sphere.center.y;
        
    %rectangle mask
    elseif storage_align2dgui.mask.align.rect.enable == 1
        storage_align2dgui.outstruct.mask.align.Type = 'rectangle';
        storage_align2dgui.outstruct.mask.align.Value.size_x = storage_align2dgui.particlesize(1);
        storage_align2dgui.outstruct.mask.align.Value.size_y = storage_align2dgui.particlesize(2);
        storage_align2dgui.outstruct.mask.align.Value.radius_x = storage_align2dgui.mask.align.rect.radius.x;
        storage_align2dgui.outstruct.mask.align.Value.radius_y = storage_align2dgui.mask.align.rect.radius.y;
        storage_align2dgui.outstruct.mask.align.Value.sigma = storage_align2dgui.mask.align.rect.sigma;
        storage_align2dgui.outstruct.mask.align.Value.center_x = storage_align2dgui.mask.align.rect.center.x;
        storage_align2dgui.outstruct.mask.align.Value.center_y = storage_align2dgui.mask.align.rect.center.y;
    else
        errordlg('Specify alignment mask!');
        return;
    end
end
storage_align2dgui.outstruct.mask.align.Value.angle_psi = '';
storage_align2dgui.outstruct.mask.align.Value.angle_phi = '';
storage_align2dgui.outstruct.mask.align.Value.angle_theta = '';

%classify1
%%%%%%%%%%%%
%ones
if storage_align2dgui.mask.classify1.off == 1
    storage_align2dgui.outstruct.mask.classify1.Apply = 0;
    storage_align2dgui.outstruct.mask.classify1.Type = 'sphere2d';
    storage_align2dgui.outstruct.mask.classify1.Value.size_x = storage_align2dgui.particlesize(1);
    storage_align2dgui.outstruct.mask.classify1.Value.size_y = storage_align2dgui.particlesize(2);
%default values
elseif storage_align2dgui.mask.classify1.defaults == 1
    storage_align2dgui.outstruct.mask.classify1.Apply = 2;
    storage_align2dgui.outstruct.mask.classify1.Type = 'sphere2d';
    storage_align2dgui.outstruct.mask.classify1.Value.size_x = storage_align2dgui.particlesize(1);
    storage_align2dgui.outstruct.mask.classify1.Value.size_y = storage_align2dgui.particlesize(2);
%real mask
else
    storage_align2dgui.outstruct.mask.classify1.Apply = 1;
    %sphere mask
    if storage_align2dgui.mask.classify1.sphere.enable == 1
        storage_align2dgui.outstruct.mask.classify1.Type = 'sphere2d';
        storage_align2dgui.outstruct.mask.classify1.Value.size_x = storage_align2dgui.particlesize(1);
        storage_align2dgui.outstruct.mask.classify1.Value.size_y = storage_align2dgui.particlesize(2);
        storage_align2dgui.outstruct.mask.classify1.Value.radius = storage_align2dgui.mask.classify1.sphere.radius;
        storage_align2dgui.outstruct.mask.classify1.Value.sigma = storage_align2dgui.mask.classify1.sphere.sigma;
        storage_align2dgui.outstruct.mask.classify1.Value.center_x = storage_align2dgui.mask.classify1.sphere.center.x;
        storage_align2dgui.outstruct.mask.classify1.Value.center_y = storage_align2dgui.mask.classify1.sphere.center.y;
        
    %rectangle mask
    elseif storage_align2dgui.mask.classify1.rect.enable == 1
        storage_align2dgui.outstruct.mask.classify1.Type = 'rectangle';
        storage_align2dgui.outstruct.mask.classify1.Value.size_x = storage_align2dgui.particlesize(1);
        storage_align2dgui.outstruct.mask.classify1.Value.size_y = storage_align2dgui.particlesize(2);
        storage_align2dgui.outstruct.mask.classify1.Value.radius_x = storage_align2dgui.mask.classify1.rect.radius.x;
        storage_align2dgui.outstruct.mask.classify1.Value.radius_y = storage_align2dgui.mask.classify1.rect.radius.y;
        storage_align2dgui.outstruct.mask.classify1.Value.sigma = storage_align2dgui.mask.classify1.rect.sigma;
        storage_align2dgui.outstruct.mask.classify1.Value.center_x = storage_align2dgui.mask.classify1.rect.center.x;
        storage_align2dgui.outstruct.mask.classify1.Value.center_y = storage_align2dgui.mask.classify1.rect.center.y;
    else
        errordlg('Specify Classification 1 mask!');
        return;
    end
end

storage_align2dgui.outstruct.mask.classify1.Value.angle_psi = '';
storage_align2dgui.outstruct.mask.classify1.Value.angle_phi = '';
storage_align2dgui.outstruct.mask.classify1.Value.angle_theta = '';


%classify2
%%%%%%%%%%%%
%ones
if storage_align2dgui.mask.classify2.off == 1
    storage_align2dgui.outstruct.mask.classify2.Apply = 0;
    storage_align2dgui.outstruct.mask.classify2.Type = 'sphere2d';
    storage_align2dgui.outstruct.mask.classify2.Value.size_x = storage_align2dgui.particlesize(1);
    storage_align2dgui.outstruct.mask.classify2.Value.size_y = storage_align2dgui.particlesize(2);
%default values
elseif storage_align2dgui.mask.classify2.defaults == 1
    storage_align2dgui.outstruct.mask.classify2.Apply = 2;
    storage_align2dgui.outstruct.mask.classify2.Type = 'sphere2d';
    storage_align2dgui.outstruct.mask.classify2.Value.size_x = storage_align2dgui.particlesize(1);
    storage_align2dgui.outstruct.mask.classify2.Value.size_y = storage_align2dgui.particlesize(2);
%real mask
else
    storage_align2dgui.outstruct.mask.classify2.Apply = 1;
    %sphere mask
    if storage_align2dgui.mask.classify2.sphere.enable == 1
        storage_align2dgui.outstruct.mask.classify2.Type = 'sphere2d';
        storage_align2dgui.outstruct.mask.classify2.Value.size_x = storage_align2dgui.particlesize(1);
        storage_align2dgui.outstruct.mask.classify2.Value.size_y = storage_align2dgui.particlesize(2);
        storage_align2dgui.outstruct.mask.classify2.Value.radius = storage_align2dgui.mask.classify2.sphere.radius;
        storage_align2dgui.outstruct.mask.classify2.Value.sigma = storage_align2dgui.mask.classify2.sphere.sigma;
        storage_align2dgui.outstruct.mask.classify2.Value.center_x = storage_align2dgui.mask.classify2.sphere.center.x;
        storage_align2dgui.outstruct.mask.classify2.Value.center_y = storage_align2dgui.mask.classify2.sphere.center.y;
        
    %rectangle mask
    elseif storage_align2dgui.mask.classify2.rect.enable == 1
        storage_align2dgui.outstruct.mask.classify2.Type = 'rectangle';
        storage_align2dgui.outstruct.mask.classify2.Value.size_x = storage_align2dgui.particlesize(1);
        storage_align2dgui.outstruct.mask.classify2.Value.size_y = storage_align2dgui.particlesize(2);
        storage_align2dgui.outstruct.mask.classify2.Value.radius_x = storage_align2dgui.mask.classify2.rect.radius.x;
        storage_align2dgui.outstruct.mask.classify2.Value.radius_y = storage_align2dgui.mask.classify2.rect.radius.y;
        storage_align2dgui.outstruct.mask.classify2.Value.sigma = storage_align2dgui.mask.classify2.rect.sigma;
        storage_align2dgui.outstruct.mask.classify2.Value.center_x = storage_align2dgui.mask.classify2.rect.center.x;
        storage_align2dgui.outstruct.mask.classify2.Value.center_y = storage_align2dgui.mask.classify2.rect.center.y;
    else
        errordlg('Specify Classification 2 mask!');
        return;
    end
end

storage_align2dgui.outstruct.mask.classify2.Value.angle_psi = '';
storage_align2dgui.outstruct.mask.classify2.Value.angle_phi = '';
storage_align2dgui.outstruct.mask.classify2.Value.angle_theta = '';


%cc_rot
%%%%%%%%%%%%

sizemask = size(tom_cart2polar(ones(storage_align2dgui.particlesize(1),storage_align2dgui.particlesize(2))));
%ones
if storage_align2dgui.mask.cc_rot.off == 1
    storage_align2dgui.outstruct.mask.cc_rot.Apply = 0;
    storage_align2dgui.outstruct.mask.cc_rot.Type = 'sphere2d';
    storage_align2dgui.outstruct.mask.cc_rot.Value.size_x = sizemask(1);
    storage_align2dgui.outstruct.mask.cc_rot.Value.size_y = sizemask(2);
%default values
elseif storage_align2dgui.mask.cc_rot.defaults == 1
    storage_align2dgui.outstruct.mask.cc_rot.Apply = 2;
    storage_align2dgui.outstruct.mask.cc_rot.Type = 'sphere2d';
    storage_align2dgui.outstruct.mask.cc_rot.Value.size_x = sizemask(1);
    storage_align2dgui.outstruct.mask.cc_rot.Value.size_y = sizemask(2);
%real mask
else
    storage_align2dgui.outstruct.mask.cc_rot.Apply = 1;
    %sphere mask
    if storage_align2dgui.mask.cc_rot.sphere.enable == 1
        storage_align2dgui.outstruct.mask.cc_rot.Type = 'sphere2d';
        storage_align2dgui.outstruct.mask.cc_rot.Value.size_x = sizemask(1);
        storage_align2dgui.outstruct.mask.cc_rot.Value.size_y = sizemask(2);
        storage_align2dgui.outstruct.mask.cc_rot.Value.radius = storage_align2dgui.mask.cc_rot.sphere.radius;
        storage_align2dgui.outstruct.mask.cc_rot.Value.sigma = storage_align2dgui.mask.cc_rot.sphere.sigma;
        storage_align2dgui.outstruct.mask.cc_rot.Value.center_x = storage_align2dgui.mask.cc_rot.sphere.center.x;
        storage_align2dgui.outstruct.mask.cc_rot.Value.center_y = storage_align2dgui.mask.cc_rot.sphere.center.y;
        
    %rectangle mask
    elseif storage_align2dgui.mask.cc_rot.rect.enable == 1
        storage_align2dgui.outstruct.mask.cc_rot.Type = 'rectangle';
        storage_align2dgui.outstruct.mask.cc_rot.Value.size_x = sizemask(1);
        storage_align2dgui.outstruct.mask.cc_rot.Value.size_y = sizemask(2);
        storage_align2dgui.outstruct.mask.cc_rot.Value.radius_x = storage_align2dgui.mask.cc_rot.rect.radius.x;
        storage_align2dgui.outstruct.mask.cc_rot.Value.radius_y = storage_align2dgui.mask.cc_rot.rect.radius.y;
        storage_align2dgui.outstruct.mask.cc_rot.Value.sigma = storage_align2dgui.mask.cc_rot.rect.sigma;
        storage_align2dgui.outstruct.mask.cc_rot.Value.center_x = storage_align2dgui.mask.cc_rot.rect.center.x;
        storage_align2dgui.outstruct.mask.cc_rot.Value.center_y = storage_align2dgui.mask.cc_rot.rect.center.y;
    else
        errordlg('Specify CC Rotation mask!');
        return;
    end
end

storage_align2dgui.outstruct.mask.cc_rot.Value.angle_psi = '';
storage_align2dgui.outstruct.mask.cc_rot.Value.angle_phi = '';
storage_align2dgui.outstruct.mask.cc_rot.Value.angle_theta = '';


%cc_trans
%%%%%%%%%%%%
%ones
if storage_align2dgui.mask.cc_trans.off == 1
    storage_align2dgui.outstruct.mask.cc_trans.Apply = 0;
    storage_align2dgui.outstruct.mask.cc_trans.Type = 'sphere2d';
    storage_align2dgui.outstruct.mask.cc_trans.Value.size_x = storage_align2dgui.particlesize(1);
    storage_align2dgui.outstruct.mask.cc_trans.Value.size_y = storage_align2dgui.particlesize(2);
%default values
elseif storage_align2dgui.mask.cc_trans.defaults == 1
    storage_align2dgui.outstruct.mask.cc_trans.Apply = 2;
    storage_align2dgui.outstruct.mask.cc_trans.Type = 'sphere2d';
    storage_align2dgui.outstruct.mask.cc_trans.Value.size_x = storage_align2dgui.particlesize(1);
    storage_align2dgui.outstruct.mask.cc_trans.Value.size_y = storage_align2dgui.particlesize(2);
%real mask
else
    storage_align2dgui.outstruct.mask.cc_trans.Apply = 1;
    %sphere mask
    if storage_align2dgui.mask.cc_trans.sphere.enable == 1
        storage_align2dgui.outstruct.mask.cc_trans.Type = 'sphere2d';
        storage_align2dgui.outstruct.mask.cc_trans.Value.size_x = storage_align2dgui.particlesize(1);
        storage_align2dgui.outstruct.mask.cc_trans.Value.size_y = storage_align2dgui.particlesize(2);
        storage_align2dgui.outstruct.mask.cc_trans.Value.radius = storage_align2dgui.mask.cc_trans.sphere.radius;
        storage_align2dgui.outstruct.mask.cc_trans.Value.sigma = storage_align2dgui.mask.cc_trans.sphere.sigma;
        storage_align2dgui.outstruct.mask.cc_trans.Value.center_x = storage_align2dgui.mask.cc_trans.sphere.center.x;
        storage_align2dgui.outstruct.mask.cc_trans.Value.center_y = storage_align2dgui.mask.cc_trans.sphere.center.y;
        
    %rectangle mask
    elseif storage_align2dgui.mask.cc_trans.rect.enable == 1
        storage_align2dgui.outstruct.mask.cc_trans.Type = 'rectangle';
        storage_align2dgui.outstruct.mask.cc_trans.Value.size_x = storage_align2dgui.particlesize(1);
        storage_align2dgui.outstruct.mask.cc_trans.Value.size_y = storage_align2dgui.particlesize(2);
        storage_align2dgui.outstruct.mask.cc_trans.Value.radius_x = storage_align2dgui.mask.cc_trans.rect.radius.x;
        storage_align2dgui.outstruct.mask.cc_trans.Value.radius_y = storage_align2dgui.mask.cc_trans.rect.radius.y;
        storage_align2dgui.outstruct.mask.cc_trans.Value.sigma = storage_align2dgui.mask.cc_trans.rect.sigma;
        storage_align2dgui.outstruct.mask.cc_trans.Value.center_x = storage_align2dgui.mask.cc_trans.rect.center.x;
        storage_align2dgui.outstruct.mask.cc_trans.Value.center_y = storage_align2dgui.mask.cc_trans.rect.center.y;
    else
        errordlg('Specify CC Translation mask!');
        return;
    end
end

storage_align2dgui.outstruct.mask.cc_trans.Value.angle_psi = '';
storage_align2dgui.outstruct.mask.cc_trans.Value.angle_phi = '';
storage_align2dgui.outstruct.mask.cc_trans.Value.angle_theta = '';



%filters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%classify
%%%%%%%%%%%%
%ones
if storage_align2dgui.filter.classify.off == 1
    storage_align2dgui.outstruct.filter.classify.Apply = 0;
    storage_align2dgui.outstruct.filter.classify.Type = 'bandpass';
%default values
elseif storage_align2dgui.filter.classify.defaults == 1
    storage_align2dgui.outstruct.filter.classify.Apply = 2;
    storage_align2dgui.outstruct.filter.classify.Type = 'bandpass';
%real mask
else
    storage_align2dgui.outstruct.filter.classify.Apply = 1;
    
    storage_align2dgui.outstruct.filter.classify.Value.times = storage_align2dgui.filter.classify.Value.times;
    %bandpass
    if storage_align2dgui.filter.classify.bandpass.enable == 1
        storage_align2dgui.outstruct.filter.classify.Type = 'bandpass';
        storage_align2dgui.outstruct.filter.classify.Value.low = storage_align2dgui.filter.classify.bandpass.low;
        storage_align2dgui.outstruct.filter.classify.Value.high = storage_align2dgui.filter.classify.bandpass.high;
        storage_align2dgui.outstruct.filter.classify.Value.smooth =storage_align2dgui.filter.classify.bandpass.smooth;
    %real space filter
    elseif storage_align2dgui.filter.classify.kernel.enable == 1
        storage_align2dgui.outstruct.filter.classify.Value.radius = storage_align2dgui.filter.classify.kernel.radius;
        storage_align2dgui.outstruct.filter.classify.Type = 'kernel';
        %circular kernel
        if storage_align2dgui.filter.classify.kernel.circ == 1
            storage_align2dgui.outstruct.filter.classify.Value.method = 'circ';
        %quadratic kernel
        elseif storage_align2dgui.filter.classify.kernel.quadr == 1
            storage_align2dgui.outstruct.filter.classify.Value.method = 'quadr';
        else
            errordlg('Specify Classification Filter!');
            return;
        end
        %fourier space
        if storage_align2dgui.filter.classify.kernel.fourier == 1
            storage_align2dgui.outstruct.filter.classify.Value.space = 'fourier';
        %real space
        elseif storage_align2dgui.filter.classify.kernel.real == 1
            storage_align2dgui.outstruct.filter.classify.Value.space = 'real';
        else
            errordlg('Specify Classification Filter!');
            return;
        end
    else
        errordlg('Specify Classification Filter!');
        return;
    end
end

%align
%%%%%%%%%%%%
%ones
if storage_align2dgui.filter.align.off == 1
    storage_align2dgui.outstruct.filter.align.Apply = 0;
    storage_align2dgui.outstruct.filter.align.Type = 'bandpass';
%default values
elseif storage_align2dgui.filter.align.defaults == 1
    storage_align2dgui.outstruct.filter.align.Apply = 2;
    storage_align2dgui.outstruct.filter.align.Type = 'bandpass';
%real mask
else
    storage_align2dgui.outstruct.filter.align.Apply = 1;
    storage_align2dgui.outstruct.filter.align.Value.times = storage_align2dgui.filter.align.Value.times;
    %bandpass
    if storage_align2dgui.filter.align.bandpass.enable == 1
        storage_align2dgui.outstruct.filter.align.Type = 'bandpass';
        storage_align2dgui.outstruct.filter.align.Value.low = storage_align2dgui.filter.align.bandpass.low;
        storage_align2dgui.outstruct.filter.align.Value.high = storage_align2dgui.filter.align.bandpass.high;
        storage_align2dgui.outstruct.filter.align.Value.smooth = storage_align2dgui.filter.align.bandpass.smooth;
        
    %real space filter
    elseif storage_align2dgui.filter.align.kernel.enable == 1
        storage_align2dgui.outstruct.filter.align.Type = 'kernel';
        storage_align2dgui.outstruct.filter.align.Value.radius = storage_align2dgui.filter.align.kernel.radius;
        %circular kernel
        if storage_align2dgui.filter.align.kernel.circ == 1
            storage_align2dgui.outstruct.filter.align.Value.method = 'circ';
        %quadratic kernel
        elseif storage_align2dgui.filter.align.kernel.quadr == 1
            storage_align2dgui.outstruct.filter.align.Value.method = 'quadr';
        else
            errordlg('Specify Alignment Filter!');
            return;
        end
        %fourier space
        if storage_align2dgui.filter.align.kernel.fourier == 1
            storage_align2dgui.outstruct.filter.align.Value.space = 'fourier';
        %real space
        elseif storage_align2dgui.filter.align.kernel.real == 1
            storage_align2dgui.outstruct.filter.align.Value.space = 'real';
        else
            errordlg('Specify Alignment Filter!');
            return;
        end
    else
        errordlg('Specify Alignment Filter!');
        return;
    end
end


%feed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stack_path = get(findobj('Tag','input_align2d_instack'),'String');
ref_path = get(findobj('Tag','input_align2d_inref'),'String');
align2d = get(findobj('Tag','input_align2d_inalign'),'String');
stack_alg_path = get(findobj('Tag','input_align2d_outstack'),'String');
ref_alg_path = get(findobj('Tag','input_align2d_outref'),'String');
filter_param = storage_align2dgui.outstruct;
if storage_align2dgui.parallelmode == 1
    parallel_param = storage_align2dgui.parstruct;
else
    parallel_param = '';
end
iterations = [str2num(get(findobj('Tag','input_align2d_iterations_alignment'),'String')) str2num(get(findobj('Tag','input_align2d_iterations_refinement'),'String'))];

demo_align = get(findobj('Tag','checkbox_align2d_demo_alignment'),'Value');
% if demo_align == 1
%     demo_align = 2;
% end
demo_multiref = get(findobj('Tag','checkbox_align2d_demo_multiref'),'Value');
% demo = demo_align + demo_multiref;

demo.developer=demo_multiref;
demo.presentation=0;
demo.show_alignment=demo_align;

tmp=get(handles.checkbox_align2d_sub_stacks,'Value');

if (tmp==1)
    sub_stack='sub_stacks';
else
    sub_stack='';
end;

if storage_align2dgui.muratflag == 0
    [align2d new_ref]=tom_av2_align_stack(stack_path,ref_path,align2d,stack_alg_path,ref_alg_path,filter_param,parallel_param,iterations,demo,sub_stack);
else
    flagstruct.stack_path = stack_path;
    flagstruct.ref_path=ref_path;
    flagstruct.align2d = align2d;
    flagstruct.stack_alg_path = stack_alg_path;
    flagstruct.ref_alg_path = ref_alg_path;
    flagstruct.filter_param = filter_param;
    flagstruct.parallel_param = parallel_param;
    flagstruct.iterations = iterations;
    flagstruct.demo = demo;
    flagstruct.sub_stack = sub_stack;
    assignin('base','alignstruct',flagstruct);
end

save(get(findobj('Tag','input_align2d_outalign'),'String'),'align2d');
msgbox('Run Finished');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display input stack in tom_av2_stackbrowser                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_align2d_stackbrowse_Callback(hObject, eventdata, handles)

alignfile = get(findobj('Tag','input_align2d_inalign'),'String');
stackfile = get(findobj('Tag','input_align2d_instack'),'String');
if isempty(stackfile)
    errordlg('Select input stack first!');
    return;
end

if isempty(alignfile)
    tom_av2_stackbrowser(stackfile);
else 
    tom_av2_stackbrowser(stackfile,alignfile);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display input reference in tom_av2_stackbrowser                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_align2d_refbrowse_Callback(hObject, eventdata, handles)

stackfile = get(findobj('Tag','input_align2d_inref'),'String');
if isempty(stackfile)
    errordlg('Select input ref first!');
    return;
end

tom_av2_stackbrowser(stackfile);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  make reference from stack                                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_align2d_stack2ref_Callback(hObject, eventdata, handles)

stackfile = get(findobj('Tag','input_align2d_instack'),'String');
if isempty(stackfile)
    errordlg('Select input stack first!');
    return;
end

refname = tom_av2_stackbrowser(stackfile,'makerefmode');
try
    set(findobj('Tag','input_align2d_inref'),'String',refname);
    loadref();
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  modify reference                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_align2d_refmodify_Callback(hObject, eventdata, handles)

stackfile = get(findobj('Tag','input_align2d_inref'),'String');
if isempty(stackfile)
    errordlg('Select input ref first!');
    return;
end

tom_av2_alignref(stackfile);
loadref();


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
%%  enable/disable filter widgets                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function toggle_filtersettings(mode)

set(findobj('Tag','input_align2d_filter_times'),'Enable',mode);
set(findobj('Tag','radio_align2d_filter_bandpass'),'Enable',mode);
set(findobj('Tag','input_align2d_filter_bandpass_low'),'Enable',mode);
set(findobj('Tag','input_align2d_filter_bandpass_high'),'Enable',mode);
set(findobj('Tag','input_align2d_filter_bandpass_smooth'),'Enable',mode);
set(findobj('Tag','radio_align2d_filter_kernel'),'Enable',mode);
set(findobj('Tag','radio_align2d_filter_kernel_real'),'Enable',mode);
set(findobj('Tag','radio_align2d_filter_kernel_fourier'),'Enable',mode);
set(findobj('Tag','radio_align2d_filter_kernel_quadratic'),'Enable',mode);
set(findobj('Tag','radio_align2d_filter_kernel_circular'),'Enable',mode);
set(findobj('Tag','input_align2d_filter_kernel_radius'),'Enable',mode);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  enable/disable mask widgets                                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function toggle_masksettings(mode)

set(findobj('Tag','radio_align2d_mask_sphere'),'Enable',mode);
set(findobj('Tag','input_align2d_mask_sphere_radius'),'Enable',mode);
set(findobj('Tag','input_align2d_mask_sphere_sigma'),'Enable',mode);
set(findobj('Tag','input_align2d_mask_sphere_centerx'),'Enable',mode);
set(findobj('Tag','input_align2d_mask_sphere_centery'),'Enable',mode);
set(findobj('Tag','radio_align2d_mask_rect'),'Enable',mode);
set(findobj('Tag','input_align2d_mask_rect_radiusx'),'Enable',mode);
set(findobj('Tag','input_align2d_mask_rect_radiusy'),'Enable',mode);
set(findobj('Tag','input_align2d_mask_rect_sigma'),'Enable',mode);
set(findobj('Tag','input_align2d_mask_rect_centerx'),'Enable',mode);
set(findobj('Tag','input_align2d_mask_rect_centery'),'Enable',mode);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  get active filter                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function t = get_filtertype()

if get(findobj('Tag','radio_align2d_filter_alignment'),'Value') == 1
    t = 'align';
elseif get(findobj('Tag','radio_align2d_filter_classification'),'Value') == 1
    t = 'classify';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  get active mask                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function t = get_masktype()

if get(findobj('Tag','radio_align2d_mask_alignment'),'Value') == 1
    t = 'align';
elseif get(findobj('Tag','radio_align2d_mask_classification1'),'Value') == 1
    t = 'classify1';
elseif get(findobj('Tag','radio_align2d_mask_classification2'),'Value') == 1
    t = 'classify2';
elseif get(findobj('Tag','radio_align2d_mask_cc_rot'),'Value') == 1
    t = 'cc_rot';
elseif get(findobj('Tag','radio_align2d_mask_cc_trans'),'Value') == 1
    t = 'cc_trans';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  update values on filter type change                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function load_filtervalues(type)

global storage_align2dgui;
    
toggle_filtersettings('on');
set(findobj('Tag','checkbox_align2d_filter_off'),'Enable','on');
set(findobj('Tag','checkbox_align2d_filter_defaults'),'Enable','on');


if storage_align2dgui.filter.(type).off == 1
    set(findobj('Tag','checkbox_align2d_filter_off'),'Value',1);
    set(findobj('Tag','checkbox_align2d_filter_defaults'),'Enable','off');
    toggle_filtersettings('off');
else
    set(findobj('Tag','checkbox_align2d_filter_off'),'Value',0);
end

if storage_align2dgui.filter.(type).defaults == 1
    set(findobj('Tag','checkbox_align2d_filter_defaults'),'Value',1);
    set(findobj('Tag','checkbox_align2d_filter_off'),'Enable','off');
    toggle_filtersettings('off');
else
    set(findobj('Tag','checkbox_align2d_filter_defaults'),'Value',0);
end

set(findobj('Tag','input_align2d_filter_times'),'String',num2str(storage_align2dgui.filter.(type).times));
set(findobj('Tag','input_align2d_filter_bandpass_low'),'String',num2str(storage_align2dgui.filter.(type).bandpass.low));
set(findobj('Tag','input_align2d_filter_bandpass_high'),'String',num2str(storage_align2dgui.filter.(type).bandpass.high));
set(findobj('Tag','input_align2d_filter_bandpass_smooth'),'String',num2str(storage_align2dgui.filter.(type).bandpass.smooth));
set(findobj('Tag','input_align2d_filter_kernel.radius'),'String',num2str(storage_align2dgui.filter.(type).kernel.radius));

if storage_align2dgui.filter.(type).bandpass.enable == 1
    set(findobj('Tag','radio_align2d_filter_bandpass'),'Value',1);
else
    set(findobj('Tag','radio_align2d_filter_bandpass'),'Value',0);
end

if storage_align2dgui.filter.(type).kernel.enable == 1
    set(findobj('Tag','radio_align2d_filter_kernel'),'Value',1);
else
    set(findobj('Tag','radio_align2d_filter_kernel'),'Value',0);
end

if storage_align2dgui.filter.(type).kernel.real == 1
    set(findobj('Tag','radio_align2d_filter_kernel_real'),'Value',1);
else
    set(findobj('Tag','radio_align2d_filter_kernel_real'),'Value',0);
end

if storage_align2dgui.filter.(type).kernel.fourier == 1
    set(findobj('Tag','radio_align2d_filter_kernel_fourier'),'Value',1);
else
    set(findobj('Tag','radio_align2d_filter_kernel_fourier'),'Value',0);
end

if storage_align2dgui.filter.(type).kernel.quadr == 1
    set(findobj('Tag','radio_align2d_filter_kernel_quadratic'),'Value',1);
else
    set(findobj('Tag','radio_align2d_filter_kernel_quadratic'),'Value',0);
end

if storage_align2dgui.filter.(type).kernel.circ == 1
    set(findobj('Tag','radio_align2d_filter_kernel_circular'),'Value',1);
else
    set(findobj('Tag','radio_align2d_filter_kernel_circular'),'Value',0);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  update values on filter type change                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function load_maskvalues(type)

global storage_align2dgui;
    
toggle_masksettings('on');
set(findobj('Tag','checkbox_align2d_mask_defaults'),'Enable','on');
set(findobj('Tag','checkbox_align2d_mask_off'),'Enable','on');

if storage_align2dgui.mask.(type).off == 1
    set(findobj('Tag','checkbox_align2d_mask_off'),'Value',1);
    set(findobj('Tag','checkbox_align2d_mask_defaults'),'Enable','off');
    toggle_masksettings('off');
else
    set(findobj('Tag','checkbox_align2d_mask_off'),'Value',0);
end

if storage_align2dgui.mask.(type).defaults == 1
    set(findobj('Tag','checkbox_align2d_mask_defaults'),'Value',1);
    set(findobj('Tag','checkbox_align2d_mask_off'),'Enable','off');
    toggle_masksettings('off');
else
    set(findobj('Tag','checkbox_align2d_mask_defaults'),'Value',0);
end

if storage_align2dgui.mask.(type).sphere.enable == 1
    set(findobj('Tag','radio_align2d_mask_sphere'),'Value',1);
else
    set(findobj('Tag','radio_align2d_mask_sphere'),'Value',0);
end

if storage_align2dgui.mask.(type).rect.enable == 1
    set(findobj('Tag','radio_align2d_mask_rect'),'Value',1);
else
    set(findobj('Tag','radio_align2d_mask_rect'),'Value',0);
end

set(findobj('Tag','input_align2d_mask_sphere_radius'),'String',num2str(storage_align2dgui.mask.(type).sphere.radius));
set(findobj('Tag','input_align2d_mask_sphere_sigma'),'String',num2str(storage_align2dgui.mask.(type).sphere.sigma));
set(findobj('Tag','input_align2d_mask_sphere_centerx'),'String',num2str(storage_align2dgui.mask.(type).sphere.center.x));
set(findobj('Tag','input_align2d_mask_sphere_centery'),'String',num2str(storage_align2dgui.mask.(type).sphere.center.y));

set(findobj('Tag','input_align2d_mask_rect_radiusx'),'String',num2str(storage_align2dgui.mask.(type).rect.radius.x));
set(findobj('Tag','input_align2d_mask_rect_radiusy'),'String',num2str(storage_align2dgui.mask.(type).rect.radius.y));
set(findobj('Tag','input_align2d_mask_rect_sigma'),'String',num2str(storage_align2dgui.mask.(type).rect.sigma));
set(findobj('Tag','input_align2d_mask_rect_centerx'),'String',num2str(storage_align2dgui.mask.(type).rect.center.x));
set(findobj('Tag','input_align2d_mask_rect_centery'),'String',num2str(storage_align2dgui.mask.(type).rect.center.y));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  load particle stack                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loadstack()

global storage_align2dgui;

filename = get(findobj('Tag','input_align2d_instack'),'String');
if ~isempty(filename)
    header = tom_reademheader(filename);
    storage_align2dgui.particlesize = header.Header.Size;
    storage_align2dgui.sampleimage = tom_emreadc(filename,'subregion',[1 1 1],[header.Header.Size(1)-1 header.Header.Size(2)-1 0]);

    string = ['size:' num2str(header.Header.Size(1)) ' x ' num2str(header.Header.Size(1))];
    string = strvcat(string,['# images in stack:' num2str(header.Header.Size(3))]);
    set(findobj('Tag','infotext_align2d_stack'),'String',string);

    tmpobj = findobj('Tag','im_align2d_sample');
    axes(tmpobj);
    tom_imagesc(storage_align2dgui.sampleimage,'noinfo');
    set(tmpobj,'Tag','im_align2d_sample');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  load particle ref                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loadref()

global storage_align2dgui;

filename = get(findobj('Tag','input_align2d_inref'),'String');
if ~isempty(filename)
    header = tom_reademheader(filename);
    storage_align2dgui.refparticlesize = header.Header.Size;
    storage_align2dgui.samplerefimage = tom_emreadc(filename,'subregion',[1 1 1],[header.Header.Size(1)-1 header.Header.Size(2)-1 0]);

    string = ['size:' num2str(header.Header.Size(1)) ' x ' num2str(header.Header.Size(1))];
    string = strvcat(string,['# images in reference:' num2str(header.Header.Size(3))]);
    set(findobj('Tag','infotext_align2d_ref'),'String',string);

    tmpobj = findobj('Tag','im_align2d_refsample');
    axes(tmpobj);
    tom_imagesc(storage_align2dgui.samplerefimage,'noinfo');
    set(tmpobj,'Tag','im_align2d_refsample');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Unused Callbacks/Create functions                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_align2d_instack_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_align2d_iterations_alignment_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_align2d_iterations_refinement_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_align2d_outref_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_align2d_inref_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_align2d_outalign_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_align2d_outstack_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_align2d_mask_rect_radiusy_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_align2d_mask_rect_radiusx_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_align2d_mask_rect_centery_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_align2d_mask_rect_centerx_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_align2d_mask_rect_sigma_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_align2d_mask_sphere_centery_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_align2d_mask_sphere_centerx_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_align2d_mask_sphere_radius_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_align2d_mask_sphere_sigma_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_align2d_filter_kernel_radius_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_align2d_filter_times_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_align2d_filter_bandpass_smooth_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_align2d_inalign_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_align2d_filter_bandpass_low_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_align2d_filter_bandpass_high_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_align2d_sub_stacks.
function checkbox_align2d_sub_stacks_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_align2d_sub_stacks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_align2d_sub_stacks


