function varargout = tom_av3_extractanglesshiftsgui(varargin)
%TOM_AV3_EXTRACTANGLESSHIFTSGUI is a GUI for tom_av3_extractanglesshifts
%
%   varargout = tom_av3_extractanglesshiftsgui(varargin)
%
%This is a frontend for tom_av3_extractanglesshifts
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
%   ... = tom_av3_extractanglesshiftsgui(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   TOM_AV3_OSCARGUI
%
%   created by AK 10/10/05
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
                   'gui_OpeningFcn', @tom_av3_extractanglesshiftsgui_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_av3_extractanglesshiftsgui_OutputFcn, ...
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
function tom_av3_extractanglesshiftsgui_OpeningFcn(hObject, eventdata, handles, varargin)

global storage_extract;

if size(varargin,2) == 3
    storage_extract.oscarvolume = varargin{1};
    set(findobj('Tag','input_oscarvolume'), 'String', storage_extract.oscarvolume);
    storage_extract.alignfile = varargin{2};
    set(findobj('Tag','input_alignfile'), 'String', storage_extract.alignfile);
    storage_extract.maskfile = '';
%    storage_extract.maskfile = varargin{3};
%    set(findobj('Tag','input_maskfile'), 'String', storage_extract.maskfile);
    display_parameters();
else
    storage_extract.oscarvolume = '';
    storage_extract.alignfile = '';
    storage_extract.maskfile = '';
end

%setappdata(0,'UseNativeSystemDialogs',0);
handles.output = hObject;
guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Output Function                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = tom_av3_extractanglesshiftsgui_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Close Function                                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function closefcn(hObject, eventdata, handles)

storage_extract = [];
delete(findobj('Tag','tom_av3_extractanglesshiftsgui'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  browse for oscar volume                                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_browseoscarvolume_Callback(hObject, eventdata, handles)

global storage_extract;

[filename, pathname] = uigetfile({'*.ccf'}, 'Pick an oscar output');
if ischar(filename)
    storage_extract.oscarvolume = [pathname filename];
    storage_extract.oscarvolume = storage_extract.oscarvolume(1:max(strfind(storage_extract.oscarvolume,'.'))-1);
    set(findobj('Tag','input_oscarvolume'),'String',storage_extract.oscarvolume);
    display_parameters();
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  directly edit oscar volume                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_oscarvolume_Callback(hObject, eventdata, handles)

global storage_extract;

storage_extract.oscarvolume = get(hObject, 'String');
display_parameters();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  browse for align file                                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_browse_alignfile_Callback(hObject, eventdata, handles)

global storage_extract;

[filename, pathname] = uigetfile({'*.mat'}, 'Pick an align file');
if ischar(filename)
    set(findobj('Tag','input_alignfile'),'String',[pathname filename]);
    storage_extract.alignfile = [pathname filename];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  directly edit align file                                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_alignfile_Callback(hObject, eventdata, handles)

global storage_extract;

storage_extract.alignfile = get(hObject, 'String');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  browse for mask file                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_browsemaskfile_Callback(hObject, eventdata, handles)

global storage_extract;

[filename, pathname] = uigetfile({'*.vol','*.em'}, 'Pick a mask file');
if ischar(filename)
    set(findobj('Tag','input_maskfile'),'String',[pathname filename]);
    storage_extract.maskfile = [pathname filename];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  directly edit mask file                                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_maskfile_Callback(hObject, eventdata, handles)

global storage_extract;

storage_extract.maskfile = get(hObject, 'String');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  generate mask file                                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_maskgen_Callback(hObject, eventdata, handles)

global storage_extract;

if isempty(storage_extract.alignfile)
    errordlg('Select particle align list first!','Error');
    return;
end

a = load(storage_extract.alignfile);
tsize = a.Align(1,1).Tomogram.Header.Size';
prompt = {'Enter radius in pixel:','Enter smoothing in pixel (outside of radius)'};
dlg_title = 'Generate mask';
num_lines = 1;
def = {num2str(round(tsize(1)./2)),'1'};
answer = inputdlg(prompt,dlg_title,num_lines,def,'on');

if ~isempty(answer)
    [filename, pathname] = uiputfile({'*.vol','*.em'}, 'Save mask as');
    if ischar(filename)
        mask=ones(tsize);
        mask=tom_spheremask(mask,str2num(answer{1}),str2num(answer{2}),[tsize(1)./2+1,tsize(2)./2+1,tsize(3)./2+1]);
        
        storage_oscar.mask.filename = [pathname,filename];
        tom_emwrite(storage_oscar.mask.filename, mask);
        set(findobj('Tag','input_maskfile'),'String',storage_oscar.mask.filename);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  extract angles and shifts and fill out align file                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_extract_Callback(hObject, eventdata, handles)

global storage_extract;

%extract the values and update alignment file
mask = tom_emreadc(storage_extract.maskfile);
a = load(storage_extract.alignfile);
Align = a.Align;
Align = tom_av3_extract_anglesshifts(storage_extract.oscarvolume, Align, mask.Value, 1);
save(storage_extract.alignfile,'Align');

%Create average volume
if get(findobj('Tag','checkbox_createavg'),'Value') == 1
    al = tom_av3_align_sum(Align);
    average = tom_av3_average(al,'sum',0,0,1);
    average = tom_emheader(average);
    tempheader = tom_reademheader(al(1,1).Filename);
    average.Header.Objectpixelsize = tempheader.Header.Objectpixelsize;
    tom_emwrite(average);
    tom_volxyz(average.Value);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display oscar parameters                                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function display_parameters()

global storage_extract;

header = tom_reademheader(strcat(storage_extract.oscarvolume,'.ccf.norm'));
vsize = header.Header.Size;

%Extract angles from Comment
angs = sscanf(char(header.Header.Comment'),'Angular range: %f %f %f %f %f %f %f %f %f');
angle_start = [angs(1) angs(4) angs(7)];
angle_end = [angs(2) angs(5) angs(8)];
angular_incr = [angs(3) angs(6) angs(9)];
set(findobj('Tag','oscar_parameters'),'String',strvcat(['Phi: from ' num2str(angle_start(1)) ' to ' num2str(angle_end(1)) ' in steps of ' num2str(angular_incr(1)) ' degrees'],['Psi: from ' num2str(angle_start(2)) ' to ' num2str(angle_end(2)) ' in steps of ' num2str(angular_incr(2)) ' degrees'],['Theta: from ' num2str(angle_start(3)) ' to ' num2str(angle_end(3)) ' in steps of ' num2str(angular_incr(3)) ' degrees']));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Unused Callbacks/Create functions                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_oscarvolume_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_alignfile_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_maskfile_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function checkbox_createavg_Callback(hObject, eventdata, handles)
