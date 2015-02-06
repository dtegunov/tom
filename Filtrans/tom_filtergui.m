function varargout = tom_filtergui(varargin)
%TOM_FILTERGUI for creating and editing filter and mask structures
%
%   varargout = tom_filtergui(varargin)
%
%SYNTAX
%   out_struct = tom_parallelsettings(type,in_struct,val_struct)
%PARAMETERS
%
%  INPUT
%   varargin            ...
%  
%  OUTPUT
%   varargout           ...
%
% INPUT
% type:         'mask' or 'filter', determines which structure is used
% in_struct:    structure with possible values for output stucture
%               in_struct.name.types = {'val1' 'val2' ... 'valn'}                
%               for each filter/mask in the output structure add one of the
%               above lines to your in_struct, replace 'name' with the
%               desired name of the filter/mask. 'val1' to 'valn' can be
%               any of the following definitions:
%
%               filters:
%                   bandpass:   bandpass filter
%                   kernel:     kernel filter (real or fourier space)
%
%               masks:
%                   sphere:     2D sphere mask
%                   sphere3d:   3D sphere mask
%                   rectangle:  2D rectangle mask
%                   cylinder3d: 3D cylinder mask
%
% val_struct:   mask/filter structure to be used for default values (optional)
%
% OUTPUT
% out_struct:   output structur, that can be used for input in
%               tom_create_mask and tom_apply_filter
%
%EXAMPLE
%   ... = tom_filtergui(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   tom_create_mask, tom_apply_filter
%
%   created by AK 16/02/06
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
                   'gui_OpeningFcn', @tom_filtergui_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_filtergui_OutputFcn, ...Value s
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
function tom_filtergui_OpeningFcn(hObject, eventdata, handles, varargin)

global storage_filterdlg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define new filters/masks here
% rock funny Hack in chek fields

storage_filterdlg.filterst.('bandpass').vals = {'times' 'low' 'high' 'smooth'}; 
storage_filterdlg.filterst.('bandpass').valsoptional = {0 0 0 1};
storage_filterdlg.filterst.('bandpass').active_defaults = {0 0 0 0 0 0 0 0 0 0 0};
storage_filterdlg.filterst.('bandpass').dropdowns = struct();

storage_filterdlg.filterst.('kernel').vals = {'times' 'radius'}; 
storage_filterdlg.filterst.('kernel').valsoptional = {0 0};
storage_filterdlg.filterst.('kernel').active_defaults = {0 0 0 0 0 0 0 0 0 0 0};
storage_filterdlg.filterst.('kernel').dropdowns.space = {'real' 'fourier'};
storage_filterdlg.filterst.('kernel').dropdowns.method = {'quadr','circ'}; 


storage_filterdlg.filterst.('sphere').vals = {'size_x' 'size_y' 'radius' 'sigma' 'center_x' 'center_y','angle'}; 
storage_filterdlg.filterst.('sphere').valsoptional = {0 0 0 1 1 1 1};
storage_filterdlg.filterst.('sphere').active_defaults = {1 1 1 0 0 0 0};
storage_filterdlg.filterst.('sphere').dropdowns = struct();

storage_filterdlg.filterst.('sphere3d').vals = {'size_x' 'size_y' 'size_z' 'radius' 'sigma' 'center_x' 'center_y' 'center_z' 'angle_phi' 'angle_psi' 'angle_theta'}; 
storage_filterdlg.filterst.('sphere3d').valsoptional = {0 0 0 0 1 1 1 1 1 1 1};
storage_filterdlg.filterst.('sphere3d').active_defaults = {1 1 1 0 0 0 0 0 0 0 0};
storage_filterdlg.filterst.('sphere3d').dropdowns = struct();

storage_filterdlg.filterst.('rectangle').vals = {'size_x' 'size_y'  'masksize_x' 'masksize_y' 'sigma' 'center_x' 'center_y' 'angle_phi' 'angle_psi' 'angle_theta'}; 
storage_filterdlg.filterst.('rectangle').valsoptional = {0 0 0 0 1 1 1 1 1 1};
storage_filterdlg.filterst.('rectangle').active_defaults = {1 1 0 0 0 0 0 0 0 0 0};
storage_filterdlg.filterst.('rectangle').dropdowns = struct();

storage_filterdlg.filterst.('cylinder3d').vals = {'size_x' 'size_y' 'size_z'  'radius' 'sigma' 'center_x' 'center_y' 'center_z','angle_phi','angle_psi','angle_theta'}; 
storage_filterdlg.filterst.('cylinder3d').valsoptional = {0 0 0 1 1 1 1 1 1 1 1};
storage_filterdlg.filterst.('cylinder3d').active_defaults = {1 1 1 0 0 0 0 0 0 0 0 };
storage_filterdlg.filterst.('cylinder3d').dropdowns = struct();

storage_filterdlg.filterst.('ellipse3d').vals = {'size_x' 'size_y' 'size_z'  'radius1' 'radius2' 'radius3' 'sigma' 'center_x' 'center_y' 'center_z','angle_phi','angle_psi','angle_theta'}; 
storage_filterdlg.filterst.('ellipse3d').valsoptional = {0 0 0 1 1 1 1 1 1 1 1 1 1};
storage_filterdlg.filterst.('ellipse3d').active_defaults = {1 1 1 0 0 0 0 0 0 0 0 };
storage_filterdlg.filterst.('ellipse3d').dropdowns = struct();

storage_filterdlg.filterst.('ellipse').vals = {'size_x' 'size_y' 'radius1' 'radius2' 'sigma' 'center_x' 'center_y' 'angle_phi' 'angle_psi' 'angle_theta'}; 
storage_filterdlg.filterst.('ellipse').valsoptional = {0 0 1 1 1 1 1 1 1 1};
storage_filterdlg.filterst.('ellipse').active_defaults = {1 1 0 0 0 0 0 0 0 0 0};
storage_filterdlg.filterst.('ellipse').dropdowns = struct();



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 3
    error('no input structure given.');
end

storage_filterdlg.keystruct = varargin{2};

fnames=fieldnames(storage_filterdlg.keystruct);
for i=1:size(fnames,1)
    if i <= 10
        set(handles.(['radiobutton_filterdlg_name_' num2str(i)]),'Visible','on','String',fnames(i));
    else
        error('Max number of names supported: 10');
    end
end

% input values given
if nargin > 5
    storage_filterdlg.outst = varargin{3};
else
    %build output structure with defaults
    for i=1:size(fnames,1)
        storage_filterdlg.outst.(fnames{i}).Apply = 2;
        storage_filterdlg.outst.(fnames{i}).Value = struct();
        storage_filterdlg.outst.(fnames{i}).Type = '';
    end
end

storage_filterdlg.outst_orig = storage_filterdlg.outst;

dispfiltertype();

loadvals();

handles.output = hObject;
guidata(hObject, handles);
uiwait(handles.figure1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Output Function                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = tom_filtergui_OutputFcn(hObject, eventdata, handles) 

global storage_filterdlg;


varargout{1} = storage_filterdlg.outst;
storage_filterdlg=[];
delete(handles.figure1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Filter name                                                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radiobutton_filterdlg_name_1_Callback(hObject, eventdata, handles)

set_filtertype();

if savevals()
    dispfiltertype();
    loadvals();
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radiobutton_filterdlg_name_2_Callback(hObject, eventdata, handles)

set_filtertype();

if savevals()
    dispfiltertype();
    loadvals();
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radiobutton_filterdlg_name_3_Callback(hObject, eventdata, handles)

set_filtertype();

if savevals()
    dispfiltertype();
    loadvals();
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radiobutton_filterdlg_name_4_Callback(hObject, eventdata, handles)

set_filtertype();

if savevals()
    dispfiltertype();
    loadvals();
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radiobutton_filterdlg_name_5_Callback(hObject, eventdata, handles)

set_filtertype();

if savevals()
    dispfiltertype();
    loadvals();
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radiobutton_filterdlg_name_6_Callback(hObject, eventdata, handles)

set_filtertype();

if savevals()
    dispfiltertype();
    loadvals();
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radiobutton_filterdlg_name_7_Callback(hObject, eventdata, handles)

set_filtertype();

if savevals()
    dispfiltertype();
    loadvals();
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radiobutton_filterdlg_name_8_Callback(hObject, eventdata, handles)

set_filtertype();

if savevals()
    dispfiltertype();
    loadvals();
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radiobutton_filterdlg_name_9_Callback(hObject, eventdata, handles)

set_filtertype();

if savevals()
    dispfiltertype();
    loadvals();
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radiobutton_filterdlg_name_10_Callback(hObject, eventdata, handles)

set_filtertype();

if savevals()
    dispfiltertype();
    loadvals();
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Filter types                                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radiobutton_filterdlg_type_1_Callback(hObject, eventdata, handles)
dispfiltervals();

function radiobutton_filterdlg_type_2_Callback(hObject, eventdata, handles)
dispfiltervals();

function radiobutton_filterdlg_type_3_Callback(hObject, eventdata, handles)
dispfiltervals();

function radiobutton_filterdlg_type_4_Callback(hObject, eventdata, handles)
dispfiltervals();

function radiobutton_filterdlg_type_5_Callback(hObject, eventdata, handles)
dispfiltervals();

function radiobutton_filterdlg_type_6_Callback(hObject, eventdata, handles)
dispfiltervals();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Exit                                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_filterdlg_exit_Callback(hObject, eventdata, handles)

if savevals()
    uiresume(handles.figure1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Reset                                                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_filterdlg_reset_Callback(hObject, eventdata, handles)

global storage_filterdlg;

storage_filterdlg.outst = storage_filterdlg.outst_orig;

dispfiltertype();
loadvals();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Checkbox filter defaults                                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkbox_filterdlg_defaults_Callback(hObject, eventdata, handles)

if get(hObject,'Value') == 1
    set(findobj('Tag','checkbox_filterdlg_off'),'Value',0);
end

checkfields();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Checkbox filter off                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkbox_filterdlg_off_Callback(hObject, eventdata, handles)

if get(hObject,'Value') == 1
    set(findobj('Tag','checkbox_filterdlg_defaults'),'Value',0);
end

checkfields();


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
%%  displays a filter's types                                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dispfiltertype()

global storage_filterdlg;

%make all filter types invisible
set(findobj('-regexp','Tag','radiobutton_filterdlg_type_*'),'Visible','off');

%find out which filter name is selected
for i=1:10
    if get(findobj('Tag',['radiobutton_filterdlg_name_' num2str(i)]),'Value') == 1
        selectedfilter = cell2mat(get(findobj('Tag',['radiobutton_filterdlg_name_' num2str(i)]),'String'));
        storage_filterdlg.oldselectedfilter = selectedfilter;
        break;
    end
end

if isempty(storage_filterdlg.outst.(selectedfilter).Type)
    set(findobj('Tag','radiobutton_filterdlg_type_1'),'Value',1);
end

%display all filter types
for i=1:size(storage_filterdlg.keystruct.(selectedfilter).types,2)
    set(findobj('Tag',['radiobutton_filterdlg_type_' num2str(i)]),'Visible','on','String',storage_filterdlg.keystruct.(selectedfilter).types(i));
end

dispfiltervals();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  displays a filter's values                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dispfiltervals()

global storage_filterdlg;

%make all filter values invisible
set(findobj('-regexp','Tag','edit_filterdlg_*'),'Visible','off','String','');
set(findobj('-regexp','Tag','text_filterdlg_*'),'Visible','off','String','');
set(findobj('-regexp','Tag','popupmenu_filterdlg_*'),'Visible','off');


%find out which filter type is selected
for i=1:6
    if get(findobj('Tag',['radiobutton_filterdlg_type_' num2str(i)]),'Value') == 1
        selectedfiltertype = cell2mat(get(findobj('Tag',['radiobutton_filterdlg_type_' num2str(i)]),'String'));
        storage_filterdlg.oldselectedfiltertype = selectedfiltertype;
        break;
    end
end


%display dropdown fields
if ~isempty(storage_filterdlg.filterst.(selectedfiltertype).dropdowns)
    fnames = fieldnames(storage_filterdlg.filterst.(selectedfiltertype).dropdowns);
    for i=1:size(fnames,1)
        if i<=4
            string = '';
            for j=1:size(storage_filterdlg.filterst.(selectedfiltertype).dropdowns.(fnames{i}),2)
                string = strvcat(string, storage_filterdlg.filterst.(selectedfiltertype).dropdowns.(fnames{i}){j});
            end
            set(findobj('Tag',['text_filterdlg_dropdown_' num2str(i)]),'Visible','on','String',fnames{i}); 
            set(findobj('Tag',['popupmenu_filterdlg_' num2str(i)]),'Visible','on','String',string);
        else
            error('Max number of dropdowns supported: 4');
        end

    end
end

%display input fields
for i=1:size(storage_filterdlg.filterst.(selectedfiltertype).vals,2)
    if i<=11
        if storage_filterdlg.filterst.(selectedfiltertype).valsoptional{i} == 1
            set(findobj('Tag',['text_filterdlg_' num2str(i)]),'Visible','on','String',[storage_filterdlg.filterst.(selectedfiltertype).vals{i} ' (optional)']); 
        else
            set(findobj('Tag',['text_filterdlg_' num2str(i)]),'Visible','on','String',storage_filterdlg.filterst.(selectedfiltertype).vals{i}); 
        end
        set(findobj('Tag',['edit_filterdlg_' num2str(i)]),'Visible','on','String','');
    else
        error('Max number of input fields supported: 11');
    end
end
checkfields();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  loads filter's values from internal structure                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loadvals()

global storage_filterdlg;

%find out which filter name is selected
for i=1:10
    if get(findobj('Tag',['radiobutton_filterdlg_name_' num2str(i)]),'Value') == 1
        selectedfilter = cell2mat(get(findobj('Tag',['radiobutton_filterdlg_name_' num2str(i)]),'String'));
        break;
    end
end

%find out which filter type is selected
for i=1:6
    if get(findobj('Tag',['radiobutton_filterdlg_type_' num2str(i)]),'Value') == 1
        selectedfiltertype = cell2mat(get(findobj('Tag',['radiobutton_filterdlg_type_' num2str(i)]),'String'));
        break;
    end
end

%filter off
if storage_filterdlg.outst.(selectedfilter).Apply == 0
    set(findobj('Tag','checkbox_filterdlg_off'),'Value',1);
    set(findobj('Tag','checkbox_filterdlg_defaults'),'Value',0);
%filter defaults
elseif storage_filterdlg.outst.(selectedfilter).Apply == 2
    set(findobj('Tag','checkbox_filterdlg_defaults'),'Value',1);
    set(findobj('Tag','checkbox_filterdlg_off'),'Value',0);
else
%load values
    set(findobj('Tag','checkbox_filterdlg_defaults'),'Value',0);
    set(findobj('Tag','checkbox_filterdlg_off'),'Value',0);
    
    %update dropdown fields
    for field=fieldnames(storage_filterdlg.filterst.(selectedfiltertype).dropdowns)'
        tag = get(findobj('String',field{1}),'Tag');
        contents = get(findobj('Tag',['popupmenu_filterdlg_' tag(end)]),'String');
        val = strmatch(storage_filterdlg.outst.(selectedfilter).Value.(field{1}),contents,'exact');
        if ~isempty(val)
            set(findobj('Tag',['popupmenu_filterdlg_' tag(end)]),'Value',val);
        else
            set(findobj('Tag',['popupmenu_filterdlg_' tag(end)]),'Value',1);
            disp('WARNING: Method field empty, inserting default');
        end
    end
        
    %update input fields
    for i=1:size(storage_filterdlg.filterst.(selectedfiltertype).vals,2)
        set(findobj('Tag',['edit_filterdlg_' num2str(i)]),'String',num2str(storage_filterdlg.outst.(selectedfilter).Value.(storage_filterdlg.filterst.(selectedfiltertype).vals{i})));
    end
    
end

checkfields();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  saves filter's values to internal structure                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function status = savevals()

global storage_filterdlg;

status = 1;

%find out which filter type was selected
selectedfiltertype = storage_filterdlg.oldselectedfiltertype;
selectedfilter = storage_filterdlg.oldselectedfilter;

defaultsflag = get(findobj('Tag','checkbox_filterdlg_defaults'),'Value');
offflag = get(findobj('Tag','checkbox_filterdlg_off'),'Value');

if defaultsflag == 0 & offflag == 0
    %check input of all fields
    for i=1:size(storage_filterdlg.filterst.(selectedfiltertype).vals,2)
        val = str2num(get(findobj('Tag',['edit_filterdlg_' num2str(i)]),'String'));
        if isempty(val) & storage_filterdlg.filterst.(selectedfiltertype).valsoptional{i} == 0
            errordlg(['value ''' cell2mat(storage_filterdlg.filterst.(selectedfiltertype).vals(i)) ''' is required!']);
            status = 0;
            %reset filter type selection
            set(findobj('-regexp','Tag','radiobutton_filterdlg_type_*'),'Value',0);
            set(findobj('String',{selectedfiltertype}),'Value',1);
            %reset filter name selection
            set(findobj('-regexp','Tag','radiobutton_filterdlg_name_*'),'Value',0);
            set(findobj('String',{selectedfilter}),'Value',1);

            break;
        end
    end
end

%no error in the input fields, save the values
if status == 1
    
    if defaultsflag == 1
        applyval = 2;
    elseif offflag == 1
        applyval = 0;
    else
        applyval = 1;
    end
   
    storage_filterdlg.outst.(selectedfilter).Apply = applyval;
    storage_filterdlg.outst.(selectedfilter).Type = selectedfiltertype;

    %get values of dropdown fields
    for field=fieldnames(storage_filterdlg.filterst.(selectedfiltertype).dropdowns)'
       tag = get(findobj('String',field{1}),'Tag');
       contents = get(findobj('Tag',['popupmenu_filterdlg_' tag(end)]),'String'); 
       storage_filterdlg.outst.(selectedfilter).Value.(field{1}) = strtrim(contents(get(findobj('Tag',['popupmenu_filterdlg_' tag(end)]),'Value'),:));
    end
    %get values of input fields
    for i=1:size(storage_filterdlg.filterst.(selectedfiltertype).vals,2)
        storage_filterdlg.outst.(selectedfilter).Value.(storage_filterdlg.filterst.(selectedfiltertype).vals{i}) = str2num(get(findobj('Tag',['edit_filterdlg_' num2str(i)]),'String'));
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Checks if fields should be enabled or disabled                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkfields()

global storage_filterdlg;

%filter enable
if get(findobj('Tag','checkbox_filterdlg_off'),'Value') == 1 | get(findobj('Tag','checkbox_filterdlg_defaults'),'Value') == 1
    flag = 'off';
%filter disable
else
    flag = 'on';
end

selectedfiltertype = storage_filterdlg.oldselectedfiltertype;

start=1;

%input fields ...HACK FB
if (strcmp(selectedfiltertype,'sphere3d') & strcmp(flag,'off'))
    start=4;
end;

if (strcmp(selectedfiltertype,'sphere') & strcmp(flag,'off'))
    start=3;
end;


if (strcmp(selectedfiltertype,'rectangle') & strcmp(flag,'off'))
    start=3;
end;

if (strcmp(flag,'off'))
    flag_o='on';
else
    flag_o='off';
end;

 def_on=get(findobj('Tag','checkbox_filterdlg_defaults'),'Value');

for i=1:11
    set(findobj('Tag',['text_filterdlg_' num2str(i)]),'Enable',flag_o); 
    set(findobj('Tag',['edit_filterdlg_' num2str(i)]),'Enable',flag_o);
end


for i=1:11
    try
        if (eval(['storage_filterdlg.filterst.' selectedfiltertype '.active_defaults{i}'])==0 |def_on==0 )
            set(findobj('Tag',['text_filterdlg_' num2str(i)]),'Enable',flag);
            set(findobj('Tag',['edit_filterdlg_' num2str(i)]),'Enable',flag);
        end;
    catch
    end;
end

for i=1:4
    set(findobj('Tag',['text_filterdlg_dropdown_' num2str(i)]),'Enable',flag);
    set(findobj('Tag',['popupmenu_filterdlg_' num2str(i)]),'Enable',flag);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Set the filter type                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function set_filtertype()

global storage_filterdlg;

%find out which filter name is selected

for i=1:10
    if get(findobj('Tag',['radiobutton_filterdlg_name_' num2str(i)]),'Value') == 1
        selectedfilter = cell2mat(get(findobj('Tag',['radiobutton_filterdlg_name_' num2str(i)]),'String'));
        break;
    end
end

%set filter type
filtertype = storage_filterdlg.outst.(selectedfilter).Type;
set(findobj('String',{filtertype}),'Value',1);

if isempty(filtertype)
    set(findobj('Tag','radiobutton_filterdlg_type_1'),'Value',1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Unused Callbacks/Create functions                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_filterdlg_8_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function popupmenu_filterdlg_1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_filterdlg_7_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_filterdlg_6_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_filterdlg_5_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_filterdlg_4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_filterdlg_3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_filterdlg_2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_filterdlg_1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function popupmenu_filterdlg_4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function popupmenu_filterdlg_3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function popupmenu_filterdlg_2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function popupmenu_filterdlg_1_Callback(hObject, eventdata, handles)
function popupmenu_filterdlg_2_Callback(hObject, eventdata, handles)
function popupmenu_filterdlg_3_Callback(hObject, eventdata, handles)
function popupmenu_filterdlg_4_Callback(hObject, eventdata, handles)
function edit_filterdlg_1_Callback(hObject, eventdata, handles)
function edit_filterdlg_2_Callback(hObject, eventdata, handles)
function edit_filterdlg_3_Callback(hObject, eventdata, handles)
function edit_filterdlg_4_Callback(hObject, eventdata, handles)
function edit_filterdlg_5_Callback(hObject, eventdata, handles)
function edit_filterdlg_6_Callback(hObject, eventdata, handles)
function edit_filterdlg_7_Callback(hObject, eventdata, handles)
function edit_filterdlg_8_Callback(hObject, eventdata, handles)
function edit_filterdlg_9_Callback(hObject, eventdata, handles)
function edit_filterdlg_9_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_filterdlg_10_Callback(hObject, eventdata, handles)
function edit_filterdlg_10_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_filterdlg_11_Callback(hObject, eventdata, handles)
function edit_filterdlg_11_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


