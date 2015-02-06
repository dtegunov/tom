function varargout = tom_av3_createalign(varargin)
%TOM_AV3_CREATEALIGN creates ...
%
%   varargout = tom_av3_createalign(varargin)
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
%   ... = tom_av3_createalign(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by AK mm/dd/yy
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
                   'gui_OpeningFcn', @tom_av3_createalign_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_av3_createalign_OutputFcn, ...
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
function tom_av3_createalign_OpeningFcn(hObject, eventdata, handles, varargin)


setappdata(0,'UseNativeSystemDialogs',0);

% Choose default command line output for tom_av3_createalign
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes tom_av3_createalign wait for user response (see UIRESUME)
uiwait(handles.figure1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Output Function                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = tom_av3_createalign_OutputFcn(hObject, eventdata, handles) 

% Get default command line output from handles structure
varargout{1} = handles.output;
delete(handles.figure1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Close Function                                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function closefcn(hObject, eventdata, handles)

setappdata(0,'UseNativeSystemDialogs',1);
guidata(hObject, handles);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Add volumes                                                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_av3_aligncreate_add_Callback(hObject, eventdata, handles)

liststring = get(findobj('Tag','input_av3_aligncreate_volumes'),'String');

[filename, pathname] = uigetfile('*.*', 'Pick input volumes','MultiSelect','on');
if iscell(filename) || ischar(filename)
    for i = 1:size(filename,2)

        if strcmp(filename{i},'.') ~= 1 && strcmp(filename{i},'..') ~= 1
            liststring = strvcat(liststring,strcat(pathname,filename{i}));
        end
    end

    oldstring = get(findobj('Tag','input_av3_aligncreate_volumes'),'String');
    set(findobj('Tag','input_av3_aligncreate_volumes'),'String',liststring);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  remove volumes                                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_av3_aligncreate_remove_Callback(hObject, eventdata, handles)

newstring  = '';
vals = get(findobj('Tag','input_av3_aligncreate_volumes'),'Value');
string = get(findobj('Tag','input_av3_aligncreate_volumes'),'String');

for i=1:size(string,1)
    if (i == vals) == zeros(1,size(vals,2))
        newstring = strvcat(newstring,string(i,:));
    end
end

set(findobj('Tag','input_av3_aligncreate_volumes'),'String',newstring,'Value',1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  browse                                                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_av3_aligncreate_browse_Callback(hObject, eventdata, handles)

[filename, pathname] = uiputfile({'*.mat'}, 'Pick a filename');
if ischar(filename)
    set(findobj('Tag','input_av3_aligncreate_outfile'),'String',[pathname filename]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  go                                                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_av3_aligncreate_go_Callback(hObject, eventdata, handles)

outfile = get(findobj('Tag','input_av3_aligncreate_outfile'),'String');
vols = get(findobj('Tag','input_av3_aligncreate_volumes'),'String');

if isempty(outfile)
    msgbox('Select output alignment file first.');
    return;
end

if isempty(vols)
    msgbox('Select input files first.');
    return;
end

handles.output = outfile;
guidata(hObject, handles);

Align = tom_av3_create_alignlist(vols,[],1);
save(outfile,'Align');
uiresume;
closefcn(hObject, eventdata, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Unused Callbacks/Create functions                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_av3_aligncreate_outfile_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_av3_aligncreate_volumes_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_av3_aligncreate_outfile_Callback(hObject, eventdata, handles)
function input_av3_aligncreate_volumes_Callback(hObject, eventdata, handles)