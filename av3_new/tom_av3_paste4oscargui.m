function varargout = tom_av3_paste4oscargui(varargin)
%TOM_AV3_PASTE4OSCARGUI is a GUI for tom_av3_paste4oscar
%
%   varargout = tom_av3_paste4oscargui(varargin)
%
%This is a frontend for tom_av3_paste4oscar
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
%   ... = tom_av3_paste4oscargui(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   TOM_AV3_OSCARGUI
%
%   created by AK 07/10/05
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
                   'gui_OpeningFcn', @tom_av3_paste4oscargui_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_av3_paste4oscargui_OutputFcn, ...
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
function tom_av3_paste4oscargui_OpeningFcn(hObject, eventdata, handles, varargin)

global storage_paste;

storage_paste.filestring = [];
storage_paste.maskstring = [];
storage_paste.outvolstring = [];
storage_paste.outalignlist = [];
storage_paste.inalignlist = [];
storage_paste.vol.x = 0;
storage_paste.vol.y = 0;
storage_paste.vol.z = 0;
storage_paste.number.x = 0;
storage_paste.number.y = 0;
storage_paste.number.z = 0;
storage_paste.border.x = 0;
storage_paste.border.y = 0;
storage_paste.border.z = 0;
storage_paste.tdim = [];
storage_paste.normflag = 0;
setappdata(0,'UseNativeSystemDialogs',0);

handles.output = hObject;
guidata(hObject, handles);
uiwait;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Output Function                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = tom_av3_paste4oscargui_OutputFcn(hObject, eventdata, handles) 

global storage_paste;
%varargout{1} = handles.output;
varargout{1} = storage_paste.outvolstring;
varargout{2} = storage_paste.outalignlist;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Close Function                                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function closefcn(hObject, eventdata, handles)

storage_paste = [];
setappdata(0,'UseNativeSystemDialogs',1);
delete(findobj('Tag','tom_av3_paste4oscargui'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Browse for input align list                                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_browseinputalignlist_Callback(hObject, eventdata, handles)

global storage_paste;

[filename, pathname] = uigetfile({'*.mat'}, 'Pick an alignlist');
if ischar(filename)
    set(findobj('Tag','input_alignlist'),'String',[pathname filename]);
    storage_paste.inalignlist = [pathname filename];

    filestring  = '';
    try
        a = load(storage_paste.inalignlist);
        Align = a.Align;
        for i = 1:size(Align,2)
            filestring = strvcat(filestring,Align(1,i).Filename);
        end
        storage_paste.filestring = filestring;
        set(findobj('Tag','input_listbox'),'String',filestring);
        status = calculate_properties();
        
        if Align(size(Align,1),1).NormFlag == 1
            set(findobj('Tag','checkbox_normflag'),'Value',1);
        else
            set(findobj('Tag','checkbox_normflag'),'Value',0);
        end
        
        if Align(size(Align,1),1).Filter(2) ~= 0
            set(findobj('Tag','checkbox_bandpass'),'Value',1);
            set(findobj('Tag','bandpass_low'),'String',num2str(Align(size(Align,1),1).Filter(1)));
            set(findobj('Tag','bandpass_high'),'String',num2str(Align(size(Align,1),1).Filter(2)));
        else
            set(findobj('Tag','checkbox_bandpass'),'Value',0);
        end
        
    catch
        errordlg('Not a valid align list!','File error');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  directly edit input align list                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_alignlist_Callback(hObject, eventdata, handles)

global storage_paste;

filename = get(hObject, 'String');
storage_paste.inalignlist = filename;
filestring  = '';
try
    a = load(storage_paste.inalignlist);
    Align = a.Align;
    for i = 1:size(Align,2)
        filestring = strvcat(filestring,Align(1,i).Filename);
    end
    storage_paste.filestring = filestring;
    set(findobj('Tag','input_listbox'),'String',filestring);
    status = calculate_properties();
catch
    errordlg('Not a valid align list!','File error');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Add input volumes to listbox                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_addlist_Callback(hObject, eventdata, handles)

global storage_paste;

liststring = get(findobj('Tag','input_listbox'),'String');

[filename, pathname] = uigetfile('*.*', 'Pick input volumes','MultiSelect','on');
if iscell(filename) | ischar(filename)
    for i = 1:size(filename,2)

        if strcmp(filename{i},'.') ~= 1 & strcmp(filename{i},'..') ~= 1
            liststring = strvcat(liststring,strcat(pathname,filename{i}));
        end
    end

    oldstring = get(findobj('Tag','input_listbox'),'String');
    set(findobj('Tag','input_listbox'),'String',liststring);
    storage_paste.filestring = liststring;

    status = calculate_properties();
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Remove input volumes from listbox                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_removelist_Callback(hObject, eventdata, handles)

global storage_paste;

newstring  = '';
vals = get(findobj('Tag','input_listbox'),'Value');
string = get(findobj('Tag','input_listbox'),'String');

for i=1:size(string,1)
    if (i == vals) == zeros(1,size(vals,2))
        newstring = strvcat(newstring,string(i,:));
    end
end

set(findobj('Tag','input_listbox'),'String',newstring,'Value',1);
storage_paste.filestring = newstring;

status = calculate_properties();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Select mask                                                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_browsem_Callback(hObject, eventdata, handles)

global storage_paste;

[filename, pathname] = uigetfile({'*.vol','*.em'}, 'Pick a filename');
if ischar(filename)
    set(findobj('Tag','maskfield'),'String',[pathname filename]);
    storage_paste.maskstring = [pathname filename];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  directly edit mask                                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function maskfield_Callback(hObject, eventdata, handles)

global storage_paste;

filename = get(hObject, 'String');
storage_paste.maskstring = filename;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Generate mask                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_maskgen_Callback(hObject, eventdata, handles)

global storage_paste;

if isempty(storage_paste.filestring)
    errordlg('Select subvolumes first!','Error');
    return;
end

prompt = {'Enter radius in pixel:','Enter smoothing in pixel (outside of radius)'};
dlg_title = 'Generate mask';
num_lines = 1;
def = {num2str(storage_paste.tdim(1)./2),'1'};
answer = inputdlg(prompt,dlg_title,num_lines,def,'on');

if ~isempty(answer)
    [filename, pathname] = uiputfile({'*.vol','*.em'}, 'Save mask as');
    if ischar(filename)

        mask=ones(storage_paste.tdim');
        mask=tom_spheremask(mask,str2num(answer{1}),str2num(answer{2}),[storage_paste.tdim(1)/2+1,storage_paste.tdim(2)./2+1,storage_paste.tdim(3)./2+1]);
        
        storage_paste.maskstring = [pathname,filename];
        tom_emwrite(storage_paste.maskstring, mask);
        set(findobj('Tag','maskfield'),'String',storage_paste.maskstring);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Select output volume                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_browseoutputvolume_Callback(hObject, eventdata, handles)

global storage_paste;

[filename, pathname] = uiputfile({'*.vol','*.em'}, 'Pick a filename');
if ischar(filename)
    set(findobj('Tag','output_volume'),'String',[pathname filename]);
    storage_paste.outvolstring = [pathname filename];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  directly edit output volume                                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output_volume_Callback(hObject, eventdata, handles)

global storage_paste;

filename = get(hObject, 'String');
storage_paste.outvolstring = filename;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  checkbox normalize                                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkbox_normflag_Callback(hObject, eventdata, handles)

global storage_paste;

storage_paste.normflag = get(hObject,'Value');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Select output alignlist                                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_alignlistcreatebrowse_Callback(hObject, eventdata, handles)

global storage_paste;

[filename, pathname] = uiputfile({'*.mat'}, 'Pick a filename');
if ischar(filename)
    if isempty(strfind(filename,'.mat')) == 1
        filename = [filename '.mat'];
    end
    set(findobj('Tag','output_alignlist'),'String',[pathname filename]);
    storage_paste.outalignlist = [pathname filename];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  directly edit out align list                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output_alignlist_Callback(hObject, eventdata, handles)

global storage_paste;

filename = get(hObject, 'String');
    if isempty(strfind(filename,'.mat')) == 1
        filename = [filename '.mat'];
    end
storage_paste.outalignlist = filename;
set(hObject,'String',filename);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Start the paste process                                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_startpaste_Callback(hObject, eventdata, handles)

global storage_paste;

if isempty(storage_paste.outvolstring)
    errordlg('Please enter a output volume file name!','File error');
else

    if isempty(storage_paste.inalignlist)
        %Create alignlist
        Align = tom_av3_create_alignlist(storage_paste.filestring,[],1);
    else
        %Append alignlist
        a = load(storage_paste.inalignlist);
        Align = tom_av3_create_alignlist(storage_paste.filestring,a.Align,1);
    end
        
    %Filter
    if get(findobj('Tag','checkbox_bandpass'),'Value') == 1
        filter = [str2num(get(findobj('Tag','bandpass_low'),'String')),str2num(get(findobj('Tag','bandpass_high'),'String'))];
    else
        filter = [0 0];
    end
    
    Align = tom_av3_paste4oscar(storage_paste.filestring, storage_paste.outvolstring,storage_paste.maskstring, Align, filter, storage_paste.normflag, 1);
    save(storage_paste.outalignlist, 'Align');
    uiresume;
    closefcn(hObject, eventdata, handles);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Calculate dimensions of new volume and do some checks              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function retval = calculate_properties()

global storage_paste;

retval = 1;

%Check if at least one subvolume has been selected
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(storage_paste.filestring)
    retval = 0;
    set(findobj('Tag','button_startpaste'),'Enable','off');
    return;
else
    set(findobj('Tag','button_startpaste'),'Enable','on');
end

%Check if all files have the same dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:size(storage_paste.filestring,1)
    header = tom_reademheader(strtrim(storage_paste.filestring(i,:)));
    
    %Use the dimensions of the first template to test the others
    if i == 1
        storage_paste.tdim = header.Header.Size;
    else
        if (header.Header.Size == storage_paste.tdim) ~= ones(3,1)
            errordlg(['Dimensions of template ' strtrim(storage_paste.filestring(i,:)) ' differ from template ' strtrim(storage_paste.filestring(i,:)) '!'],'Dimensions mismatch');
            retval = 0;
            set(findobj('Tag','button_startpaste'),'Enable','off');
            return;
        end
    end
end

%Calculate the dimensions of the new template
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
storage_paste.number.x = floor(size(storage_paste.filestring,1)^(1./3));
storage_paste.number.y = storage_paste.number.x;
storage_paste.size.x = storage_paste.number.x*storage_paste.tdim(1)+storage_paste.tdim(1);
storage_paste.size.y = storage_paste.number.x*storage_paste.tdim(2)+storage_paste.tdim(2);
storage_paste.number.z = ceil(size(storage_paste.filestring,1)./storage_paste.number.x^2);
storage_paste.size.z = storage_paste.number.z*storage_paste.tdim(3)+storage_paste.tdim(3);

storage_paste.border.x = floor(storage_paste.tdim(1)./2);
storage_paste.border.y = floor(storage_paste.tdim(2)./2);
storage_paste.border.z = floor(storage_paste.tdim(3)./2);

set(findobj('Tag','bandpass_low'),'String','0');
set(findobj('Tag','bandpass_high'),'String',num2str(storage_paste.tdim(1)./2));

%Display the properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
emptyvols = storage_paste.number.x.*storage_paste.number.y.*storage_paste.number.z - size(storage_paste.filestring,1);
propstring = strvcat(['Dimensions of new volume: ' num2str(storage_paste.size.x) ' x ' num2str(storage_paste.size.y) ' x ' num2str(storage_paste.size.z)],['Number of subvolumes: ', num2str(storage_paste.number.x) ' x ' num2str(storage_paste.number.y) ' x ' num2str(storage_paste.number.z)],['Empty volumes: ' num2str(emptyvols)],['Particle dimensions: ' num2str(storage_paste.tdim(1)) ' x ' num2str(storage_paste.tdim(2)) ' x ' num2str(storage_paste.tdim(3))]);
set(findobj('Tag','outvol_properties'),'String',propstring);

if retval == 1
    set(findobj('Tag','button_startpaste'),'Enable','on');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Unused Callbacks/Create functions                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_listbox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_listbox_Callback(hObject, eventdata, handles)
function output_volume_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function output_alignlist_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_alignlist_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function maskfield_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function checkbox_bandpass_Callback(hObject, eventdata, handles)
function bandpass_low_Callback(hObject, eventdata, handles)
function bandpass_low_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function bandpass_high_Callback(hObject, eventdata, handles)
function bandpass_high_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end