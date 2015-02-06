function varargout = tom_av2_mergealign(varargin)
%GUI for merging 2D particle alignment files and creating stacks from the
%alignment structure
%
%SYNTAX
%tom_av2_mergealign
%
%
%SEE ALSO
%
%
%Copyright (c) 2005
%TOM toolbox for Electron Tomography
%Max-Planck-Institute for Biochemistry
%Dept. Molecular Structural Biology
%82152 Martinsried, Germany
%http://www.biochem.mpg.de/tom
%
%Created: 14/02/06 AK

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_av2_mergealign_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_av2_mergealign_OutputFcn, ...
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
function tom_av2_mergealign_OpeningFcn(hObject, eventdata, handles, varargin)

global storage_alignmerge;


set(findobj('Tag','listbox_inputclasses'),'Enable','off');
set(findobj('Tag','listbox_outputclasses'),'Enable','off');
set(findobj('Tag','button_merge_move'),'Enable','off');
set(findobj('Tag','button_merge_merge'),'Enable','off');

storage_alignmerge.aligns = '';
storage_alignmerge.liststring = '';
storage_alignmerge.inalign = '';
storage_alignmerge.outalign = '';
storage_alignmerge.classes = {};
storage_alignmerge.colors = {};
storage_alignmerge.radii = {};

setappdata(0,'UseNativeSystemDialogs',0);

handles.output = hObject;
guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Output Function                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = tom_av2_mergealign_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Add alignment file                                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_merge_addalign_Callback(hObject, eventdata, handles)

global storage_alignmerge;

liststring = get(findobj('Tag','listbox_loadaligns'),'String');

[filename, pathname] = uigetfile('*.mat', 'Pick input alignment files','MultiSelect','on');
if iscell(filename) | ischar(filename)
    for i = 1:size(filename,2)

        if strcmp(filename{i},'.') ~= 1 & strcmp(filename{i},'..') ~= 1
            liststring = strvcat(liststring,strcat(pathname,filename{i}));
        end
    end

    oldstring = get(findobj('Tag','listbox_loadaligns'),'String');
    set(findobj('Tag','listbox_loadaligns'),'String',liststring);
    storage_alignmerge.filestring = liststring;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Delete alignment file                                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_merge_deletealign_Callback(hObject, eventdata, handles)

global storage_alignmerge;

newstring  = '';
vals = get(findobj('Tag','listbox_loadaligns'),'Value');
string = get(findobj('Tag','listbox_loadaligns'),'String');

for i=1:size(string,1)
    if (i == vals) == zeros(1,size(vals,2))
        newstring = strvcat(newstring,string(i,:));
    end
end

set(findobj('Tag','listbox_loadaligns'),'String',newstring,'Value',1);
storage_alignmerge.filestring = newstring;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Load alignment files                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_merge_loadaligns_Callback(hObject, eventdata, handles)

global storage_alignmerge;

string = get(findobj('Tag','listbox_loadaligns'),'String');

for i=1:size(string,1)
    s = load(strtrim(string(i,:)));
    [classes,colors,radii] = get_classes(s.align2d);
    try align2d = rmfield(s.align2d,{'x' 'y'});s.align2d = align2d;end
    tmpalign = struct();
    for j=1:size(classes,2)
        storage_alignmerge.classes{size(storage_alignmerge.classes,2)+1} = classes{j};
        storage_alignmerge.radii{size(storage_alignmerge.radii,2)+1} = radii{j};
        storage_alignmerge.colors{size(storage_alignmerge.colors,2)+1} = colors{j};
    end
    for k=1:size(s.align2d,2)
        match = strmatch(s.align2d(k).class,classes,'exact');
        try
            tmpalign(match).align(1,size(tmpalign(match).align,2)+1) = s.align2d(1,k);
        catch
            tmpalign(match).align(1,1) = s.align2d(1,k);
        end
    end
    storage_alignmerge.inalign = [storage_alignmerge.inalign tmpalign];
end

string = '';

for i=1:size(storage_alignmerge.classes,2)
    string = strvcat(string,[cell2mat(storage_alignmerge.classes(i)) ' (radius: ' num2str(cell2mat(storage_alignmerge.radii(i))) ')']);
end

set(findobj('Tag','listbox_inputclasses'),'String',string);

set(findobj('Tag','listbox_inputclasses'),'Enable','on');
set(findobj('Tag','listbox_outputclasses'),'Enable','on');
set(findobj('Tag','button_merge_move'),'Enable','on');
set(findobj('Tag','button_merge_merge'),'Enable','on');

set(findobj('Tag','button_merge_addalign'),'Enable','off');
set(findobj('Tag','button_merge_deletealign'),'Enable','off');
set(findobj('Tag','button_merge_loadaligns'),'Enable','off');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Move class                                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_merge_move_Callback(hObject, eventdata, handles)

global storage_alignmerge;

%update listboxes
outputstring = get(findobj('Tag','listbox_outputclasses'),'String');
newstring  = '';
vals = get(findobj('Tag','listbox_inputclasses'),'Value');
string = get(findobj('Tag','listbox_inputclasses'),'String');
for i=1:size(string,1)
    if (i == vals) == zeros(1,size(vals,2))
        newstring = strvcat(newstring,string(i,:));
    else
        outputstring = strvcat(outputstring,string(i,:)); 
    end
end
set(findobj('Tag','listbox_inputclasses'),'String',newstring,'Value',1);
set(findobj('Tag','listbox_outputclasses'),'String',outputstring,'Value',1);

%update output alignment structure
storage_alignmerge.outalign = [storage_alignmerge.outalign storage_alignmerge.inalign(vals)];

%reorder input alignment structure
for i=1:size(storage_alignmerge.inalign,2)
    if i>vals
        storage_alignmerge.inalign(i-1).align = [];
        storage_alignmerge.inalign(i-1).align = storage_alignmerge.inalign(i).align;
    end
end
storage_alignmerge.inalign = storage_alignmerge.inalign(1:end-1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Merge Class                                                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_merge_merge_Callback(hObject, eventdata, handles)

global storage_alignmerge;

%update listboxes
newstring  = '';
vals = get(findobj('Tag','listbox_inputclasses'),'Value');
string = get(findobj('Tag','listbox_inputclasses'),'String');
for i=1:size(string,1)
    if (i == vals) == zeros(1,size(vals,2))
        newstring = strvcat(newstring,string(i,:));
    end
end
set(findobj('Tag','listbox_inputclasses'),'String',newstring,'Value',1);

%get output class to merge
outval = get(findobj('Tag','listbox_outputclasses'),'Value');

%update output alignment structure
for i=1:size(storage_alignmerge.inalign(vals).align,2)
    storage_alignmerge.inalign(vals).align(i).class = storage_alignmerge.outalign(outval).align(1).class;
    storage_alignmerge.inalign(vals).align(i).color = storage_alignmerge.outalign(outval).align(1).color;
end
storage_alignmerge.outalign(outval).align = [storage_alignmerge.outalign(outval).align storage_alignmerge.inalign(vals).align];

%reorder input alignment structure
for i=1:size(storage_alignmerge.inalign,2)
    if i>vals
        storage_alignmerge.inalign(i-1).align = [];
        storage_alignmerge.inalign(i-1).align = storage_alignmerge.inalign(i).align;
    end
end
storage_alignmerge.inalign = storage_alignmerge.inalign(1:end-1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Browse output alignment file                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_alignmerge_browseoutputalign_Callback(hObject, eventdata, handles)

[filename, pathname] = uiputfile({'*.mat'}, 'Save alignment file as');
if ischar(filename)
    set(findobj('Tag','alignmerge_outalignfile'),'String',[pathname '/' filename]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Browse output stack file                                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_alignmerge_browseoutputstack_Callback(hObject, eventdata, handles)

[filename, pathname] = uiputfile({'*.em'}, 'Save stack file as');
if ischar(filename)
    set(findobj('Tag','alignmerge_outstackfile'),'String',[pathname '/' filename]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Exit                                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_merge_exit_Callback(hObject, eventdata, handles)

global storage_alignmerge;

outfile = get(findobj('Tag','alignmerge_outalignfile'),'String');
outstack = get(findobj('Tag','alignmerge_outstackfile'),'String');

if ~isempty(outfile)
    align2d = '';
    for i=1:size(storage_alignmerge.outalign,2)
        align2d = [align2d storage_alignmerge.outalign(i).align];
    end
    save(outfile,'align2d');

    %create stack
    if ~isempty(outstack)
        normflag = get(findobj('Tag','checkbox_alignmerge_aligned'),'Value');
        alignedflag = get(findobj('Tag','alignmerge_normed'),'Value');
        tom_av2_createstack(outfile, '', outstack, outfile, normflag, alignedflag);
    end
    
else
    errordlg('Specify output alignment file first!');
end
setappdata(0,'UseNativeSystemDialogs',1);
close(gcf);


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
%%  Get classes in alignment file                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [newclasses,newcolors,newradii] = get_classes(align)
newclasses = {};
newcolors = {};
newradii = {};
for i=1:size(align,2)
    if strmatch(align(1,i).class,newclasses,'exact')
    else
        newclasses{size(newclasses,2)+1} = align(1,i).class;
        newcolors{size(newcolors,2)+1} = align(1,i).color;
        newradii{size(newradii,2)+1} = align(1,i).radius;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Unused Callbacks/Create functions                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function listbox_outputclasses_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function listbox_inputclasses_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function listbox_loadaligns_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function listbox_loadaligns_Callback(hObject, eventdata, handles)
function listbox_outputclasses_Callback(hObject, eventdata, handles)
function listbox_inputclasses_Callback(hObject, eventdata, handles)
function alignmerge_outalignfile_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function alignmerge_outalignfile_Callback(hObject, eventdata, handles)
function alignmerge_outstackfile_Callback(hObject, eventdata, handles)
function alignmerge_outstackfile_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function checkbox_alignmerge_aligned_Callback(hObject, eventdata, handles)
function alignmerge_normed_Callback(hObject, eventdata, handles)
function alignmerge_normed_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end