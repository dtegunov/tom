function varargout = tom_kerdensomgui(varargin)
%TOM_KERDENSOMGUI creates ...
%
%   varargout = tom_kerdensomgui(varargin)
%
%      TOM_KERDENSOMGUI, by itself, creates a new TOM_KERDENSOMGUI or raises the existing
%      singleton*.
%
%      H = TOM_KERDENSOMGUI returns the handle to a new TOM_KERDENSOMGUI or the handle to
%      the existing singleton*.
%
%      TOM_KERDENSOMGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TOM_KERDENSOMGUI.M with the given input arguments.
%
%      TOM_KERDENSOMGUI('Property','Value',...) creates a new TOM_KERDENSOMGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before tom_kerdensomgui_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property
%      application
%      stop.  All inputs are passed to tom_kerdensomgui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
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
%   GUIDE, GUIDATA, GUIHANDLES
%
%   created by ... (author date)
%   updated by ...
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
% Edit the above text to modify the response to help tom_kerdensomgui
%
% Last Modified by GUIDE v2.5 29-Mar-2007 17:07:27
%
% Begin initialization code - DO NOT EDIT

% error(nargchk(0, 1, nargin, 'struct'))

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_kerdensomgui_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_kerdensomgui_OutputFcn, ...
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
function tom_kerdensomgui_OpeningFcn(hObject, eventdata, handles, varargin)

handles = show_netaxes(handles,'off');

axes(handles.kerdensom_browser);
axis off;
handles = unset_radiobuttons(handles);
handles.rectangles = [];
handles.selectedpoints = [];
handles.error = [];
handles.storage.thresholds = [];
handles.selectedclasses = [];
handles.lines = [];

handles = change_selectpanel(handles,'off');
handles = change_amirapanel(handles,'off');
handles = change_visupanel(handles,'off');
handles = change_calcpanel(handles,'off');
handles = change_netpanel(handles,'off');
handles.avgcache = struct();
handles.storage.net = struct();


handles.showiteration = 1;

% Choose default command line output for tom_kerdensomgui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% uiwait(handles.figure1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Output Function                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = tom_kerdensomgui_OutputFcn(hObject, eventdata, handles) 

% Get default command line output from handles structure
varargout{1} = handles.output;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  input alignment file                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function kerdensom_alignmentfile_Callback(hObject, eventdata, handles)

handles = loadfile(handles);

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  browse input alignment file                                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function kerdensom_button_browsealignmentfile_Callback(hObject, eventdata, handles)

[filename, pathname] = uigetfile({'*.mat'}, 'Pick an alignment file');
if ischar(filename)
    set(handles.kerdensom_alignmentfile,'String',[pathname '/' filename]);
    handles = loadfile(handles);
end

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  generate alignment file                                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function kerdensom_button_generatealignmentfile_Callback(hObject, eventdata, handles)

filename = tom_av3_createalign();

if ischar(filename)
   set(handles.kerdensom_alignmentfile,'String',filename); 
   handles = loadfile(handles);
end

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  load alignment file                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function kerdensom_button_loadalignmentfile_Callback(hObject, eventdata, handles)

contents = get(handles.kerdensom_run,'String');
run = str2num(contents(get(handles.kerdensom_run,'Value')));

set(handles.kerdensom_run,'Enable','off');
set(handles.kerdensom_binning,'Enable','off');
binning = str2num(get(handles.kerdensom_binning,'String'));
handles.storage.binning = binning;
handles.storage.run = run;

h = waitbar(0,'Reading files...');

%1d
if handles.storage.dims == 1
    handles.storage.stackfilename = get(handles.kerdensom_alignmentfile,'String');
    header = tom_reademheader(handles.storage.stackfilename);
    handles.storage.Header = header.Header;
    handles.storage.imstack = tom_emreadc(handles.storage.stackfilename);
    handles.storage.imstack = double(handles.storage.imstack.Value);
%2d
elseif handles.storage.dims == 2
    handles.storage.align = handles.storage.alignallruns(run,:);   
    [pathstr, name, ext] = fileparts(get(handles.kerdensom_alignmentfile,'String'));
    handles.storage.stackfilename = [pathstr '/' name '.em'];
    header = tom_reademheader(handles.storage.stackfilename);
    handles.storage.Header = header.Header;
    handles.storage.imstack = tom_reshape_stack(handles.storage.stackfilename,'',binning);
    
%3d    
elseif handles.storage.dims == 3
    handles.storage.align = tom_av3_align_sum(handles.storage.alignallruns(1:run,:));
    for i=1:size(handles.storage.align,2)
        handles.storage.align(1,i).Shift.X = handles.storage.align(1,i).Shift.X ./ 2^binning;
        handles.storage.align(1,i).Shift.Y = handles.storage.align(1,i).Shift.Y ./ 2^binning;        
        handles.storage.align(1,i).Shift.Z = handles.storage.align(1,i).Shift.Z ./ 2^binning;
    end
    
    filecell = {};
    for i=1:size(handles.storage.align,2)
        filecell{i} = handles.storage.align(i).Filename;
    end
    header = tom_reademheader(filecell{1});
    handles.storage.Header = header.Header;
    in_struct.rotmask.types = {'sphere3d','cylinder3d'};
    maskstruct = tom_filtergui('mask',in_struct);
    maskstruct.Apply=2;
    [handles.storage.imstack handles.storage.ang_out] = tom_reshape_stack(filecell,handles.storage.align,binning,'',0,maskstruct);
   
else
    errordlg('Alignment file not valid!');
    close(h);
    return;
end
waitbar(0.5,h,'Removing mean...');
%[handles.storage.imstack handles.storage.mean]=tom_rm_mean(handles.storage.imstack);
handles.cachevalid = 0;
close(h);

handles = change_calcpanel(handles,'on');
handles.storage.thresholds = [];
handles.selectedclasses = [];

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Execute calculation                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_kerdensomgui_view_Callback(hObject, eventdata, handles)

imstack = handles.storage.imstack;
for i=1:size(handles.storage.imstack,1)
        imstack(i,:)=tom_norm(imstack(i,:),'mean0+1std');
end

tom_stattool(imstack,'statmatrix');

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Execute calculation                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function kerdensom_button_calculate_Callback(hObject, eventdata, handles)

contents = get(handles.kerdensom_run,'String');
run = str2num(contents(get(handles.kerdensom_run,'Value')));
binning = str2num(get(handles.kerdensom_binning,'String'));


contents = get(handles.popupmenu_algo,'String');
algo = contents{get(handles.popupmenu_algo,'Value')};


h = waitbar(0,'Calculating kerdensom...');

outputMapHeight = str2num(get(handles.kerdensom_input_gridheight,'String'));
outputMapWidth = str2num(get(handles.kerdensom_input_gridwidth,'String'));
initialReg = str2num(get(handles.kerdensom_input_initialreg,'String'));
finalReg = str2num(get(handles.kerdensom_input_finalreg,'String'));
regSteps = str2num(get(handles.kerdensom_input_regsteps,'String'));
som_Nsteps = str2num(get(handles.kerdensom_input_som_Nsteps,'String'));
initialMap = [];

handles.storage.gridsize = [outputMapWidth outputMapHeight];
handles.storage.regSteps = regSteps;



    if isequal(algo,'SOM')
        handles.storage.result = tom_som(handles.storage.imstack, [outputMapHeight, outputMapWidth], regSteps,1);
    else
        if get(handles.wedge_it,'Value')==1
           handles.storage.result = kerdensom3(handles.storage.imstack',outputMapHeight,  outputMapWidth, initialReg, finalReg, regSteps, som_Nsteps, initialMap,2,handles.storage.ang_out);
        else
           handles.storage.result = kerdensom(handles.storage.imstack',outputMapHeight,  outputMapWidth, initialReg, finalReg, regSteps, som_Nsteps, initialMap,2);
        end;
   end;


    



    
for i=1:regSteps
    
    if handles.storage.dims == 2
        handles.storage.result(i).outputMap=reshape(handles.storage.result(i).outputMap,[handles.storage.Header.Size(1)./ 2^binning handles.storage.Header.Size(2)./ 2^binning outputMapWidth.*outputMapHeight]);
    else
        handles.storage.result(i).outputMap=reshape(handles.storage.result(i).outputMap,[handles.storage.Header.Size(1)./ 2^binning handles.storage.Header.Size(2)./ 2^binning handles.storage.Header.Size(3)./ 2^binning outputMapWidth.*outputMapHeight]);
    end


    if handles.storage.dims == 3
        for ii=1:size(handles.storage.result(1).assignVtoX,1)
            handles.storage.align(i,ii).Class = handles.storage.result(i).assignVtoX(ii);
        end
    else
        for ii=1:size(handles.storage.result(1).assignVtoX,1)
            handles.storage.align(i,ii).Class = handles.storage.result(i).assignVtoX(ii);
        end
    end
end



string = '';
for i=1:regSteps
    string = strvcat(string,num2str(i));
end
set(handles.popupmenu_kerdensom_iteration,'String',string,'Value',1);
handles.showiteration = 1;
set(handles.edit_movie_stop,'String',num2str(regSteps));
close(h);
handles = init_avgcache(handles);
handles = change_visupanel(handles,'on');

set(handles.edit_origin,'String',[num2str(floor(outputMapWidth./2+1)) ',' num2str(floor(outputMapHeight./2+1))]);

handles = calc_net(handles);
handles.storage.thresholds = [];
handles.selectedclasses = [];
handles.rectangles = [];
handles.selectedpoints = [];
guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Load calculation                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function kerdensom_button_loadcalc_Callback(hObject, eventdata, handles)


[filename, pathname] = uigetfile('*.mat', 'Load Calculation');
if ischar(filename)
    s = load([pathname '/' filename]);

    if isfield(s,'kerdensom') ~= 1
        errordlg('This is not a kerdensom calculation file.','File Error');
    else
        handles.storage = s.kerdensom;
        try
            set(handles.kerdensom_alignmentfile,'String',handles.storage.stackfilename);
        end
        set(handles.kerdensom_binning,'String',num2str(handles.storage.binning));
        %update runs dropdown field
        liststring = '';
        for i=1:size(handles.storage.alignallruns,1)
            liststring = strvcat(liststring,num2str(i));
        end
        set(handles.kerdensom_run,'String',liststring);
        set(handles.kerdensom_run,'Value',handles.storage.run);
        set(handles.kerdensom_run,'Enable','off');
        set(handles.kerdensom_binning,'Enable','off');
        string = '';
        for i=1:size(handles.storage.result,2)
            string = strvcat(string,num2str(i));
        end
        set(handles.popupmenu_kerdensom_iteration,'String',string,'Value',1);
        handles.showiteration = 1;
        set(handles.edit_movie_stop,'String',num2str(size(handles.storage.result,2)));
        
        if isfield(s,'avgcache')
            handles.avgcache = s.avgcache;
        else
            handles = init_avgcache(handles);
        end
    end
end
handles = unset_radiobuttons(handles);

if handles.storage.dims == 3
    handles = change_amirapanel(handles,'on');
else
    handles = change_amirapanel(handles,'off');
end

set(handles.edit_origin,'String',[num2str(floor(handles.storage.gridsize(1)./2+1)) ',' num2str(floor(handles.storage.gridsize(2)./2+1))]);
handles = calc_net(handles);
handles.selectedclasses = [];
handles.rectangles = [];
handles.selectedpoints = [];
handles = change_visupanel(handles,'on');
handles = change_calcpanel(handles,'on');

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Save calculation                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function kerdensom_button_savecalc_Callback(hObject, eventdata, handles)

[filename, pathname] = uiputfile('*.mat', 'Save Calculation as');
if ischar(filename)
    
    if get(handles.checkbox_kerdensom_saveavgcache,'Value') == 1
        avgcache = handles.avgcache;
    end
    kerdensom = handles.storage;
    if ~exist('avgcache');
        avgcache = [];
    end
    save([pathname '/' filename],'-v7.3','kerdensom','avgcache');
end

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display codevectors                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function display_codevectors_Callback(hObject, eventdata, handles)

handles = show_netaxes(handles,'off');
set(handles.disp_gallery,'Value',0);
set(handles.display_codevectors,'Value',1);
set(handles.disp_colormap,'Value',0);
set(handles.display_colormap3d,'Value',0);
set(handles.display_net,'Value',0);
set(handles.display_variancemap,'Value',0);
handles = change_selectpanel(handles,'on');
handles = change_netpanel(handles,'off');
handles = disp_codvecs(handles);

handles.rectangles = [];
handles.selectedpoints = [];

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display gallery                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function disp_gallery_Callback(hObject, eventdata, handles)

handles = show_netaxes(handles,'off');
set(handles.disp_gallery,'Value',1);
set(handles.display_codevectors,'Value',0);
set(handles.disp_colormap,'Value',0);
set(handles.display_colormap3d,'Value',0);
set(handles.display_net,'Value',0);
set(handles.display_variancemap,'Value',0);
handles = change_selectpanel(handles,'on');
handles = change_netpanel(handles,'off');
handles = disp_avggallery(handles);

handles.rectangles = [];
handles.selectedpoints = [];

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display colormap                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function disp_colormap_Callback(hObject, eventdata, handles)

handles = show_netaxes(handles,'off');
set(handles.disp_gallery,'Value',0);
set(handles.display_codevectors,'Value',0);
set(handles.disp_colormap,'Value',1);
set(handles.display_colormap3d,'Value',0);
set(handles.display_net,'Value',0);
set(handles.display_variancemap,'Value',0);
handles = change_selectpanel(handles,'off');
handles = change_netpanel(handles,'off');
handles = disp_colormap(handles);

handles.rectangles = [];
handles.selectedpoints = [];

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display colormap 3d                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function display_colormap3d_Callback(hObject, eventdata, handles)

handles = show_netaxes(handles,'off');
set(handles.disp_gallery,'Value',0);
set(handles.display_codevectors,'Value',0);
set(handles.disp_colormap,'Value',0);
set(handles.display_colormap3d,'Value',1);
set(handles.display_variancemap,'Value',0);
set(handles.display_net,'Value',0);

handles = change_selectpanel(handles,'off');
handles = change_netpanel(handles,'off');
handles = disp_colormap3d(handles);

handles.rectangles = [];
handles.selectedpoints = [];

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display distance map                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function display_net_Callback(hObject, eventdata, handles)

handles = show_netaxes(handles,'on');
set(handles.disp_gallery,'Value',0);
set(handles.display_codevectors,'Value',0);
set(handles.disp_colormap,'Value',0);
set(handles.display_colormap3d,'Value',0);
set(handles.display_net,'Value',1);
set(handles.display_variancemap,'Value',0);
handles = change_selectpanel(handles,'on');

handles = change_netpanel(handles,'on');

handles = disp_distmap(handles);

handles.rectangles = [];
handles.selectedpoints = [];

handles = plot_lines(handles);

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display variance map                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function display_variancemap_Callback(hObject, eventdata, handles)

handles = show_netaxes(handles,'off');
set(handles.disp_gallery,'Value',0);
set(handles.display_codevectors,'Value',0);
set(handles.disp_colormap,'Value',0);
set(handles.display_colormap3d,'Value',0);
set(handles.display_net,'Value',0);
set(handles.display_variancemap,'Value',1);
handles = change_selectpanel(handles,'on');

handles = disp_variancemap(handles);

handles.rectangles = [];
handles.selectedpoints = [];

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  distance map metric                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function popupmenu_kerdensom_metric_Callback(hObject, eventdata, handles)

handles = calc_net(handles);
handles = disp_distmap(handles);
handles = plot_lines(handles);
guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  distance map origin                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edit_origin_Callback(hObject, eventdata, handles)

handles = calc_net(handles);
handles = disp_distmap(handles);
handles = plot_lines(handles);
guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  selection 2 stackbrowser                                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_selection2stackbrowser_Callback(hObject, eventdata, handles)

handles = saveselection(handles,[pwd '/tmp.mat']);

if (handles.error == 1)
    handles.error = [];
    guidata(hObject, handles);
    return;
end

if handles.storage.dims == 2 || handles.storage.dims == 1
    tom_av2_stackbrowser([pwd '/tmp.em'],[pwd '/tmp.mat']);
elseif handles.storage.dims == 3 
    tom_av3_stackbrowser([pwd '/tmp.mat']);
end

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  save selection                                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_saveselection_Callback(hObject, eventdata, handles)

handles = saveselection(handles);

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  select classes                                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_select_Callback(hObject, eventdata, handles)

if get(handles.display_net,'Value') == 1
    
    if get(handles.checkbox_net_classavgs,'Value') == 1
        type = 'avg';
    else
        type = 'cod';
    end
    
    axes(handles.kerdensom_net);
    coords = ginput(1);
    
    %determine which particle was picked
    [class,node] = tom_nearestpoint(coords, handles.storage.net(handles.showiteration).net);
    %[xxx,miniview] = tom_nearestpoint(node, handles.miniviewpoints);
    
    coords = ginput(1);
    axes2 = get(gca,'Tag');
    xlim = get(handles.kerdensom_net,'Xlim');
    ylim = get(handles.kerdensom_net,'Ylim');
    
    %calculate line coordinates
    if ~isempty(strfind(axes2,'axes_kerdensom_right')) == 1
        linestop = [xlim(2) handles.ypoints(str2num(axes2(end)))];
    elseif ~isempty(strfind(axes2,'axes_kerdensom_left')) == 1
        linestop = [xlim(1) handles.ypoints(str2num(axes2(end)))];
    elseif ~isempty(strfind(axes2,'axes_kerdensom_top')) == 1            
        linestop = [handles.xpoints(str2num(axes2(end))) ylim(2)];
    elseif ~isempty(strfind(axes2,'axes_kerdensom_bottom')) == 1
        linestop = [handles.xpoints(str2num(axes2(end))) ylim(1)];
    end

    l = size(handles.lines,2)+1;
    handles.lines(l).class = class;
    handles.lines(l).endpoint = linestop;
    handles.lines(l).type = type;
    handles.lines(l).axes = axes2;

    handles = plot_lines(handles);
    
else
    axes(handles.kerdensom_browser);
    coords = ginput(1);

    if get(handles.disp_gallery,'Value') == 1
        binning = 0;
    else
        binning = str2num(get(handles.kerdensom_binning,'String'));
    end


    %determine which particle was picked
    numcols = str2double(get(handles.kerdensom_input_gridwidth,'String'));
    row = ceil(coords(2)./handles.storage.Header.Size(2).*2^binning);
    column = ceil(coords(1)./handles.storage.Header.Size(1).*2^binning);
    classno = (row-1).*numcols+column-1;

    %draw rectangle
    handles = disp_mark(handles,classno,binning);
    handles.selectedclasses = [handles.selectedclasses classno];

    %select particles in this class in the alignment structure
    for i=1:size(handles.storage.align,2)
        if handles.storage.align(handles.showiteration,i).Class == classno+1
            handles.selectedpoints = [handles.selectedpoints i];
        end
    end
end

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  reset selection                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_reset_Callback(hObject, eventdata, handles)

for i=1:length(handles.rectangles)
   try
       delete(handles.rectangles(i)); 
   end
end

handles.rectangles = [];
handles.selectedpoints = [];
handles.selectedclasses = [];


if get(handles.display_net,'Value') == 1
    handles.lines = [];
    handles = show_netaxes(handles,'off');
    handles = show_netaxes(handles,'on');
    handles = disp_distmap(handles);
end

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   set mass                                                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_setmass_Callback(hObject, eventdata, handles)

molmass = str2num(get(handles.edit_mass,'String'));
pixelsize = str2num(get(handles.edit_pixelsize,'String'));


if isempty(molmass)
    errordlg('No molmass given');
    return;
end

if isempty(pixelsize)
    errordlg('No pixelsize given');
    return;
end


for i=1:size(handles.selectedpoints,2)
    al(i) = handles.storage.align(1,handles.selectedpoints(i));
end

if exist('al','var') && ~isempty(al(1,1).Tomogram)
    volume = tom_av3_average(al,'sum',0,al(1,end).Class,1);
    if get(handles.checkbox_inverse,'Value') == 1
        volume = -volume;
    end
    
    handles.storage.thresholds(handles.selectedclasses+1)=tom_calc_isosurface(tom_norm(volume,1),molmass,pixelsize,0.01);
end

if get(handles.disp_gallery,'Value') == 1
    binning = 0;
else
    binning = str2num(get(handles.kerdensom_binning,'String'));
end
handles = disp_particlenumbers(handles,binning);

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  matrix 2 amira                                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_matrix2amira_Callback(hObject, eventdata, handles)

host = get(handles.edit_kerdensom_host,'String');
voltexflag = get(handles.checkbox_kerdensom_voltex,'Value');

if get(handles.checkbox_inverse,'Value') == 1
    inverse = 1;
else
    inverse = 0;
end

if voltexflag == 1
    tom_amira_load_seismic(host);
else
    if length(handles.storage.thresholds) < handles.storage.gridsize(1).*handles.storage.gridsize(2)
        handles.storage.thresholds(handles.storage.gridsize(1).*handles.storage.gridsize(2)) = 0;
    end
end



if get(handles.disp_gallery,'Value') == 1
    
    if voltexflag == 1
        map = handles.avgcache(handles.showiteration).avg;
        for i=1:size(map,4)
            map(:,:,:,i) = tom_filter(map(:,:,:,i),2);
            map(:,:,:,i) = tom_limit(tom_norm(map(:,:,:,i),'3std')+.5,.1,.8);

        end
        tom_amira_displayisosurfacematrix(map, handles.storage.result(handles.showiteration).histogram, handles.storage.gridsize, inverse, 0.5,'cool', voltexflag, host);
    else
        map = handles.avgcache(handles.showiteration).avg;
        tom_amira_displayisosurfacematrix(map, handles.storage.result(handles.showiteration).histogram, handles.storage.gridsize, inverse, handles.storage.thresholds,'cool', voltexflag, host);
    end
    
    
else
    if voltexflag == 1
        map = handles.storage.result(handles.showiteration).outputMap;
        for i=1:size(map,4)
            map(:,:,:,i) = tom_filter(map(:,:,:,i),2);
            map(:,:,:,i) = tom_limit(tom_norm(map(:,:,:,i),'3std')+.5,.1,.8);
        end
        tom_amira_displayisosurfacematrix(map, handles.storage.result(handles.showiteration).histogram, handles.storage.gridsize, inverse, 0.5,'cool', voltexflag, host);
    else
        for  j=1:size(handles.storage.result(handles.showiteration).outputMap,4)
            map(:,:,:,j) = tom_norm(handles.storage.result(handles.showiteration).outputMap(:,:,:,j),1);
        end
        tom_amira_displayisosurfacematrix(map, handles.storage.result(handles.showiteration).histogram, handles.storage.gridsize, inverse, handles.storage.thresholds,'cool', voltexflag, host);
    end
    
end

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  selection 2 dspcub                                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_kerdensom_selection2dspcub_Callback(hObject, eventdata, handles)

if get(handles.checkbox_inverse,'Value') == 1
    inverse = 1;
else
    inverse = 0;
end

for i=1:size(handles.selectedpoints,2)
 %for i=1:350
    al(i) = handles.storage.align(1,handles.selectedpoints(i));
end

if get(handles.disp_gallery,'Value') == 1
    if exist('al','var') && ~isempty(al(1,1).Tomogram)
%    averg = tom_av3_average(al,'sum',0,al(1,end).Class,1);
    for i=1:size(al,2);
        al(1,i).Class=0;
    end;
    averg = tom_av3_average(al,'sum',-1,0,1);
    end
else
    if handles.storage.dims == 3
        averg = zeros(size(handles.storage.result(1).outputMap,1),size(handles.storage.result(1).outputMap,2),size(handles.storage.result(1).outputMap,3));
       for i=1:size(handles.selectedclasses,2)
            averg = averg + handles.storage.result(handles.showiteration).outputMap(:,:,:,handles.selectedclasses(i)+1);
       end
    else
        averg = zeros(size(handles.storage.result(1).outputMap,1),size(handles.storage.result(1).outputMap,2));
        averg = averg + handles.storage.result(handles.showiteration).outputMap(:,:,handles.selectedclasses(i)+1);
    end
end

if exist('averg')
    if inverse == 1
        averg = -averg;
    end

    figure;tom_dspcub(averg,1);
end    

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  selection 2 tom_volxyz                                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_kerdensom_selection2volxyz_Callback(hObject, eventdata, handles)

if get(handles.checkbox_inverse,'Value') == 1
    inverse = 1;
else
    inverse = 0;
end

for i=1:size(handles.selectedpoints,2)
    al(i) = handles.storage.align(1,handles.selectedpoints(i));

end

if get(handles.disp_gallery,'Value') == 1
    if exist('al','var') && ~isempty(al(1,1).Tomogram)
    averg = tom_av3_average(al,'sum',0,al(1,end).Class,1);
    end
else
    if handles.storage.dims == 3
        averg = zeros(size(handles.storage.result(1).outputMap,1),size(handles.storage.result(1).outputMap,2),size(handles.storage.result(1).outputMap,3));
       for i=1:size(handles.selectedclasses,2)
            averg = averg + handles.storage.result(handles.showiteration).outputMap(:,:,:,handles.selectedclasses(i)+1);
       end
    else
        averg = zeros(size(handles.storage.result(1).outputMap,1),size(handles.storage.result(1).outputMap,2));
        averg = averg + handles.storage.result(handles.showiteration).outputMap(:,:,handles.selectedclasses(i)+1);
    end
end

if exist('averg')
    if inverse == 1
        averg = -averg;
    end

    tom_volxyz(averg);
end    

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  selection 2 amira                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_selection2amira_Callback(hObject, eventdata, handles)

host = get(handles.edit_kerdensom_host,'String');
voltexflag = get(handles.checkbox_kerdensom_voltex,'Value');
    
if get(handles.checkbox_inverse,'Value') == 1
    inverse = 1;
else
    inverse = 0;
end

if length(handles.storage.thresholds) < handles.storage.gridsize(1).*handles.storage.gridsize(2)
    handles.storage.thresholds(handles.storage.gridsize(1).*handles.storage.gridsize(2)) = 0;
end

for i=1:size(handles.selectedpoints,2)
    al(i) = handles.storage.align(1,handles.selectedpoints(i));
end

if get(handles.disp_gallery,'Value') == 1
    if exist('al','var') && ~isempty(al(1,1).Tomogram)
    averg = tom_av3_average(al,'sum',0,al(1,end).Class,1);
    
    end
else
    if handles.storage.dims == 3
        averg = zeros(size(handles.storage.result(1).outputMap,1),size(handles.storage.result(1).outputMap,2),size(handles.storage.result(1).outputMap,3));
       for i=1:size(handles.selectedclasses,2)
            averg = averg + handles.storage.result(handles.showiteration).outputMap(:,:,:,i);
       end
    else
        averg = zeros(size(handles.storage.result(1).outputMap,1),size(handles.storage.result(1).outputMap,2));
        averg = averg + handles.storage.result(handles.showiteration).outputMap(:,:,i);
    end
end

if exist('averg')
    if voltexflag == 1
        averg = tom_filter(averg,2);
        averg = tom_limit(tom_norm(averg,'3std')+.5,.1,.8);
    else
        averg = tom_norm(averg,1);
    end
    
    if inverse == 1
        tom_emwritec([pwd filesep 'ktmp.em'],-averg);
    else
        tom_emwritec([pwd filesep 'ktmp.em'],averg);
    end
    tom_amira_loadfile([pwd filesep 'ktmp.em'],'avg', [], host);

    if voltexflag == 0
        tom_amira_createisosurface('avg',handles.storage.thresholds(handles.selectedclasses(1)), 'iavg',[],eye(4),[],host);
    else
        tom_amira_load_seismic(host);
        tom_amira_createvoltex('avg',[], 'iavg',[],eye(4),[],host);
    end
end

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display iteration                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function popupmenu_kerdensom_iteration_Callback(hObject, eventdata, handles)

contents = get(hObject,'String');
run = str2num(contents(get(hObject,'Value'),:));

handles.showiteration = run;

if get(handles.display_codevectors,'Value') == 1
    handles = disp_codvecs(handles);drawnow;
    %movie of avg gallery
elseif get(handles.disp_gallery,'Value') == 1
    handles = disp_avggallery(handles);drawnow;
    %movie of colormap
elseif get(handles.disp_colormap,'Value') == 1
    handles = disp_colormap(handles);drawnow;
    %movie of colormap 3d
elseif get(handles.display_colormap3d,'Value') == 1
    handles = disp_colormap3d(handles);drawnow;
    %movie of net
elseif get(handles.display_net,'Value') == 1
    handles = disp_distmap(handles);
    handles = plot_lines(handles);drawnow;
    %movie of variance map
elseif get(handles.display_variancemap,'Value') == 1
    handles = disp_variancemap(handles);
end

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  movie play                                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_kerdensom_movie_Callback(hObject, eventdata, handles)

pt = str2double(get(handles.edit_kerdensom_movietime,'String'));
pt = 1./pt;

start = str2double(get(handles.edit_movie_start,'String'));
stop = str2double(get(handles.edit_movie_stop,'String'));

counter = 0;

for i=start:stop
    if isempty(handles.avgcache) || size(handles.avgcache,2) < i || isempty(handles.avgcache(i).isvalid) || handles.avgcache(i).isvalid == 0
    else
        counter = counter + 1;
    end
end

%generate avg cache
if get(handles.disp_gallery,'Value') == 1 && counter < size(handles.storage.result,2)
    for i=1:size(handles.storage.result,2)
        set(handles.popupmenu_kerdensom_iteration,'Value',i);
        handles.showiteration = i;
        handles = generate_avgs(handles);
        handles.avgcache(i).avg = handles.avg;
        handles.avgcache(i).isvalid = 1;
    end
end


for i=start:stop

    %movie of output Map
    if get(handles.display_codevectors,'Value') == 1
        set(handles.popupmenu_kerdensom_iteration,'Value',i);
        handles.showiteration = i;
        handles = disp_codvecs(handles);drawnow;
        pause(pt);
    %movie of avg gallery
    elseif get(handles.disp_gallery,'Value') == 1
        set(handles.popupmenu_kerdensom_iteration,'Value',i);
        handles.showiteration = i;
        handles = disp_avggallery(handles);drawnow;
        pause(pt);
    %movie of colormap    
    elseif get(handles.disp_colormap,'Value') == 1
        set(handles.popupmenu_kerdensom_iteration,'Value',i);
        handles.showiteration = i;
        handles = disp_colormap(handles);drawnow;
        pause(pt);
    %movie of colormap 3d   
    elseif get(handles.display_colormap3d,'Value') == 1
        set(handles.popupmenu_kerdensom_iteration,'Value',i);
        handles.showiteration = i;
        handles = disp_colormap3d(handles);drawnow;
        pause(pt);
    %movie of net   
    elseif get(handles.display_net,'Value') == 1
        set(handles.popupmenu_kerdensom_iteration,'Value',i);
        handles.showiteration = i;
        handles = disp_distmap(handles);
        handles = plot_lines(handles);drawnow;
        pause(pt);
    %movie of variance map
    elseif get(handles.display_variancemap,'Value') == 1
        set(handles.popupmenu_kerdensom_iteration,'Value',i);
        handles.showiteration = i;
        handles = disp_variancemap(handles);drawnow;
        pause(pt);

    end
    
end    
    
guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  movie save                                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_kerdensomsavemovie_Callback(hObject, eventdata, handles)

fps = str2double(get(handles.edit_kerdensom_movietime,'String'));
fps = fps;
start = str2double(get(handles.edit_movie_start,'String'));
stop = str2double(get(handles.edit_movie_stop,'String'));
counter = 0;

for i=start:stop
    if isempty(handles.avgcache) || size(handles.avgcache,2) < i || isempty(handles.avgcache(i).isvalid) || handles.avgcache(i).isvalid == 0
    else
        counter = counter + 1;
    end
end

[filename, pathname] = uiputfile({'*.avi'}, 'Save as avi file');
if ~ischar(filename)
    return;
end

set(gcf,'DoubleBuffer','on');
mov = avifile([pathname '/' filename]);
mov.Fps=fps;
if isequal(computer,'PCWIN') || isequal(computer,'PCWIN64')
    mov.Compression='Cinepak';
else
    mov.Compression='None';
end
% 'Cinepak', 'Indeo3', 'Indeo5', 'MSVC',', 'RLE', 'None'
mov.Quality=100;

%generate avg cache
if get(handles.disp_gallery,'Value') == 1 && counter < size(handles.storage.result,2)
    for i=1:size(handles.storage.result,2)
        set(handles.popupmenu_kerdensom_iteration,'Value',i);
        handles.showiteration = i;
        handles = generate_avgs(handles);
        handles.avgcache(i).avg = handles.avg;
        handles.avgcache(i).isvalid = 1;
    end
end


for i=start:stop

    %movie of output Map
    if get(handles.display_codevectors,'Value') == 1
        set(handles.popupmenu_kerdensom_iteration,'Value',i);
        handles.showiteration = i;
        handles = disp_codvecs(handles);drawnow;
        F = getframe(gca);
        mov = addframe(mov,F);
    %movie of avg gallery
    elseif get(handles.disp_gallery,'Value') == 1
        set(handles.popupmenu_kerdensom_iteration,'Value',i);
        handles.showiteration = i;
        handles = disp_avggallery(handles);drawnow;
        F = getframe(gca);
        mov = addframe(mov,F);
    %movie of colormap    
    elseif get(handles.disp_colormap,'Value') == 1
        set(handles.popupmenu_kerdensom_iteration,'Value',i);
        handles.showiteration = i;
        handles = disp_colormap(handles);drawnow;
        F = getframe(gca);
        mov = addframe(mov,F);
    %movie of colormap 3d   
    elseif get(handles.display_colormap3d,'Value') == 1
        set(handles.popupmenu_kerdensom_iteration,'Value',i);
        handles.showiteration = i;
        handles = disp_colormap3d(handles);drawnow;
        F = getframe(gca);
        mov = addframe(mov,F);
    %movie of net 
    elseif get(handles.display_net,'Value') == 1
        set(handles.popupmenu_kerdensom_iteration,'Value',i);
        handles.showiteration = i;
        handles = disp_distmap(handles);drawnow;
        F = getframe(gca);
        mov = addframe(mov,F);    
    %movie of variance map
    elseif get(handles.display_variancemap,'Value') == 1
        set(handles.popupmenu_kerdensom_iteration,'Value',i);
        handles.showiteration = i;
        handles = disp_variancemap(handles);drawnow;
        F = getframe(gca);
        mov = addframe(mov,F);    
    end
    
end       

mov = close(mov);    

msgbox('Movie written');

guidata(hObject, handles);


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
%%  make net parameters  active/inactive                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = change_netpanel(handles,val)

set(handles.edit_origin,'Visible',val);
set(handles.popupmenu_kerdensom_metric,'Visible',val);
set(handles.text_metric,'Visible',val);
set(handles.text_origin,'Visible',val);
set(handles.checkbox_net_classavgs,'Visible',val);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  make calculation panel active/inactive                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = change_calcpanel(handles,val)


set(handles.popupmenu_algo,'Enable',val);
set(handles.kerdensom_input_gridheight,'Enable',val);
set(handles.kerdensom_input_gridwidth,'Enable',val);
set(handles.kerdensom_input_initialreg,'Enable',val);
set(handles.kerdensom_input_finalreg,'Enable',val);
set(handles.kerdensom_input_regsteps,'Enable',val);
set(handles.kerdensom_input_som_Nsteps,'Enable',val);
set(handles.kerdensom_button_calculate,'Enable',val);
%set(handles.kerdensom_button_loadcalc,'Enable',val);
set(handles.kerdensom_button_savecalc,'Enable',val);
set(handles.checkbox_kerdensom_saveavgcache,'Enable',val);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  make select panel active/inactive                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = change_selectpanel(handles,val)

set(handles.button_select,'Enable',val);
set(handles.button_reset,'Enable',val);
set(handles.button_selection2stackbrowser,'Enable',val);
set(handles.button_saveselection,'Enable',val);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  make Amira panel active/inactive                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = change_amirapanel(handles,val)

set(handles.button_matrix2amira,'Enable',val);
set(handles.button_selection2amira,'Enable',val);
set(handles.checkbox_inverse,'Enable',val);
set(handles.button_kerdensom_selection2dspcub,'Enable',val);
set(handles.button_kerdensom_selection2volxyz,'Enable',val);
set(handles.checkbox_kerdensom_voltex,'Enable',val);
set(handles.edit_kerdensom_host,'Enable',val);
set(handles.edit_pixelsize,'Enable',val);
set(handles.edit_mass,'Enable',val);
set(handles.button_setmass,'Enable',val);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  make visualization panel active/inactive                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = change_visupanel(handles,val)

set(handles.display_codevectors,'Enable',val);
set(handles.disp_gallery,'Enable',val);
set(handles.disp_colormap,'Enable',val);
set(handles.display_colormap3d,'Enable',val);
set(handles.display_net,'Enable',val);
set(handles.display_variancemap,'Enable',val);
set(handles.popupmenu_kerdensom_iteration,'Enable',val);
set(handles.edit_kerdensom_movietime,'Enable',val);
set(handles.edit_movie_start,'Enable',val);
set(handles.edit_movie_stop,'Enable',val);
set(handles.button_kerdensom_movie,'Enable',val);
set(handles.button_kerdensomsavemovie,'Enable',val);
try
    if handles.storage.dims == 3
        handles = change_amirapanel(handles,val);
    else
        handles = change_amirapanel(handles,'off');
    end
catch
    handles = change_amirapanel(handles,'off');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  save selection to file                                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = saveselection(handles, name)

if handles.storage.dims == 1
    [filename, pathname] = uigetfile({'*.mat'}, 'Load alignment file');
    if ~isempty(filename)
        s = load([pathname '/' filename]);
        handles.storage.align = s.align2d;
        [pathstr, name, ext] = fileparts([pathname '/' filename]);
        stackfile2 = [pathstr '/' name '.em'];
    end
end

try;handles.storage.align = rmfield(handles.storage.align,'selectedtowrite');end;
indd = handles.selectedpoints;



if isempty(indd)
    errordlg('No particles selected!');
    handles.error = 1;
    return;
end


al_tmp = handles.storage.align(1,:);
for i=1:length(indd)
    al_tmp(indd(i)).selectedtowrite = 1;
end

if nargin == 1
    [filename, pathname] = uiputfile({'*.mat'}, 'Save as mat file');
    filename = [pathname '/' filename];
else
    filename = name;
end

if ischar(filename)
    [pathstr, name, ext] = fileparts(filename);
    outstackfile = [pathstr '/' name '.em'];
    outalignfile = [pathstr '/' name '.mat'];
    
    if handles.storage.dims == 2 
        tom_av2_filterstack(handles.storage.stackfilename,al_tmp,'selectedtowrite','==',1,outstackfile,outalignfile);
    elseif handles.storage.dims ==1
        tom_av2_filterstack(stackfile2,al_tmp,'selectedtowrite','==',1,outstackfile,outalignfile);
    else
        al = handles.storage.align;
        lauf=1;
        j = size(al,1);
        for k=1:size(al,2)
            if al_tmp(1,k).selectedtowrite == 1
               for h=1:j
                   Align(h,lauf) = al(h,k);
               end
               lauf=lauf+1; 
            end
        end
        try;Align = rmfield(Align,'selectedtowrite');end
        save(outalignfile,'Align');
    end
end
try;handles.storage.align = rmfield(handles.storage.align,'selectedtowrite');end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display marks                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = disp_mark(handles,number,binning)

if nargin < 3
    binning = str2num(get(handles.kerdensom_binning,'String'));
end

%numrows = str2double(get(handles.kerdensom_input_gridheight,'String'));
numcols = str2double(get(handles.kerdensom_input_gridwidth,'String'));

pos_y = floor(number ./ numcols) * handles.storage.Header.Size(2)./2^binning;
pos_x = (mod(number, numcols)) * handles.storage.Header.Size(1)./2^binning;

width = handles.storage.Header.Size(1)./2^binning;
height = handles.storage.Header.Size(2)./2^binning;
handles.rectangles = [handles.rectangles rectangle('Position',[pos_x, pos_y, width, height],'EdgeColor',[1 0 0],'LineStyle','-','LineWidth',4)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  unset all radio buttons                                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = unset_radiobuttons(handles)

set(handles.disp_gallery,'Value',0);
set(handles.display_codevectors,'Value',0);
set(handles.disp_colormap,'Value',0);
set(handles.display_colormap3d,'Value',0);
set(handles.display_net,'Value',0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  load file                                                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = loadfile(handles)

filename = get(handles.kerdensom_alignmentfile,'String');

try
    s = load(filename);
catch
    if tom_isemfile(filename)
        handles.storage.dims = 1;
    else
        errordlg('File is not readable or does not exist.');
    end
    return;
end

%3d
if isfield(s,'Align')

    handles.storage.dims = 3;
    handles.storage.alignallruns = s.Align;
    
%2d
elseif isfield(s,'align2d')
    
    handles.storage.dims = 2;
    handles.storage.alignallruns = s.align2d;
    
else
    errordlg('Could not determine if file is 2D or 3D alignment file.');
    return;
end

%update runs dropdown field
liststring = '';
for i=1:size(handles.storage.alignallruns,1)
    liststring = strvcat(liststring,num2str(i));
end
set(handles.kerdensom_run,'String',liststring);
set(handles.kerdensom_run,'Enable','on');
set(handles.kerdensom_binning,'Enable','on');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display codvecs                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = disp_codvecs(handles)

axes(handles.kerdensom_browser);
if handles.storage.dims == 2
    tom_dspcub(squeeze(handles.storage.result(handles.showiteration).outputMap));
else
    start = size(handles.storage.result(handles.showiteration).outputMap,3)./2-3;
    stop = size(handles.storage.result(handles.showiteration).outputMap,3)./2+3;
    tom_dspcub(squeeze(sum(handles.storage.result(handles.showiteration).outputMap(:,:,start:stop,:),3)));
end

handles = disp_particlenumbers(handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display average gallery                                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = disp_avggallery(handles)

axes(handles.kerdensom_browser);


if isempty(handles.avgcache) || size(handles.avgcache,2) < handles.showiteration || isempty(handles.avgcache(handles.showiteration).isvalid) || handles.avgcache(handles.showiteration).isvalid == 0
    handles = generate_avgs(handles);
    handles.avgcache(handles.showiteration).avg = handles.avg;
    handles.avgcache(handles.showiteration).isvalid = 1;
else
    %if ~isfield(handles,'avg')
    %    handles.avg
    %end
end

if handles.storage.dims == 2
    tom_dspcub(squeeze(handles.avgcache(handles.showiteration).avg));
else
    start = size(handles.avgcache(handles.showiteration).avg,3)./2-3;
    stop = size(handles.avgcache(handles.showiteration).avg,3)./2+3;
    tom_dspcub(squeeze(sum(handles.avgcache(handles.showiteration).avg(:,:,start:stop,:),3)));
end


handles = disp_particlenumbers(handles,0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display colormap                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = disp_colormap(handles)

%generate colormap
figure;cm=colormap('hsv');close(gcf);

histogram = tom_norm(handles.storage.result(handles.showiteration).histogram,63)+1;

height = str2double(get(handles.kerdensom_input_gridheight,'String'));
width = str2double(get(handles.kerdensom_input_gridwidth,'String'));

handles.storage.colormap = zeros(1,1,height.* width);

for i=1:height.*width
    handles.storage.colormap(1,1,i) = histogram(i);
end
handles.storage.colormap = reshape(handles.storage.colormap,[height,width]);
axes(handles.kerdensom_browser);
imagesc(handles.storage.colormap');colormap('hot');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display colormap 3d                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = disp_colormap3d(handles)

%generate colormap
figure;cm=colormap('hsv');close(gcf);

histogram = tom_norm(handles.storage.result(handles.showiteration).histogram,63)+1;

height = str2double(get(handles.kerdensom_input_gridheight,'String'));
width = str2double(get(handles.kerdensom_input_gridwidth,'String'));

handles.storage.colormap = zeros(1,1,height.* width);

for i=1:height.*width
    handles.storage.colormap(1,1,i) = histogram(i);
end
handles.storage.colormap = reshape(handles.storage.colormap,[height,width]);
axes(handles.kerdensom_browser);
h = bar3(handles.storage.colormap');colormap('hot');

for i = 1:length(h)
    zdata = ones(6*length(h),4);
    k = 1;
    for j = 0:6:(6*length(h)-6)
        zdata(j+1:j+6,:) = handles.storage.colormap(k,i);
        k = k+1;
    end
    set(h(i),'Cdata',zdata)
end
handles.colorbar = colorbar;
shading interp;
for i = 1:length(h)
    zdata = get(h(i),'Zdata');
    set(h(i),'Cdata',zdata)
    set(h,'EdgeColor','k')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display variance map                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = disp_variancemap(handles)

height = str2double(get(handles.kerdensom_input_gridheight,'String'));
width = str2double(get(handles.kerdensom_input_gridwidth,'String'));

vari = handles.storage.result(1).outputMap;
if handles.storage.dims == 2
    numclasses = size(handles.storage.result(1).outputMap,3);
else
    numclasses = size(handles.storage.result(1).outputMap,4);
end

for i=1:numclasses
    m=tom_calc_variance(handles.storage.imstack,handles.storage.result(handles.showiteration).assignVtoX,i);
    if handles.storage.dims == 2
        vari(:,:,i) = tom_norm(reshape(m,size(handles.storage.result(1).outputMap,1),size(handles.storage.result(1).outputMap,2))+100,'phase');
        handles.storage.valnos(i) = mean(mean(vari(:,:,i)));
    else
        vari(:,:,:,i) = tom_norm(reshape(m,size(handles.storage.result(1).outputMap,1),size(handles.storage.result(1).outputMap,2),size(handles.storage.result(1).outputMap,3))+100,'phase');
        handles.storage.valnos(i) = mean(mean(mean(vari(:,:,:,i))));
    end
end

axes(handles.kerdensom_browser);

if handles.storage.dims == 2
    tom_dspcub(tom_norm(tom_norm(vari,1),'3std'));
else
    start = size(vari,3)./2-3;
    stop = size(vari,3)./2+3;

    tom_dspcub(tom_norm(tom_norm(squeeze(sum(vari(:,:,start:stop,:),3)),1),'3std'));
end

handles = disp_particlenumbers(handles);
handles = disp_varnumbers(handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display net                                                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = disp_distmap(handles)

height = str2double(get(handles.kerdensom_input_gridheight,'String'));
width = str2double(get(handles.kerdensom_input_gridwidth,'String'));

try;
    delete(handles.colorbar);
end;
axes(handles.kerdensom_net);
hold off;
reset(handles.kerdensom_net);
cla;
set(gcf,'DoubleBuffer','on');
%set(gca,'XTick',[],'YTick',[],'XLim',[handles.storage.net(handles.showiteration).xmin+handles.storage.net(handles.showiteration).xmin.*.3 handles.storage.net(handles.showiteration).xmax+handles.storage.net(handles.showiteration).xmax.*.3],'YLim',[handles.storage.net(handles.showiteration).ymin+handles.storage.net(handles.showiteration).ymin.*.3 handles.storage.net(handles.showiteration).ymax+handles.storage.net(handles.showiteration).ymax.*.3]);
set(gca,'XTick',[],'YTick',[],'XLim',[-.2 1.2],'YLim',[-.2 1.2]);
tom_plot_grid(handles.storage.net(handles.showiteration).net,[height width], 'rect',gca, 1, handles.storage.result(handles.showiteration).histogram);
handles = generate_miniviewpoints(handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  calculate net                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = calc_net(handles)

height = str2double(get(handles.kerdensom_input_gridheight,'String'));
width = str2double(get(handles.kerdensom_input_gridwidth,'String'));

s = get(handles.edit_origin,'String');
[T,R] = strtok(s,',');
origin(1) = str2num(T);
origin(2) = str2num(R);



contents = get(handles.popupmenu_kerdensom_metric,'String');
metric = contents(get(handles.popupmenu_kerdensom_metric,'Value'));

h = waitbar(0,'Calculating net...');

for i=1:handles.storage.regSteps
    if handles.storage.dims == 3
        map = reshape(handles.storage.result(i).outputMap,handles.storage.Header.Size(1)./2^handles.storage.binning.*handles.storage.Header.Size(2)./2^handles.storage.binning.*handles.storage.Header.Size(3)./2^handles.storage.binning,[]);
        dimensions = [handles.storage.Header.Size(1)./2^handles.storage.binning,handles.storage.Header.Size(2)./2^handles.storage.binning,handles.storage.Header.Size(3)./2^handles.storage.binning];
    elseif handles.storage.dims ==2 
        map = reshape(handles.storage.result(i).outputMap,handles.storage.Header.Size(1)./2^handles.storage.binning.*handles.storage.Header.Size(2)./2^handles.storage.binning,[]);
        dimensions = [handles.storage.Header.Size(1)./2^handles.storage.binning,handles.storage.Header.Size(2)./2^handles.storage.binning];
    else
        map = handles.storage.result(i).outputMap;
    end
        
    handles.storage.net(i).net = tom_codevecs2distmap(map',[height width],'rect',metric,dimensions,origin);
    handles.storage.net(i).net = tom_norm(handles.storage.net(i).net,1);
    handles.storage.net(i).xmin = min(handles.storage.net(i).net(:,1));
    handles.storage.net(i).xmax = max(handles.storage.net(i).net(:,1));
    handles.storage.net(i).ymin = min(handles.storage.net(i).net(:,2));
    handles.storage.net(i).ymax = max(handles.storage.net(i).net(:,2));
    waitbar(i./handles.storage.regSteps,h,'Calculating net...');
end

close(h);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  generate gallery                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = generate_avgs(handles)


binning = str2double(get(handles.kerdensom_binning,'String'));    
if handles.storage.dims == 3
    handles.avg = zeros(size(handles.storage.result(1).outputMap,1).*2^binning,size(handles.storage.result(1).outputMap,2).*2^binning,size(handles.storage.result(1).outputMap,3).*2^binning, length(handles.storage.result(1).histogram));
else
    handles.avg = zeros(size(handles.storage.result(1).outputMap,1).*2^binning,size(handles.storage.result(1).outputMap,2).*2^binning,size(handles.storage.result(1).outputMap,3));
end    

h = waitbar(0,'Generating class averages...');
numparticles = size(handles.storage.align,2);

if handles.storage.dims == 2
    for i=1:numparticles
        part=tom_emreadc(handles.storage.stackfilename,'subregion',[1 1 i],[handles.storage.Header.Size(1)-1 handles.storage.Header.Size(2)-1 0]);
        handles.avg(:,:,handles.storage.align(handles.showiteration,i).Class)=handles.avg(:,:,handles.storage.align(handles.showiteration,i).Class)+part.Value;
        waitbar(i./numparticles,h,['Generating class averages... ( ' num2str(i) ' of ' num2str(numparticles) ' done.)']);
    end
else
    for i=1:length(handles.storage.result(1).histogram)
        lauf = 1;
        for ii=1:numparticles
            if handles.storage.align(handles.showiteration,ii).Class == i
                al(lauf) = handles.storage.align(1,ii); 
                al(lauf).CCC = 1;
                lauf = lauf + 1;
            end
        end
        if exist('al','var') && ~isempty(al(1,1).Tomogram)
           handles.avg(:,:,:,i) = tom_av3_average(al,'sum',0,al(1,end).Class,1);
        end
        
        clear('al');
    end
end
close(h);

if handles.storage.dims == 2
    for i=1:length(handles.storage.result(1).histogram)
        handles.avg(:,:,i) = tom_norm(handles.avg(:,:,i),1);
    end
else
    for i=1:length(handles.storage.result(1).histogram)
        handles.avg(:,:,:,i) = tom_norm(handles.avg(:,:,:,i),1);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display particle numbers                                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = disp_particlenumbers(handles,binning);

if nargin == 1
    binning = str2double(get(handles.kerdensom_binning,'String'));
end

numrows = str2double(get(handles.kerdensom_input_gridheight,'String'));
numcols = str2double(get(handles.kerdensom_input_gridwidth,'String'));
%fontsize = handles.storage.Header.Size(1)./numcols./1.3;
fontsize=16;
ii = 1;

for row=1:numrows
    pos_y = (row-1).*handles.storage.Header.Size(2)./2^binning+handles.storage.Header.Size(2)./2^binning.*0.9;
    for column=1:numcols
        val = handles.storage.result(handles.showiteration).histogram(ii);
        pos_x = (column-1).*handles.storage.Header.Size(1)./2^binning+handles.storage.Header.Size(1)./2^binning.*0.65;
        t = text(pos_x,pos_y,sprintf('%i',val));
        set(t,'Tag','string','Color',[1 1 0],'FontSize',fontsize,'FontWeight','bold');
        try
        if isfield(handles.storage,'thresholds') && ~isempty(handles.storage.thresholds(ii)) && handles.storage.thresholds(ii) ~= 0
            pos_x = (column-1).*handles.storage.Header.Size(1)./2^binning+handles.storage.Header.Size(1)./2^binning.*0.1;
            t = text(pos_x,pos_y,sprintf('%0.2g',handles.storage.thresholds(ii)));
            set(t,'Tag','string','Color',[1 1 0],'FontSize',fontsize,'FontWeight','bold');
        end
        end
        ii = ii + 1;        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display variance numbers                                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = disp_varnumbers(handles,binning);

if nargin == 1
    binning = str2double(get(handles.kerdensom_binning,'String'));
end

numrows = str2double(get(handles.kerdensom_input_gridheight,'String'));
numcols = str2double(get(handles.kerdensom_input_gridwidth,'String'));
%fontsize = handles.storage.Header.Size(1)./numcols./1.3;
fontsize=12;
ii = 1;

for row=1:numrows
    pos_y = (row-1).*handles.storage.Header.Size(2)./2^binning+handles.storage.Header.Size(2)./2^binning.*0.1;
    for column=1:numcols
        val = handles.storage.valnos(ii);
        ii = ii + 1;
        pos_x = (column-1).*handles.storage.Header.Size(1)./2^binning+handles.storage.Header.Size(1)./2^binning.*0.4;
        t = text(pos_x,pos_y,sprintf('%0.3g',val));
        set(t,'Tag','string','Color',[1 1 0],'FontSize',fontsize,'FontWeight','bold');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  initialize empty avg cache                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = init_avgcache(handles)

handles.avgcache = struct();


for i=1:size(handles.storage.result,2)
    handles.avgcache(i).avg = [];
    handles.avgcache(i).isvalid = 0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  enable / disable net axes                                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = show_netaxes(handles,val)

set(handles.kerdensom_net,'Visible',val,'XTick',[],'YTick',[]);
set(handles.axes_kerdensom_right1,'Visible',val,'XTick',[],'YTick',[]);axes(handles.axes_kerdensom_right1);cla;set(gca,'Tag','axes_kerdensom_right1');
set(handles.axes_kerdensom_right2,'Visible',val,'XTick',[],'YTick',[]);axes(handles.axes_kerdensom_right2);cla;set(gca,'Tag','axes_kerdensom_right2');
set(handles.axes_kerdensom_right3,'Visible',val,'XTick',[],'YTick',[]);axes(handles.axes_kerdensom_right3);cla;set(gca,'Tag','axes_kerdensom_right3');
set(handles.axes_kerdensom_right4,'Visible',val,'XTick',[],'YTick',[]);axes(handles.axes_kerdensom_right4);cla;set(gca,'Tag','axes_kerdensom_right4');
set(handles.axes_kerdensom_right5,'Visible',val,'XTick',[],'YTick',[]);axes(handles.axes_kerdensom_right5);cla;set(gca,'Tag','axes_kerdensom_right5');
set(handles.axes_kerdensom_right6,'Visible',val,'XTick',[],'YTick',[]);axes(handles.axes_kerdensom_right6);cla;set(gca,'Tag','axes_kerdensom_right6');
set(handles.axes_kerdensom_left1,'Visible',val,'XTick',[],'YTick',[]);axes(handles.axes_kerdensom_left1);cla;set(gca,'Tag','axes_kerdensom_left1');
set(handles.axes_kerdensom_left2,'Visible',val,'XTick',[],'YTick',[]);axes(handles.axes_kerdensom_left2);cla;set(gca,'Tag','axes_kerdensom_left2');
set(handles.axes_kerdensom_left3,'Visible',val,'XTick',[],'YTick',[]);axes(handles.axes_kerdensom_left3);cla;set(gca,'Tag','axes_kerdensom_left3');
set(handles.axes_kerdensom_left4,'Visible',val,'XTick',[],'YTick',[]);axes(handles.axes_kerdensom_left4);cla;set(gca,'Tag','axes_kerdensom_left4');
set(handles.axes_kerdensom_left5,'Visible',val,'XTick',[],'YTick',[]);axes(handles.axes_kerdensom_left5);cla;set(gca,'Tag','axes_kerdensom_left5');
set(handles.axes_kerdensom_left6,'Visible',val,'XTick',[],'YTick',[]);axes(handles.axes_kerdensom_left6);cla;set(gca,'Tag','axes_kerdensom_left6');
set(handles.axes_kerdensom_top1,'Visible',val,'XTick',[],'YTick',[]);axes(handles.axes_kerdensom_top1);cla;set(gca,'Tag','axes_kerdensom_top1');
set(handles.axes_kerdensom_top2,'Visible',val,'XTick',[],'YTick',[]);axes(handles.axes_kerdensom_top2);cla;set(gca,'Tag','axes_kerdensom_top2');
set(handles.axes_kerdensom_top3,'Visible',val,'XTick',[],'YTick',[]);axes(handles.axes_kerdensom_top3);cla;set(gca,'Tag','axes_kerdensom_top3');
set(handles.axes_kerdensom_top4,'Visible',val,'XTick',[],'YTick',[]);axes(handles.axes_kerdensom_top4);cla;set(gca,'Tag','axes_kerdensom_top4');
set(handles.axes_kerdensom_top5,'Visible',val,'XTick',[],'YTick',[]);axes(handles.axes_kerdensom_top5);cla;set(gca,'Tag','axes_kerdensom_top5');
set(handles.axes_kerdensom_top6,'Visible',val,'XTick',[],'YTick',[]);axes(handles.axes_kerdensom_top6);cla;set(gca,'Tag','axes_kerdensom_top6');
set(handles.axes_kerdensom_bottom1,'Visible',val,'XTick',[],'YTick',[]);axes(handles.axes_kerdensom_bottom1);cla;set(gca,'Tag','axes_kerdensom_bottom1');
set(handles.axes_kerdensom_bottom2,'Visible',val,'XTick',[],'YTick',[]);axes(handles.axes_kerdensom_bottom2);cla;set(gca,'Tag','axes_kerdensom_bottom2');
set(handles.axes_kerdensom_bottom3,'Visible',val,'XTick',[],'YTick',[]);axes(handles.axes_kerdensom_bottom3);cla;set(gca,'Tag','axes_kerdensom_bottom3');
set(handles.axes_kerdensom_bottom4,'Visible',val,'XTick',[],'YTick',[]);axes(handles.axes_kerdensom_bottom4);cla;set(gca,'Tag','axes_kerdensom_bottom4');
set(handles.axes_kerdensom_bottom5,'Visible',val,'XTick',[],'YTick',[]);axes(handles.axes_kerdensom_bottom5);cla;set(gca,'Tag','axes_kerdensom_bottom5');
set(handles.axes_kerdensom_bottom6,'Visible',val,'XTick',[],'YTick',[]);axes(handles.axes_kerdensom_bottom6);cla;set(gca,'Tag','axes_kerdensom_bottom6');

if isequal(val,'on')
    axes(handles.kerdensom_browser);
    cla;
    reset(handles.kerdensom_browser);
    set(handles.kerdensom_browser,'Visible','off','XTick',[],'YTick',[]);
else
    axes(handles.kerdensom_net);
    cla;
    reset(handles.kerdensom_net);
    set(handles.kerdensom_net,'Visible','off');
    set(handles.kerdensom_browser,'Visible','on','XTick',[],'YTick',[]);
end
drawnow;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  generate miniviewpoints                                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = generate_miniviewpoints(handles)

handles.miniviewpoints = [];  
handles.xlim = get(handles.kerdensom_net,'Xlim');
handles.ylim = get(handles.kerdensom_net,'Ylim');

xabs = (handles.xlim(2)-handles.xlim(1))./7;
yabs = (handles.ylim(2)-handles.ylim(1))./7;

for i=1:6
    handles.xpoints(i) = handles.xlim(1)+(i).*xabs;
    handles.ypoints(7-i) = handles.ylim(1)+(i).*yabs;
end

for i=1:6
    handles.miniviewpoints = cat(1,handles.miniviewpoints,[handles.xlim(1),handles.ypoints(i)]);
    handles.miniviewpoints = cat(1,handles.miniviewpoints,[handles.xlim(2),handles.ypoints(i)]);
    handles.miniviewpoints = cat(1,handles.miniviewpoints,[handles.xpoints(i),handles.ylim(1)]);
    handles.miniviewpoints = cat(1,handles.miniviewpoints,[handles.xpoints(i),handles.ylim(2)]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  plot miniviews                                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = plot_lines(handles)

for i=1:size(handles.lines,2)

    %draw line from node to miniview
    axes(handles.kerdensom_net);
    linestart = handles.storage.net(handles.showiteration).net(handles.lines(i).class,:);
    line([linestart(1) handles.lines(i).endpoint(1)],[linestart(2) handles.lines(i).endpoint(2)],'Color',[0 0 0],'LineWidth',2);

    
    %draw miniview
    axes(handles.(handles.lines(i).axes));
    if isequal(handles.lines(i).type,'avg')
        if ~isfield(handles,'avgcache') || isempty(handles.avgcache) || size(handles.avgcache,2) < handles.showiteration || isempty(handles.avgcache(handles.showiteration).isvalid) || handles.avgcache(handles.showiteration).isvalid == 0
            handles = generate_avgs(handles);
            handles.avgcache(handles.showiteration).avg = handles.avg;
            handles.avgcache(handles.showiteration).isvalid = 1;
        end
        if handles.storage.dims == 2
            im = handles.avgcache(handles.showiteration).avg(:,:,handles.lines(i).class);
        else
            start = size(handles.avgcache(handles.showiteration).avg,3)./2-5;
            stop = size(handles.avgcache(handles.showiteration).avg,3)./2+5;
            im = sum(handles.avgcache(handles.showiteration).avg(:,:,start:stop,handles.lines(i).class),3);
        end
    else
        if handles.storage.dims == 2
            im = handles.storage.result(handles.showiteration).outputMap(:,:,handles.lines(i).class);
        else
            start = size(handles.storage.result(handles.showiteration).outputMap,3)./2-5;
            stop = size(handles.storage.result(handles.showiteration).outputMap,3)./2+5;
            im = sum(handles.storage.result(handles.showiteration).outputMap(:,:,start:stop,handles.lines(i).class),3);
        end
    end

    tom_imagesc(im,'noinfo');axis off;
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Unused Callbacks/Create functions                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function kerdensom_alignmentfile_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function kerdensom_input_som_Nsteps_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function kerdensom_input_som_Nsteps_Callback(hObject, eventdata, handles)
function kerdensom_run_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function kerdensom_binning_Callback(hObject, eventdata, handles)
function kerdensom_binning_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function kerdensom_input_gridheight_Callback(hObject, eventdata, handles)
function kerdensom_input_gridheight_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function kerdensom_input_gridwidth_Callback(hObject, eventdata, handles)
function kerdensom_input_gridwidth_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function kerdensom_input_initialreg_Callback(hObject, eventdata, handles)
function kerdensom_input_initialreg_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function kerdensom_input_finalreg_Callback(hObject, eventdata, handles)
function kerdensom_input_finalreg_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function kerdensom_input_regsteps_Callback(hObject, eventdata, handles)
function kerdensom_input_regsteps_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function kerdensom_run_Callback(hObject, eventdata, handles)
function checkbox_inverse_Callback(hObject, eventdata, handles)
function popupmenu_kerdensom_iteration_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_kerdensom_movietime_Callback(hObject, eventdata, handles)
function edit_kerdensom_movietime_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function popupmenu_algo_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function popupmenu_algo_Callback(hObject, eventdata, handles)
function popupmenu4_Callback(hObject, eventdata, handles)
function popupmenu4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_kerdensom_host_Callback(hObject, eventdata, handles)
function edit_kerdensom_host_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function checkbox_kerdensom_saveavgcache_Callback(hObject, eventdata, handles)
function edit_movie_start_Callback(hObject, eventdata, handles)
function edit_movie_start_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_movie_stop_Callback(hObject, eventdata, handles)
function edit_movie_stop_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function checkbox_kerdensom_voltex_Callback(hObject, eventdata, handles)
function popupmenu_kerdensom_metric_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_origin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_mass_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_mass_Callback(hObject, eventdata, handles)
function edit_pixelsize_Callback(hObject, eventdata, handles)
function edit_pixelsize_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function checkbox_net_classavgs_Callback(hObject, eventdata, handles)


% --- Executes on button press in wedge_it.
function wedge_it_Callback(hObject, eventdata, handles)
% hObject    handle to wedge_it (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of wedge_it


