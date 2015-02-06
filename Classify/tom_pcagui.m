
function varargout = tom_pcagui(varargin)
%TOM_PCAGUI is a GUI for principal component analysis and clustering
%
%GUI for principal component analysis and clustering
%
%   varargout = tom_pcagui(varargin)
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
%   ... = tom_pcagui(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by AK 06/08/06
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

gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_pcagui_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_pcagui_OutputFcn, ...
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
function tom_pcagui_OpeningFcn(hObject, eventdata, handles, varargin)

if nargin > 3
    set(handles.input_pcagui_inputfile,'String',varargin{1});
    handles = loadfile(handles);
end

%delete all ticks on all axes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(handles.axes_pcagui_mainview,'XTick',[],'YTick',[]);

set(handles.axes_pcagui_extrema_ul,'XTick',[],'YTick',[]);
set(handles.axes_pcagui_extrema_ur,'XTick',[],'YTick',[]);
set(handles.axes_pcagui_extrema_ll,'XTick',[],'YTick',[]);
set(handles.axes_pcagui_extrema_lr,'XTick',[],'YTick',[]);

set(handles.axes_pcagui_right1,'XTick',[],'YTick',[]);
set(handles.axes_pcagui_right2,'XTick',[],'YTick',[]);
set(handles.axes_pcagui_right3,'XTick',[],'YTick',[]);
set(handles.axes_pcagui_right4,'XTick',[],'YTick',[]);
set(handles.axes_pcagui_right5,'XTick',[],'YTick',[]);
set(handles.axes_pcagui_right6,'XTick',[],'YTick',[]);
set(handles.axes_pcagui_right7,'XTick',[],'YTick',[]);
set(handles.axes_pcagui_left1,'XTick',[],'YTick',[]);
set(handles.axes_pcagui_left2,'XTick',[],'YTick',[]);
set(handles.axes_pcagui_left3,'XTick',[],'YTick',[]);
set(handles.axes_pcagui_left4,'XTick',[],'YTick',[]);
set(handles.axes_pcagui_left5,'XTick',[],'YTick',[]);
set(handles.axes_pcagui_left6,'XTick',[],'YTick',[]);
set(handles.axes_pcagui_left7,'XTick',[],'YTick',[]);
set(handles.axes_pcagui_top1,'XTick',[],'YTick',[]);
set(handles.axes_pcagui_top2,'XTick',[],'YTick',[]);
set(handles.axes_pcagui_top3,'XTick',[],'YTick',[]);
set(handles.axes_pcagui_top4,'XTick',[],'YTick',[]);
set(handles.axes_pcagui_top5,'XTick',[],'YTick',[]);
set(handles.axes_pcagui_top6,'XTick',[],'YTick',[]);
set(handles.axes_pcagui_top7,'XTick',[],'YTick',[]);
set(handles.axes_pcagui_bottom1,'XTick',[],'YTick',[]);
set(handles.axes_pcagui_bottom2,'XTick',[],'YTick',[]);
set(handles.axes_pcagui_bottom3,'XTick',[],'YTick',[]);
set(handles.axes_pcagui_bottom4,'XTick',[],'YTick',[]);
set(handles.axes_pcagui_bottom5,'XTick',[],'YTick',[]);
set(handles.axes_pcagui_bottom6,'XTick',[],'YTick',[]);
set(handles.axes_pcagui_bottom7,'XTick',[],'YTick',[]);

handles.cluster = struct();
handles.storage.eigenvalues = [];
handles.storage.onedstackfile = '';
handles.storage.onedalignfile = '';

handles.labels = [];
handles.points = [];

handles.line.dummy = [];
handles.point.dummy = [];

handles.dispclassno = '';

% Choose default command line output for tom_pcagui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes tom_pcagui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Output Function                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = tom_pcagui_OutputFcn(hObject, eventdata, handles) 

% Get default command line output from handles structure
varargout{1} = handles.output;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  input file                                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_pcagui_inputfile_Callback(hObject, eventdata, handles)

handles = loadfile(handles);

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  browse input file                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_pcagui_browseinputfile_Callback(hObject, eventdata, handles)

[filename, pathname] = uigetfile({'*.mat'}, 'Pick an alignment file');
if ischar(filename)
    set(handles.input_pcagui_inputfile,'String',[pathname '/' filename]);
    handles = loadfile(handles);
end

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  generate input file                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_pcagui_generate_inputfile_Callback(hObject, eventdata, handles)

filename = tom_av3_createalign();

if ischar(filename)
   set(handles.input_pcagui_inputfile,'String',filename); 
   handles = loadfile(handles);
end

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  load input file                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_pcagui_loadfile_Callback(hObject, eventdata, handles)

contents = get(handles.popupmenu_pcagui_run,'String');
run = str2num(contents(get(handles.popupmenu_pcagui_run,'Value')));

contents = get(handles.popupmenu_pcagui_space,'String');
space = contents(get(handles.popupmenu_pcagui_space,'Value'));

set(handles.popupmenu_pcagui_run,'Enable','off');
set(handles.input_pcagui_binning,'Enable','off');
binning = str2num(get(handles.input_pcagui_binning,'String'));

h = waitbar(0,'Reading files...');

normflag = get(handles.checkbox_pcagui_phasenorm,'Value');

%1d
if handles.storage.dims == 1
    set(handles.panel3d,'Visible','off');
    handles.storage.stackfilename = get(handles.input_pcagui_inputfile,'String');
    header = tom_reademheader(handles.storage.stackfilename);
    handles.storage.Header = header.Header;
    handles.storage.imstack = tom_emreadc(handles.storage.stackfilename);
    handles.storage.imstack = double(handles.storage.imstack.Value);
    set(handles.panel3d,'Visible','off');
    set(handles.panel1d,'Visible','on');
%2d
elseif handles.storage.dims == 2
    handles.storage.align = handles.storage.alignallruns(run,:);   
    [pathstr, name, ext] = fileparts(get(handles.input_pcagui_inputfile,'String'));
    handles.storage.stackfilename = [pathstr '/' name '.em'];
    header = tom_reademheader(handles.storage.stackfilename);
    handles.storage.Header = header.Header;
    if strcmp(space,'pixel space')
        if normflag == 1
            in_struct.mask.types = {'sphere','rectangle'};
            maskstruct = tom_filtergui('mask',in_struct);
            handles.storage.imstack = tom_reshape_stack(handles.storage.stackfilename,'',str2num(get(handles.input_pcagui_binning,'String')),'',normflag,maskstruct);
        else
            handles.storage.imstack = tom_reshape_stack(handles.storage.stackfilename,'',str2num(get(handles.input_pcagui_binning,'String')));
        end
    end
    set(handles.panel3d,'Visible','off');
    set(handles.panel1d,'Visible','off');
    
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
    
    if (tom_isemfile(filecell{1}))
        header = tom_reademheader(filecell{1});
    else
       tmppp=tom_spiderread(filecell{1});
       header=tom_emheader(tmppp.Value);
    end;
    
    
    handles.storage.Header = header.Header;
    if strcmp(space,'pixel space')
        if normflag == 1
            in_struct.mask.types = {'sphere3d','cylinder3d'};
            in_struct.rotmask.types = {'sphere3d','cylinder3d'};
            maskstruct = tom_filtergui('mask',in_struct);
        else
            in_struct.rotmask.types = {'sphere3d','cylinder3d'};
            maskstruct = tom_filtergui('mask',in_struct);
        end
        maskstruct.Apply=2;
        handles.storage.imstack = tom_reshape_stack(filecell,handles.storage.align,str2num(get(handles.input_pcagui_binning,'String')),'',normflag,maskstruct);
    end
    set(handles.panel3d,'Visible','on');
    set(handles.panel1d,'Visible','off');
    
else
    errordlg('Alignment file not valid!');
    close(h);
    return;
end
if strcmp(space,'pixel space')
    waitbar(0.5,h,'Removing mean...');
    [handles.storage.imstack handles.storage.mean]=tom_rm_mean(handles.storage.imstack);
end

close(h);

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  calculate pca                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_pcagui_calcpca_Callback(hObject, eventdata, handles)

contents = get(handles.popupmenu_pcagui_run,'String');
run = str2num(contents(get(handles.popupmenu_pcagui_run,'Value')));
binning = str2num(get(handles.input_pcagui_binning,'String'));

contents = get(handles.popupmenu_pcagui_space,'String');
space = contents(get(handles.popupmenu_pcagui_space,'Value'));

if strcmp(space,'particle space')
    if str2num(get(handles.input_pcagui_numeigenvectors,'String')) >= size(handles.storage.align,2)
        error('Number of eigenvectors cannot exceed number of particles');
    end
end

if strcmp(space,'pixel space')
    h = waitbar(0,'Calculating pca...');
    [handles.storage.scores,handles.storage.coefs, eigenvalues]=tom_calc_pca(handles.storage.imstack,str2num(get(handles.input_pcagui_numeigenvectors,'String')),'pca');
    close(h);
else
    if handles.storage.dims == 2
        in_struct.mask.types = {'sphere','rectangle'};
        maskstruct = tom_filtergui('mask',in_struct);
    else
        in_struct.mask.types = {'sphere3d','cylinder3d'};
        maskstruct = tom_filtergui('mask',in_struct);
    end
    maskstruct.mask.Value = handles.storage.Header.Size./2^binning;
    if get(handles.checkbox_pcagui_cpca,'Value') == 1
        [handles.storage.scores,handles.storage.coefs, eigenvalues]=tom_calc_pca('',str2num(get(handles.input_pcagui_numeigenvectors,'String')),'cpca',handles.storage.align,1,0,30,maskstruct);
    else
        [handles.storage.scores,handles.storage.coefs, eigenvalues]=tom_calc_pca('',str2num(get(handles.input_pcagui_numeigenvectors,'String')),'pca',handles.storage.align,1,0,30,maskstruct);
    end
end


for i=1:size(eigenvalues,1)
    handles.storage.eigenvalues = [handles.storage.eigenvalues;eigenvalues(i,i)];
end

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  save pca                                                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_pcagui_savecalc_Callback(hObject, eventdata, handles)

[filename, pathname] = uiputfile('*.mat','Save pca calculation');

if ischar(filename)
    pcacalc = struct();
    pcacalc.scores = handles.storage.scores;
    pcacalc.coefs = handles.storage.coefs;
    pcacalc.filename = get(handles.input_pcagui_inputfile,'String');
    save([pathname '/' filename],'pcacalc');
    disp('PCA Calculation saved.');
end

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  load pca                                                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_pcagui_pcacalcload_Callback(hObject, eventdata, handles)

[filename, pathname] = uigetfile('*.mat','Load pca calculation');

if ischar(filename)
    p = load([pathname '/' filename]);
    if ~isfield(p,'pcacalc')
        errordlg('This is not a saved pca calculation!');
        return;
    end
    if strcmp(p.pcacalc.filename,get(handles.input_pcagui_inputfile,'String')) == 0
        button = questdlg('PCA Calculation was performed on other data file, continue loading?','Warning','yes','no','no');
        if strcmp(button,'no')
            disp('Cancelled.');
            return;
        end
    end
    handles.storage.scores = p.pcacalc.scores;
    handles.storage.coefs = p.pcacalc.coefs;
    disp('PCA calculation loaded.');
end

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display radio buttons                                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radiobutton_pcagui_dispeigenimages_Callback(hObject, eventdata, handles)


if get(hObject,'Value') == 1
    set(handles.radiobutton_pcagui_dispparticles,'Value',0);
end
guidata(hObject, handles);

function radiobutton_pcagui_dispparticles_Callback(hObject, eventdata, handles)

if get(hObject,'Value') == 1
    set(handles.radiobutton_pcagui_dispeigenimages,'Value',0);
end
guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  pick particles                                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_pcagui_pickparticles_Callback(hObject, eventdata, handles)

xlim = get(handles.axes_pcagui_mainview,'Xlim');
ylim = get(handles.axes_pcagui_mainview,'Ylim');

xabs = (xlim(2)-xlim(1))./8;
yabs = (ylim(2)-ylim(1))./8;

for i=1:7
    xpoints(i) = xlim(1)+(i).*xabs;
    ypoints(8-i) = ylim(1)+(i).*yabs;
end

if get(hObject,'Value') == 1
    mousebutton = 'normal';
    while strcmp(mousebutton,'normal') == 1
        set(handles.output,'Pointer','crosshair');
        axes1 = '';
        while strcmp(axes1,'axes_pcagui_mainview') == 0
            waitforbuttonpress;
            point1 = get(gca,'CurrentPoint');
            axes1 = get(gca,'Tag');
        end
        
        axes2 = '';
        while ~(~isempty(strfind(axes2,'axes_pcagui_right')) == 1 || ~isempty(strfind(axes2,'axes_pcagui_left')) == 1 || ~isempty(strfind(axes2,'axes_pcagui_top')) == 1 || ~isempty(strfind(axes2,'axes_pcagui_bottom')) == 1)
            waitforbuttonpress;
            point2 = get(gca,'CurrentPoint');
            axes2 = get(gca,'Tag');
            mousebutton = get(gcf,'SelectionType');
        end

        %display eigen image
        if get(handles.radiobutton_pcagui_dispeigenimages,'Value') == 1
            [handles, im_out] = get_eigenimage(handles,[point1(1),point1(3)]);
            pointtype = 'eigen';
            
        %display orig image
        else
            [handles, im_out, coords_out] = get_origimage(handles,[point1(1),point1(3)]);
            point1(1) = coords_out(1);
            point1(3) = coords_out(2);
            pointtype = 'orig';
        end
        
        %calculate line coordinates
        linestart = [point1(1),point1(3)];
        if ~isempty(strfind(axes2,'axes_pcagui_right')) == 1
            linestop = [xlim(2) ypoints(str2num(axes2(end)))];
        elseif ~isempty(strfind(axes2,'axes_pcagui_left')) == 1
            linestop = [xlim(1) ypoints(str2num(axes2(end)))];
        elseif ~isempty(strfind(axes2,'axes_pcagui_top')) == 1            
            linestop = [xpoints(str2num(axes2(end))) ylim(2)];
        elseif ~isempty(strfind(axes2,'axes_pcagui_bottom')) == 1
            linestop = [xpoints(str2num(axes2(end))) ylim(1)];
        end
        
        %draw image
        axes(handles.(axes2));
        if handles.storage.dims > 1
            imagesc(im_out');colormap gray;%axis ij;
        else
            plot(im_out);
        end
        set(handles.(axes2),'XTick',[],'YTick',[],'Tag',axes2,'Userdata',{pointtype [point1(1) point1(3)]});
        [b c d] = strread(axes2,'%s %s %s','delimiter','_');
        set(handles.(['text_pcagui_' d{1}]),'String',pointtype);
        
        %render line
        axes(handles.axes_pcagui_mainview);
        hold on;
        if isfield(handles.line, axes2)
            try;delete(handles.line.(axes2));end;
            try;delete(handles.point.(axes2));end;
        end
        handles.point.(axes2) = plot(point1(1),point1(3),'Color',[0 0 0],'Marker','o'); 
        handles.line.(axes2) = line([linestart(1) linestop(1)],[linestart(2) linestop(2)],'Color',[0 0 0],'LineWidth',1);

        hold off;
        set(handles.axes_pcagui_mainview,'Tag','axes_pcagui_mainview');

    end
    set(handles.output,'Pointer','arrow');
else
    
end

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display particle 3D                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_pcagui_show3d_Callback(hObject, eventdata, handles)

axes2 = '';
while ~(~isempty(strfind(axes2,'axes_pcagui_right')) == 1 || ~isempty(strfind(axes2,'axes_pcagui_left')) == 1 || ~isempty(strfind(axes2,'axes_pcagui_top')) == 1 || ~isempty(strfind(axes2,'axes_pcagui_bottom')) == 1)
    waitforbuttonpress;
    axes2 = get(gca,'Tag');
    mousebutton = get(gcf,'SelectionType');
end

userdata = get(handles.(axes2), 'Userdata');
try
    pointtype = userdata{1};
catch
    return;
end

%eigen image
if isequal(pointtype,'eigen')
    [handles, im_out] = get_eigenimage(handles,userdata{2},1);
%orig image
elseif isequal(pointtype,'orig')
    [handles, im_out] = get_origimage(handles,userdata{2},1);
%average
else
    if handles.storage.dims==2
        im_out=zeros(handles.storage.Header.Size(1),handles.storage.Header.Size(2));
    end

    indd = handles.storage.selectedpoints;

    for i=1:length(indd)
        if (handles.storage.dims == 3)
            al(i) = handles.storage.align(indd(i)); 
        else
            part=tom_emreadc(handles.storage.stackfilename,'subregion',[1 1 indd(i)],[handles.storage.Header.Size(1)-1 handles.storage.Header.Size(1)-1 0]);
            part=part.Value;
            im_out=im_out+part;
        end
    
    end

    if handles.storage.dims == 3
        im_out = tom_av3_average(al,'sum',0,0,1);
    end
end



if get(handles.radiobutton_pcagui_volxyz,'Value') == 1
    tom_volxyz(im_out);
else
    figure;tom_dspcub(im_out);
end

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display particle 3D in dspcub                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radiobutton_pcagui_dspcub_Callback(hObject, eventdata, handles)

if get(hObject,'Value') == 1
    set(handles.radiobutton_pcagui_volxyz,'Value',0);
end

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display particle 3D  in volxyz                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radiobutton_pcagui_volxyz_Callback(hObject, eventdata, handles)

if get(hObject,'Value') == 1
    set(handles.radiobutton_pcagui_dspcub,'Value',0);
end

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  plot pca                                                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_pcagui_disppca_Callback(hObject, eventdata, handles)

delete(get(handles.axes_pcagui_right1,'Children'));
delete(get(handles.axes_pcagui_right2,'Children'));
delete(get(handles.axes_pcagui_right3,'Children'));
delete(get(handles.axes_pcagui_right4,'Children'));
delete(get(handles.axes_pcagui_right5,'Children'));
delete(get(handles.axes_pcagui_right6,'Children'));
delete(get(handles.axes_pcagui_right7,'Children'));

delete(get(handles.axes_pcagui_left1,'Children'));
delete(get(handles.axes_pcagui_left2,'Children'));
delete(get(handles.axes_pcagui_left3,'Children'));
delete(get(handles.axes_pcagui_left4,'Children'));
delete(get(handles.axes_pcagui_left5,'Children'));
delete(get(handles.axes_pcagui_left6,'Children'));
delete(get(handles.axes_pcagui_left7,'Children'));

delete(get(handles.axes_pcagui_top1,'Children'));
delete(get(handles.axes_pcagui_top2,'Children'));
delete(get(handles.axes_pcagui_top3,'Children'));
delete(get(handles.axes_pcagui_top4,'Children'));
delete(get(handles.axes_pcagui_top5,'Children'));
delete(get(handles.axes_pcagui_top6,'Children'));
delete(get(handles.axes_pcagui_top7,'Children'));

delete(get(handles.axes_pcagui_bottom1,'Children'));
delete(get(handles.axes_pcagui_bottom2,'Children'));
delete(get(handles.axes_pcagui_bottom3,'Children'));
delete(get(handles.axes_pcagui_bottom4,'Children'));
delete(get(handles.axes_pcagui_bottom5,'Children'));
delete(get(handles.axes_pcagui_bottom6,'Children'));
delete(get(handles.axes_pcagui_bottom7,'Children'));

set(handles.axes_pcagui_right1,'Userdata','');
set(handles.axes_pcagui_right2,'Userdata','');
set(handles.axes_pcagui_right3,'Userdata','');
set(handles.axes_pcagui_right4,'Userdata','');
set(handles.axes_pcagui_right5,'Userdata','');
set(handles.axes_pcagui_right6,'Userdata','');
set(handles.axes_pcagui_right7,'Userdata','');

set(handles.axes_pcagui_left1,'Userdata','');
set(handles.axes_pcagui_left2,'Userdata','');
set(handles.axes_pcagui_left3,'Userdata','');
set(handles.axes_pcagui_left4,'Userdata','');
set(handles.axes_pcagui_left5,'Userdata','');
set(handles.axes_pcagui_left6,'Userdata','');
set(handles.axes_pcagui_left7,'Userdata','');

set(handles.axes_pcagui_top1,'Userdata','');
set(handles.axes_pcagui_top2,'Userdata','');
set(handles.axes_pcagui_top3,'Userdata','');
set(handles.axes_pcagui_top4,'Userdata','');
set(handles.axes_pcagui_top5,'Userdata','');
set(handles.axes_pcagui_top6,'Userdata','');
set(handles.axes_pcagui_top7,'Userdata','');

set(handles.axes_pcagui_bottom1,'Userdata','');
set(handles.axes_pcagui_bottom2,'Userdata','');
set(handles.axes_pcagui_bottom3,'Userdata','');
set(handles.axes_pcagui_bottom4,'Userdata','');
set(handles.axes_pcagui_bottom5,'Userdata','');
set(handles.axes_pcagui_bottom6,'Userdata','');
set(handles.axes_pcagui_bottom7,'Userdata','');

set(handles.text_pcagui_right1,'String','');
set(handles.text_pcagui_right2,'String','');
set(handles.text_pcagui_right3,'String','');
set(handles.text_pcagui_right4,'String','');
set(handles.text_pcagui_right5,'String','');
set(handles.text_pcagui_right6,'String','');
set(handles.text_pcagui_right7,'String','');

set(handles.text_pcagui_left1,'String','');
set(handles.text_pcagui_left2,'String','');
set(handles.text_pcagui_left3,'String','');
set(handles.text_pcagui_left4,'String','');
set(handles.text_pcagui_left5,'String','');
set(handles.text_pcagui_left6,'String','');
set(handles.text_pcagui_left7,'String','');

set(handles.text_pcagui_top1,'String','');
set(handles.text_pcagui_top2,'String','');
set(handles.text_pcagui_top3,'String','');
set(handles.text_pcagui_top4,'String','');
set(handles.text_pcagui_top5,'String','');
set(handles.text_pcagui_top6,'String','');
set(handles.text_pcagui_top7,'String','');

set(handles.text_pcagui_bottom1,'String','');
set(handles.text_pcagui_bottom2,'String','');
set(handles.text_pcagui_bottom3,'String','');
set(handles.text_pcagui_bottom4,'String','');
set(handles.text_pcagui_bottom5,'String','');
set(handles.text_pcagui_bottom6,'String','');
set(handles.text_pcagui_bottom7,'String','');

handles = rmfield(handles,'line');
handles = rmfield(handles,'point');
handles.line.dummy = [];
handles.point.dummy = [];

axes(handles.axes_pcagui_mainview);
index1=str2num(get(handles.input_pcagui_eigenvector1,'String'));
index2=str2num(get(handles.input_pcagui_eigenvector2,'String'));
index3=str2num(get(handles.input_pcagui_eigenvector3,'String'));

%2 eigenvectors
if isempty(index3)
    handles.points = plot(handles.storage.scores(index1,:),handles.storage.scores(index2,:),'r+');
    set(handles.axes_pcagui_mainview,'XTick',[],'YTick',[],'Tag','axes_pcagui_mainview');

%create extrema eigen images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlim = get(handles.axes_pcagui_mainview,'XLim');
ylim = get(handles.axes_pcagui_mainview,'YLim');

[handles, im_out] = get_eigenimage(handles,[xlim(1) ylim(2)]);
axes(handles.axes_pcagui_extrema_ul);
if handles.storage.dims > 1
    imagesc(im_out');axis ij;colormap gray;
else
    plot(im_out);
end
set(handles.axes_pcagui_extrema_ul,'XTick',[],'YTick',[],'Tag','axes_pcagui_extrema_ul');

[handles, im_out] = get_eigenimage(handles,[xlim(2) ylim(2)]);
axes(handles.axes_pcagui_extrema_ur);
if handles.storage.dims > 1
    imagesc(im_out');axis ij;colormap gray;
else
    plot(im_out);
end

set(handles.axes_pcagui_extrema_ur,'XTick',[],'YTick',[],'Tag','axes_pcagui_extrema_ur');

[handles, im_out] = get_eigenimage(handles,[xlim(1) ylim(1)]);
axes(handles.axes_pcagui_extrema_ll);
if handles.storage.dims > 1
    imagesc(im_out');axis ij;colormap gray;
else
    plot(im_out);
end

set(handles.axes_pcagui_extrema_ll,'XTick',[],'YTick',[],'Tag','axes_pcagui_extrema_ll');

[handles, im_out] = get_eigenimage(handles,[xlim(2) ylim(2)]);
axes(handles.axes_pcagui_extrema_lr);
if handles.storage.dims > 1
    imagesc(im_out');axis ij;colormap gray;
else
    plot(im_out);
end
set(handles.axes_pcagui_extrema_lr,'XTick',[],'YTick',[],'Tag','axes_pcagui_extrema_lr');

else
%3 eigenvectors
    figure;
    handles.points = plot3(handles.storage.scores(index1,:),handles.storage.scores(index2,:),handles.storage.scores(index3,:),'r+');grid on;
    set(gca,'Tag','3dcluster');
end

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  plot eigenimage gallery                                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_pcagui_eigenvecgallery_Callback(hObject, eventdata, handles)

[handles, im_out] = get_eigenvectorgallery(handles);

if handles.storage.dims == 2
    figure;tom_dspcub(im_out);
elseif handles.storage.dims == 3
    figure;tom_imagesc(im_out);
end

figure;
bar(handles.storage.eigenvalues);

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  save eigenimage gallery                                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_pcagui_savegallery_Callback(hObject, eventdata, handles)

[handles, im_out] = get_eigenvectorgallery(handles);

tom_emwrite(im_out);

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  select in plot                                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_pcagui_selectinplot_Callback(hObject, eventdata, handles)

try;delete(handles.selectedpointmarkers);end;
handles.selectedpointmarkers = [];

eigvecs(1) = str2num(get(handles.input_pcagui_eigenvector1,'String'));
eigvecs(2) = str2num(get(handles.input_pcagui_eigenvector2,'String'));

axes(handles.axes_pcagui_mainview);

[xi,yi]=getline(handles.axes_pcagui_mainview, 'closed');

% get coordinates
for i=1:size(eigvecs,2)
    index=eigvecs(i);
    tmp(:,i)=handles.storage.scores(index,:);
end;

in = inpolygon(tmp(:,1),tmp(:,2),xi,yi);
indd=find(in);
handles.storage.selectedpoints = indd;
hold on; handles.selectedpointmarkers = plot(tmp(indd,1),tmp(indd,2),'go'); hold off;

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  save selected points                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_pcagui_saveselection_Callback(hObject, eventdata, handles)

handles = saveselection(handles);

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  show average of selected points                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_selectionaverage_Callback(hObject, eventdata, handles)

handles = generate_miniviewpoints(handles);  
handles = disp_average(handles);

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display selection in stackbrowser                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_pcagui_avg2stackbrowser_Callback(hObject, eventdata, handles)

handles = saveselection(handles,[pwd '/tmp.mat']);

if handles.storage.dims == 2 || handles.storage.dims == 1
    tom_av2_stackbrowser([pwd '/tmp.em'],[pwd '/tmp.mat']);
elseif handles.storage.dims == 3 
    tom_av3_stackbrowser([pwd '/tmp.mat']);
end


guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  reset selection                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_resetselection_Callback(hObject, eventdata, handles)

try;delete(handles.selectedpointmarkers);end;
handles.selectedpointmarkers = [];
guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  start clustering gui                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_pcagui_clustering_Callback(hObject, eventdata, handles)

inp.scores = handles.storage.scores;
inp.no_classes = 2;
index3 = str2num(get(handles.input_pcagui_eigenvector3,'String'))

if isempty(index3)
    inp.no_eigenvectors = [str2num(get(handles.input_pcagui_eigenvector1,'String')) str2num(get(handles.input_pcagui_eigenvector2,'String'))];
else
    inp.no_eigenvectors = [str2num(get(handles.input_pcagui_eigenvector1,'String')) str2num(get(handles.input_pcagui_eigenvector2,'String')) str2num(get(handles.input_pcagui_eigenvector3,'String'))];
end

handles.cluster = tom_clustergui(inp);



if isempty(handles.cluster)
    return;
end

set(handles.input_pcagui_eigenvector1,'String',num2str(handles.cluster.eigenvectors(1)));
set(handles.input_pcagui_eigenvector2,'String',num2str(handles.cluster.eigenvectors(2)));

if ~isempty(index3)
    set(handles.input_pcagui_eigenvector3,'String',num2str(handles.cluster.eigenvectors(3)));
end

%color code points
if isempty(index3)
    axes(handles.axes_pcagui_mainview);
else
    axes(findobj('Tag','3dcluster'));
end

hold on;
for i=1:size(handles.points,2)
    try;delete(handles.points(i));
    catch
        cla; 
        handles.points = []; handles.labels = [];
    end;
end

for i=1:size(handles.labels,2)
    try;delete(handles.labels(i));end;
end
handles.labels = [];
handles.points = [];
index1=str2num(get(handles.input_pcagui_eigenvector1,'String'));
index2=str2num(get(handles.input_pcagui_eigenvector2,'String'));
index3=str2num(get(handles.input_pcagui_eigenvector3,'String'));
colortable = tom_colorpalette(handles.cluster.noclasses);

for i=1:handles.cluster.noclasses
    idx = find(handles.cluster.classes==i);
    if isempty(index3)
        han = plot(handles.storage.scores(index1,idx),handles.storage.scores(index2,idx),'Marker','+','LineStyle','none','MarkerEdgeColor',colortable(i,:));
        handles.points = [handles.points,han];
        han = text(handles.cluster.centroid(i,1),handles.cluster.centroid(i,2),num2str(i),'BackgroundColor',[.7 .9 .7],'FontWeight','Bold');
        handles.labels = [handles.labels,han];
    else
        han = plot3(handles.storage.scores(index1,idx),handles.storage.scores(index2,idx),handles.storage.scores(index3,idx),'Marker','+','LineStyle','none','MarkerEdgeColor',colortable(i,:));
        handles.points = [handles.points,han];
        han = text(handles.cluster.centroid(i,1),handles.cluster.centroid(i,2),handles.cluster.centroid(i,3),num2str(i),'BackgroundColor',[.7 .9 .7],'FontWeight','Bold');
        handles.labels = [handles.labels,han];
    end
    
end

legstring = '';
for i=1:handles.cluster.noclasses
    legstring = [legstring, ',''', num2str(i) ''''];
end
eval(['legend(handles.points, ' legstring(2:end) ');']);
hold off;

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  select cluster class                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_pcagui_selectclass_Callback(hObject, eventdata, handles)

handles = select_clusterclass(handles,str2num(get(handles.input_pcagui_clusterclassno,'String')));

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  cluster gallery                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_pcagui_clustergallery_Callback(hObject, eventdata, handles)


if get(handles.checkbox_pcagui_miniview,'Value') == 1
    handles = generate_miniviewpoints(handles); 
    guidata(hObject, handles);
    for i=1:handles.cluster.noclasses
        handles = select_clusterclass(handles,i);
        guidata(hObject, handles);
        handles = disp_average(handles);
        guidata(hObject, handles);
    end
    
else
    
    if handles.storage.dims == 2
        avg = zeros(handles.storage.Header.Size(1),handles.storage.Header.Size(2),handles.cluster.noclasses);
        for i=1:handles.cluster.noclasses
            indd = find(handles.cluster.classes==i);
            for j=1:length(indd)
                part=tom_emreadc(handles.storage.stackfilename,'subregion',[1 1 indd(j)],[handles.storage.Header.Size(1)-1 handles.storage.Header.Size(1)-1 0]);
                part=part.Value;
                avg(:,:,i)=avg(:,:,i)+part;
            end
            avg(:,:,i) = avg(:,:,i)./length(indd);
        end

        tom_emwrite('tmp_pcaguigallery.em',avg);
        tom_av2_stackbrowser('tmp_pcaguigallery.em','makerefmode');

    else
        msgbox('Not implemented, use miniview');
    end
end

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  operate                                                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_operate_select_Callback(hObject, eventdata, handles)

mousebutton = 'normal';
i=1;
while strcmp(mousebutton,'normal') == 1
    waitforbuttonpress;
    axes2{i} = get(gca,'Tag');
    mousebutton = get(gcf,'SelectionType');
    i=i+1;
end

%get images

for i=1:size(axes2,2)
    userdata = get(handles.(axes2{i}), 'Userdata');
    try
        pointtype = userdata{1};
    catch
        return;
    end

    %eigen image
    if isequal(pointtype,'eigen')
        [handles, im_out{i}] = get_eigenimage(handles,userdata{2},1);
    %orig image
    elseif isequal(pointtype,'orig')
        [handles, im_out{i}] = get_origimage(handles,userdata{2},1);
    %average
    else
        if handles.storage.dims==2
            im_out{i}=zeros(handles.storage.Header.Size(1),handles.storage.Header.Size(2));
        end

        indd = userdata{3};

        for j=1:length(indd)
            if (handles.storage.dims == 3)
                al(j) = handles.storage.align(indd(j));
            else
                part=tom_emreadc(handles.storage.stackfilename,'subregion',[1 1 indd(j)],[handles.storage.Header.Size(1)-1 handles.storage.Header.Size(1)-1 0]);
                part=part.Value;
                im_out{i}=im_out{i}+part;
            end

        end

        if handles.storage.dims == 3
            im_out{i} = tom_av3_average(al,'sum',0,0,1);
        end
    end
    
end

norm_string='''mean0+1std''';
%operate on images
operator = get(handles.tom_pcagui_operator,'String');
evstring = '';
for i=1:size(axes2,2)
    evstring = [evstring 'tom_norm(im_out{' num2str(i) '},'  norm_string ')'];
    %evstring = ['' evstring ',' norm_string ')'];
    if i<size(axes2,2)
        evstring = [evstring operator];
    end
end
out = eval(evstring);    

%display result
if handles.storage.dims == 2
    figure;tom_imagesc(out);
else
    if get(handles.radiobutton_pcagui_dspcub,'Value') == 1
        figure;tom_dspcub(out);
    else
        tom_volxyz(out);
    end
end
  

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  show 1d particles                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_pcagui_show1dparticles_Callback(hObject, eventdata, handles)

%read particles from stack
if get(handles.checkbox_pcagui_1dstack,'Value') == 1
    
    if isempty(handles.storage.onedstackfile)
        [filename, pathname] = uigetfile({'*.em'}, 'Load stack file');
        if ~isempty(filename)
            handles.storage.onedstackfile = [pathname '/' filename];
        end
    end
    handles = disp_average(handles, handles.storage.onedstackfile);
    
    
%extract particles directly from original image files
else
    if isempty(handles.storage.onedalignfile)
        [filename, pathname] = uigetfile({'*.mat'}, 'Load alignment file');
        if ~isempty(filename)
            handles.storage.onedalignfile = [pathname '/' filename];
        end
    end
    s = load(handles.storage.onedalignfile);
    
    %filter selected particles
    indd = handles.storage.selectedpoints;
    for i=1:length(indd)
        s.align2d(1,indd(i)).selectedtowrite = 1;
    end
    al_out = [];
    lauf = 1;
    for k=1:size(s.align2d,2)
        if s.align2d(1,k).selectedtowrite == 1
           try;al_out(lauf) = s.align2d(1,k);catch;al_out = s.align2d(1,k);end;
           lauf=lauf+1; 
        end
    end
    stack = tom_av2_createstack(al_out, '', '', '', 1, 0, s.align2d(1,1).radius);
    handles = generate_miniviewpoints(handles);
    handles = disp_generic(handles, sum(stack,3));
end

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
%%  load file                                                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = loadfile(handles)

filename = get(handles.input_pcagui_inputfile,'String');

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
set(handles.popupmenu_pcagui_run,'String',liststring);
set(handles.popupmenu_pcagui_run,'Enable','on');
set(handles.input_pcagui_binning,'Enable','on');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  save selection to file                                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = saveselection(handles, name)

if handles.storage.dims == 1
   points = handles.storage.selectedpoints;
   save('points.mat','points');
   return;
end

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
indd = handles.storage.selectedpoints;
for i=1:length(indd)
    handles.storage.align(1,indd(i)).selectedtowrite = 1;
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
        tom_av2_filterstack(handles.storage.stackfilename,handles.storage.align,'selectedtowrite','==',1,outstackfile,outalignfile);
    elseif handles.storage.dims ==1
        tom_av2_filterstack(stackfile2,handles.storage.align,'selectedtowrite','==',1,outstackfile,outalignfile);
    else
        al = handles.storage.align;
        lauf=1;
        j = size(al,1);
        for k=1:size(al,2)
            if al(1,k).selectedtowrite == 1
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
handles.storage.align = rmfield(handles.storage.align,'selectedtowrite');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  calculate eigenimage                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [handles, im_out] = get_eigenimage(handles,coords,full3dflag);

if nargin == 2
    full3dflag = 0;
end

X = coords(1);
Y = coords(2);

v_tmp=[X Y];

eigvecs(1) = str2num(get(handles.input_pcagui_eigenvector1,'String'));
eigvecs(2) = str2num(get(handles.input_pcagui_eigenvector2,'String'));

 for i=1:size(eigvecs,2)
     index=eigvecs(i);
     tmp(:,i)=handles.storage.coefs(:,index);
 end;

img_tmp=tmp*v_tmp';

input_binning = str2num(get(handles.input_pcagui_binning,'String'));

if handles.storage.dims==3
    sz=handles.storage.Header.Size;
    sz=round(sz./(2^input_binning));   
    im_out = reshape(tom_rm_mean(img_tmp,handles.storage.mean),sz(1),sz(2),sz(3));
    if full3dflag == 0
        im_out = sum(im_out(:,:,sz(3)./2-3:sz(3)./2+3),3);
    end
elseif handles.storage.dims==2
    sz=round(handles.storage.Header.Size');
    sz=round(sz./(2^input_binning));
    im_out = reshape(tom_rm_mean(img_tmp,handles.storage.mean),sz(1),sz(2));
else
    im_out = tom_rm_mean(img_tmp,handles.storage.mean);
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  calculate original image                                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [handles, im_out, coords_out] = get_origimage(handles,coords,full3dflag);

if nargin == 2
    full3dflag = 0;
end

X = coords(1);
Y = coords(2);

v_tmp=[X Y];

eigvecs(1) = str2num(get(handles.input_pcagui_eigenvector1,'String'));
eigvecs(2) = str2num(get(handles.input_pcagui_eigenvector2,'String'));

for i=1:size(eigvecs,2)
    index=eigvecs(i);
    tmp(:,i)=handles.storage.scores(index,:);
end;

[pointidx, coords_out] = tom_nearestpoint(v_tmp,tmp);


if handles.storage.dims==3
    part=tom_emreadc(handles.storage.align(pointidx).Filename);
    im_out = part.Value;
    if full3dflag == 0
        im_out = sum(im_out(:,:,part.Header.Size(3)./2-4:part.Header.Size(3)./2+4),3);
    end
elseif handles.storage.dims ==2

    im_out=tom_emreadc(handles.storage.stackfilename,'subregion',[1 1 pointidx],[handles.storage.Header.Size(1)-1 handles.storage.Header.Size(2)-1 0]);
    im_out = im_out.Value;
else
    im_out=tom_emreadc(handles.storage.stackfilename,'subregion',[pointidx 1 1],[0 handles.storage.Header.Size(2)-1 0]);
    im_out = im_out.Value;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  get eigenvector gallery                                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [handles, im_out] = get_eigenvectorgallery(handles)

eigvecs_gal(1) = str2num(get(handles.input_eigengallery_start,'String'));
eigvecs_gal(2) = str2num(get(handles.input_eigengallery_stop,'String'));

if handles.storage.dims==3
   
    input_binning = str2num(get(handles.input_pcagui_binning,'String'));
    
    sz=handles.storage.Header.Size;
    sz=round(sz./(2^input_binning)); 
    
    start_y=1;
    for i=eigvecs_gal(1):eigvecs_gal(2)
        index=i;
        vol_tmp=reshape(handles.storage.coefs(:,index),sz(1),sz(2),sz(3));
        start_x=1;
        for ii=1:size(vol_tmp,3)
           im_out(start_x:start_x+sz(1)-1,start_y:(start_y+sz(2)-1))=vol_tmp(:,:,ii);
           start_x=start_x+sz(1);
        end;
        
        start_y=start_y+sz(2);    
    end;

elseif handles.storage.dims == 2
    sz=sqrt(size(handles.storage.coefs,1));
    for i=eigvecs_gal(1):eigvecs_gal(2)
        index=i;
        im_out(:,:,i)=reshape(handles.storage.coefs(:,index),sz,sz);
    end;

else
    figure;l=1;
    for i=eigvecs_gal(1):eigvecs_gal(2)
        subplot(ceil(eigvecs_gal(2)-eigvecs_gal(1)./5),5,l);plot(handles.storage.coefs(:,i));l=l+1;axis tight;axis off;
    end
    im_out = [];
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  calculate centroid of point selection                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [coords, handles] = calc_centroid(handles)

indd = handles.storage.selectedpoints;

eigvecs(1) = str2num(get(handles.input_pcagui_eigenvector1,'String'));
eigvecs(2) = str2num(get(handles.input_pcagui_eigenvector2,'String'));

% get coordinates
for i=1:size(eigvecs,2)
    index=eigvecs(i);
    tmp(:,i)=handles.storage.scores(index,indd);
end;

coords(1) = sum(tmp(:,1))./length(indd);
coords(2) = sum(tmp(:,2))./length(indd);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  select cluster class                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = select_clusterclass(handles,no)
index1=str2num(get(handles.input_pcagui_eigenvector1,'String'));
index2=str2num(get(handles.input_pcagui_eigenvector2,'String'));
index3=str2num(get(handles.input_pcagui_eigenvector3,'String'));

try;delete(handles.selectedpointmarkers);end;
handles.selectedpointmarkers = [];
handles.storage.selectedpoints = [];
for i=no
    handles.storage.selectedpoints = [handles.storage.selectedpoints, find(handles.cluster.classes==i)];
end
hold on;
if isempty(index3)
    handles.selectedpointmarkers = plot(handles.storage.scores(index1,handles.storage.selectedpoints),handles.storage.scores(index2,handles.storage.selectedpoints),'go');
else
    axes(findobj('Tag','3dcluster'));
    hold on;
    handles.selectedpointmarkers = plot3(handles.storage.scores(index1,handles.storage.selectedpoints),handles.storage.scores(index2,handles.storage.selectedpoints),handles.storage.scores(index3,handles.storage.selectedpoints),'go');
end
hold off;
handles.dispclassno = ['class ' num2str(no)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  generate miniviewpoints                                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = generate_miniviewpoints(handles)

handles.miniviewpoints = [];  
handles.xlim = get(handles.axes_pcagui_mainview,'Xlim');
handles.ylim = get(handles.axes_pcagui_mainview,'Ylim');

xabs = (handles.xlim(2)-handles.xlim(1))./8;
yabs = (handles.ylim(2)-handles.ylim(1))./8;

for i=1:7
    handles.xpoints(i) = handles.xlim(1)+(i).*xabs;
    handles.ypoints(8-i) = handles.ylim(1)+(i).*yabs;
end

for i=1:7
    handles.miniviewpoints = cat(1,handles.miniviewpoints,[handles.xlim(1),handles.ypoints(i)]);
    handles.miniviewpoints = cat(1,handles.miniviewpoints,[handles.xlim(2),handles.ypoints(i)]);
    handles.miniviewpoints = cat(1,handles.miniviewpoints,[handles.xpoints(i),handles.ylim(1)]);
    handles.miniviewpoints = cat(1,handles.miniviewpoints,[handles.xpoints(i),handles.ylim(2)]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display average selection                                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = disp_average(handles, filename)

if nargin == 1
    if handles.storage.dims ==1 
        filename = handles.storage.stackfilename;
    elseif handles.storage.dims == 2
        filename = handles.storage.stackfilename;
    elseif handles.storage.dims == 3
        filename = handles.storage.alignallruns(1).Filename;
    end
end

header = tom_reademheader(filename);

if handles.storage.dims == 2 || nargin == 2
    avg=zeros(header.Header.Size(1),header.Header.Size(2));
elseif handles.storage.dims == 1 && nargin ~= 2
    avg = zeros(1,header.Header.Size(2));
end

indd = handles.storage.selectedpoints;

idx=0;

for i=1:length(indd)
    if handles.storage.dims == 3
        al(i) = handles.storage.align(indd(i)); 
    elseif handles.storage.dims == 2 || nargin == 2
        part=tom_emreadc(filename,'subregion',[1 1 indd(i)],[header.Header.Size(1)-1 header.Header.Size(2)-1 0]);
        part=part.Value;
        avg=avg+part;
        idx=idx+1;
    else
        part=tom_emreadc(filename,'subregion',[indd(i) 1 1],[0 handles.storage.Header.Size(2)-1 0]);
        part=part.Value;
        avg=avg+part;
        idx=idx+1;
    end
    
end

if idx>0
    avg=avg./idx;
end;

if handles.storage.dims == 3
    avg = tom_av3_average(al,'sum',0,0,1);
end

if get(handles.checkbox_pcagui_miniview,'Value') == 1
    
    [startcoords, handles] = calc_centroid(handles);
    if isempty(find(handles.miniviewpoints>-1e99))
        return;
    end 
    [idx,endcoords] = tom_nearestpoint(startcoords,handles.miniviewpoints);
    
    if endcoords(1) == handles.xlim(1)
        %left
        orient = 'left';
        no = find(handles.ypoints==endcoords(2));
    elseif endcoords(1) == handles.xlim(2)
        %right
        orient = 'right';
        no = find(handles.ypoints==endcoords(2));
    elseif endcoords(2) == handles.ylim(1)
        %bottom
        orient = 'bottom';
        no = find(handles.xpoints==endcoords(1));
    elseif endcoords(2) == handles.ylim(2)
        %top
        orient = 'top';
        no = find(handles.xpoints==endcoords(1));
    end
    
    axes2 = ['axes_pcagui_' orient num2str(no)];
    handles.miniviewpoints(idx,:) = [-9e100 -9e100];
    
    
    %draw image
    if handles.storage.dims == 3
        part=tom_reademheader(handles.storage.align(1).Filename);
        avg = sum(avg(:,:,part.Header.Size(3)./2-4:part.Header.Size(3)./2+4),3);
    end

    axes(handles.(axes2));
    if handles.storage.dims > 1 || nargin == 2
        imagesc(avg');colormap gray;%axis ij;
    else
        plot(avg);
    end
    set(handles.(axes2),'XTick',[],'YTick',[],'Tag',axes2,'Userdata',{'avg' [startcoords(1) startcoords(2)] [handles.storage.selectedpoints]});
    [b c d] = strread(axes2,'%s %s %s','delimiter','_');
    if ~isempty(handles.dispclassno)
        set(handles.(['text_pcagui_' d{1}]),'String',strvcat('avg',handles.dispclassno));
        handles.dispclassno = '';
    else
        set(handles.(['text_pcagui_' d{1}]),'String','avg');
    end
        
    %render line
    axes(handles.axes_pcagui_mainview);
    hold on;
    if isfield(handles.line, axes2)
        try;delete(handles.line.(axes2));end;
        try;delete(handles.point.(axes2));end;
    end
    handles.point.(axes2) = plot(startcoords(1),startcoords(2),'Color',[0 0 0],'Marker','o');
    handles.line.(axes2) = line([startcoords(1) endcoords(1)],[startcoords(2) endcoords(2)],'Color',[0 0 0],'LineWidth',1);
    hold off;
    set(handles.axes_pcagui_mainview,'Tag','axes_pcagui_mainview');

else

    if (handles.storage.dims==3)
        tom_volxyz(avg);
    end;
    if (handles.storage.dims==2)    
        figure; tom_imagesc(avg);
    end;
    if (handles.storage.dims==1)    
        figure; plot(avg);
    end;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display image generic                                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = disp_generic(handles, avg)

if get(handles.checkbox_pcagui_miniview,'Value') == 1
    
    [startcoords, handles] = calc_centroid(handles);
    if isempty(find(handles.miniviewpoints>-1e99))
        return;
    end 
    [idx,endcoords] = tom_nearestpoint(startcoords,handles.miniviewpoints);
    
    if endcoords(1) == handles.xlim(1)
        %left
        orient = 'left';
        no = find(handles.ypoints==endcoords(2));
    elseif endcoords(1) == handles.xlim(2)
        %right
        orient = 'right';
        no = find(handles.ypoints==endcoords(2));
    elseif endcoords(2) == handles.ylim(1)
        %bottom
        orient = 'bottom';
        no = find(handles.xpoints==endcoords(1));
    elseif endcoords(2) == handles.ylim(2)
        %top
        orient = 'top';
        no = find(handles.xpoints==endcoords(1));
    end
    
    axes2 = ['axes_pcagui_' orient num2str(no)];
    handles.miniviewpoints(idx,:) = [-9e100 -9e100];
    
    
    %draw image
    if handles.storage.dims == 3
        part=tom_reademheader(handles.storage.align(1).Filename);
        avg = sum(avg(:,:,part.Header.Size(3)./2-4:part.Header.Size(3)./2+4),3);
    end

    axes(handles.(axes2));
    if handles.storage.dims > 1 || nargin == 2
        imagesc(avg');colormap gray;%axis ij;
    else
        plot(avg);
    end
    set(handles.(axes2),'XTick',[],'YTick',[],'Tag',axes2,'Userdata',{'avg' [startcoords(1) startcoords(2)] [handles.storage.selectedpoints]});
    [b c d] = strread(axes2,'%s %s %s','delimiter','_');
    set(handles.(['text_pcagui_' d{1}]),'String','avg');

    %render line
    axes(handles.axes_pcagui_mainview);
    hold on;
    if isfield(handles.line, axes2)
        try;delete(handles.line.(axes2));end;
        try;delete(handles.point.(axes2));end;
    end
    handles.point.(axes2) = plot(startcoords(1),startcoords(2),'Color',[0 0 0],'Marker','o');
    handles.line.(axes2) = line([startcoords(1) endcoords(1)],[startcoords(2) endcoords(2)],'Color',[0 0 0],'LineWidth',1);
    hold off;
    set(handles.axes_pcagui_mainview,'Tag','axes_pcagui_mainview');

else

    if (handles.storage.dims==3)
        tom_volxyz(avg);
    else
        figure; tom_imagesc(avg);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Unused Callbacks/Create functions                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function slider1_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function slider2_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function input_pcagui_inputfile_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function popupmenu_pcagui_binning_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function popupmenu_pcagui_run_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function popupmenu_pcagui_run_Callback(hObject, eventdata, handles)
function input_pcagui_numeigenvectors_Callback(hObject, eventdata, handles)
function input_pcagui_numeigenvectors_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_pcagui_binning_Callback(hObject, eventdata, handles)
function input_pcagui_binning_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_pcagui_eigenvector1_Callback(hObject, eventdata, handles)
function input_pcagui_eigenvector1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_pcagui_eigenvector2_Callback(hObject, eventdata, handles)
function input_pcagui_eigenvector2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function dummy(a,b)
function input_eigengallery_stop_Callback(hObject, eventdata, handles)
function input_eigengallery_stop_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_eigengallery_start_Callback(hObject, eventdata, handles)
function input_eigengallery_start_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function togglebutton2_Callback(hObject, eventdata, handles)
function checkbox_pcagui_phasenorm_Callback(hObject, eventdata, handles)
function input_pcagui_clusterclassno_Callback(hObject, eventdata, handles)
function input_pcagui_clusterclassno_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function checkbox_pcagui_miniview_Callback(hObject, eventdata, handles)
function input_pcagui_eigenvector3_Callback(hObject, eventdata, handles)
function input_pcagui_eigenvector3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_pcagui_loadcalc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function tom_pcagui_operator_Callback(hObject, eventdata, handles)
function tom_pcagui_operator_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function popupmenu_pcagui_space_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function checkbox_pcagui_cpca_Callback(hObject, eventdata, handles)
function popupmenu_pcagui_space_Callback(hObject, eventdata, handles)
function checkbox_pcagui_1dstack_Callback(hObject, eventdata, handles)
function input_pcagui_loadcalc_Callback(hObject, eventdata, handles)
function pushbutton_pcagui_browsecalc_Callback(hObject, eventdata, handles)
