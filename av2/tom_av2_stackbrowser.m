function varargout = tom_av2_stackbrowser(varargin)
%TOM_AV2_STACKBROWSER is a GUI for browsing 2D particle stacks
%
%   varargout = tom_av2_stackbrowser(varargin)
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
%   ... = tom_av2_stackbrowser(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by .. 01/12/06
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
                   'gui_OpeningFcn', @tom_av2_stackbrowser_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_av2_stackbrowser_OutputFcn, ...
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
function tom_av2_stackbrowser_OpeningFcn(hObject, eventdata, handles, varargin)

global storage_av2_stackbrowser;

storage_av2_stackbrowser.stack = [];
storage_av2_stackbrowser.align = [];
storage_av2_stackbrowser.numcols = 10;
storage_av2_stackbrowser.numimages = 100;
storage_av2_stackbrowser.pos = 1;
storage_av2_stackbrowser.imsize = [0 0 0];
storage_av2_stackbrowser.page = 0;
storage_av2_stackbrowser.display.marks = 1;
storage_av2_stackbrowser.display.ccc = 1;
storage_av2_stackbrowser.display.quality = 1;
storage_av2_stackbrowser.display.ref_classno = 1;
storage_av2_stackbrowser.quality = 1;
storage_av2_stackbrowser.fontsize = 20;
storage_av2_stackbrowser.bandpass.enable = 0;
storage_av2_stackbrowser.bandpass.low = 0;
storage_av2_stackbrowser.bandpass.high = 0;
storage_av2_stackbrowser.makerefmode = 0;
storage_av2_stackbrowser.pcamode = 0;
storage_av2_stackbrowser.sn_show_parts=0;
storage_av2_stackbrowser.stack_in_mem=[];
storage_av2_stackbrowser.incre_load=1;
storage_av2_stackbrowser.sh_part_filt_kern=2;
storage_av2_stackbrowser.sh_picks=0;
storage_av2_stackbrowser.stack_is_compressed=0;
storage_av2_stackbrowser.stackfile_org=[];
storage_av2_stackbrowser.alignmentfile_org=[];
storage_av2_stackbrowser.align_org=[];

set(findobj('Tag','input_quality'),'String',num2str(storage_av2_stackbrowser.quality));

axes(findobj('Tag','imchoosebox'));axis off;

%load stack if given as input parameter
if nargin > 3
    storage_av2_stackbrowser.stackfile = varargin{1};
    set(findobj('Tag','input_particlestack'),'String',varargin{1});
    
    if size(varargin,2) == 1 || (size(varargin,2) == 1 && strcmp(varargin{2},'makerefmode') == 1)
        storage_av2_stackbrowser.align = tom_av2_create_alignfromstack(storage_av2_stackbrowser.stackfile);
        storage_av2_stackbrowser.alignmentfile = '';
    else
        if strcmp(varargin{2},'makerefmode') == 1
            storage_av2_stackbrowser.alignmentfile = '';
            storage_av2_stackbrowser.makerefmode = 1;
            set(findobj('Tag','panel_filter'),'Visible','off');
            set(findobj('Tag','panel_save'),'Visible','off');
            set(findobj('Tag','display_ccc'),'Value',0);
            set(findobj('Tag','display_quality'),'Value',0);
            set(findobj('Tag','button_stackbrowse_makeref'),'Visible','on');
            storage_av2_stackbrowser.display.ccc = 0;
            storage_av2_stackbrowser.display.quality = 0;
            storage_av2_stackbrowser.align = tom_av2_create_alignfromstack(storage_av2_stackbrowser.stackfile);
        elseif strcmp(varargin{2},'pcamode') == 1
            storage_av2_stackbrowser.alignmentfile = varargin{3};
            s = load(storage_av2_stackbrowser.alignmentfile);
            storage_av2_stackbrowser.align = s.align2d;
            set(findobj('Tag','input_alignmentfile'),'String',varargin{3});
            storage_av2_stackbrowser.pcamode = 1;
            set(findobj('Tag','button_av2stackbrowser_pca'),'Visible','on');
        else
            storage_av2_stackbrowser.alignmentfile = varargin{2};
            s = load(storage_av2_stackbrowser.alignmentfile);
            storage_av2_stackbrowser.align = s.align2d;
            set(findobj('Tag','input_alignmentfile'),'String',varargin{2});
        end
    end
    
    for i=1:size(storage_av2_stackbrowser.align,2)
         storage_av2_stackbrowser.align(i).partno = i;
    end

    if ~isfield(storage_av2_stackbrowser.align,'ref_class')
        for i=1:size(storage_av2_stackbrowser.align,2)
            storage_av2_stackbrowser.align(i).ref_class = 0;
        end
    end
    

    if storage_av2_stackbrowser.makerefmode == 0 && storage_av2_stackbrowser.pcamode == 0
        %select all particles by default
        for i=1:size(storage_av2_stackbrowser.align,2)
            storage_av2_stackbrowser.align(i).selected = 1;
        end
    else
        for i=1:size(storage_av2_stackbrowser.align,2)
            storage_av2_stackbrowser.align(i).selected = 0;
        end
    end
    
    loadstack(1);
    numpages = floor(size(storage_av2_stackbrowser.align,2)./storage_av2_stackbrowser.numimages);
    if numpages > 1
        %set(findobj('Tag','slider_particles'),'Max',numpages,'SliderStep',[1./numpages 1./numpages]);
      %  set(findobj('Tag','slider_particles'),'Max',numpages,'SliderStep',[1./numpages 1./numpages]);
    else
        set(findobj('Tag','slider_particles'),'Enable','off');
    end
    set(findobj('Tag','text_pageinfo'),'String',['page 1 of ' num2str(ceil(size(storage_av2_stackbrowser.align,2)./storage_av2_stackbrowser.numimages))]);
    set(findobj('Tag','button_loadstack'),'Enable','off');

    render_image();
end

handles.output = hObject;

guidata(hObject, handles);

if storage_av2_stackbrowser.makerefmode == 1
    uiwait(handles.figure1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Output Function                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = tom_av2_stackbrowser_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

global storage_av2_stackbrowser;

if storage_av2_stackbrowser.makerefmode == 1
    close(handles.figure1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  directy input particle stack                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_particlestack_Callback(hObject, eventdata, handles)

global storage_av2_stackbrowser;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  browse particle stack                                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_particlestackbrowse_Callback(hObject, eventdata, handles)

global storage_av2_stackbrowser;

[filename, pathname] = uigetfile({'*.em'}, 'Pick an image stack');
if ischar(filename)
    if tom_isemfile([pathname '/' filename])
        set(findobj('Tag','input_particlestack'),'String',[pathname '/' filename]);
        storage_av2_stackbrowser.stackfile = [pathname '/' filename];
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  directy input alignment file                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_alignmentfile_Callback(hObject, eventdata, handles)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  browse alignment file                                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_alignmentfilebrowse_Callback(hObject, eventdata, handles)

global storage_av2_stackbrowser;

[filename, pathname] = uigetfile({'*.mat'}, 'Pick an alignment file');
if ischar(filename)
    if exist([pathname '/' filename],'file')
        set(findobj('Tag','input_alignmentfile'),'String',[pathname '/' filename]);
        storage_av2_stackbrowser.alignmentfile = [pathname '/' filename];
    end
end

s = load(storage_av2_stackbrowser.alignmentfile);
storage_av2_stackbrowser.align = s.align2d;

%select all particles by default
for i=1:size(storage_av2_stackbrowser.align,2)
    storage_av2_stackbrowser.align(i).selected = 1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  delete selected particles                                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_particlesdelete_Callback(hObject, eventdata, handles)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  invert selection                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_invertselection_Callback(hObject, eventdata, handles)

global storage_av2_stackbrowser;

if (get(handles.chk_inv_all,'Value')==0)
    
    %draw marks
    delete(findobj('-regexp','Tag','rect_*'));
    if storage_av2_stackbrowser.display.marks == 1
        numrows = storage_av2_stackbrowser.numimages ./ storage_av2_stackbrowser.numcols;
        for row=1:numrows
            pos_y = (row-1).*storage_av2_stackbrowser.imsize(2)+storage_av2_stackbrowser.imsize(2).*0.1;
            for column=1:storage_av2_stackbrowser.numcols
                abspartnumber = storage_av2_stackbrowser.page.*storage_av2_stackbrowser.numimages+(row-1).*storage_av2_stackbrowser.numcols+column;
                if abspartnumber <= size(storage_av2_stackbrowser.align,2)
                    %mark particle green if selected and red if not selected
                    if storage_av2_stackbrowser.align(abspartnumber).selected == 1
                        storage_av2_stackbrowser.align(abspartnumber).selected = 0;
                        color = [1 0 0];
                    else
                        storage_av2_stackbrowser.align(abspartnumber).selected = 1;
                        color = [0 1 0];
                    end
                    pos_x = (column-1).*storage_av2_stackbrowser.imsize(1)+storage_av2_stackbrowser.imsize(1).*0.1;
                    r = rectangle('Position',[pos_x,pos_y,storage_av2_stackbrowser.imsize(1).*0.1,storage_av2_stackbrowser.imsize(2).*0.1],'FaceColor',color,'Tag',['rect_' abspartnumber]);
                end
            end
        end
    end
    
else
    for i=1:length(storage_av2_stackbrowser.align)
        if storage_av2_stackbrowser.align(i).selected == 1
            storage_av2_stackbrowser.align(i).selected = 0;
        else
            storage_av2_stackbrowser.align(i).selected = 1;
        end
    end;
    render_image();
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  slider callback                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function slider_particles_Callback(hObject, eventdata, handles)

global storage_av2_stackbrowser;

set(findobj('Tag','text_pageinfo'),'String',['page ' num2str(round(get(hObject,'Value')+1)) ' of ' num2str(ceil(size(storage_av2_stackbrowser.align,2)./storage_av2_stackbrowser.numimages))]);

position = round(get(hObject,'Value')).*storage_av2_stackbrowser.numimages+1;
storage_av2_stackbrowser.page = round(get(hObject,'Value'));
loadstack(position);
drawnow;
render_image();
drawnow;
set(findobj('Tag','input_nocols'),'Enable','off');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  load stack                                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_loadstack_Callback(hObject, eventdata, handles)

global storage_av2_stackbrowser;

if isempty(storage_av2_stackbrowser.stackfile)
    errordlg('Load stack file first.');
    return;   
end

if isempty(storage_av2_stackbrowser.align)
    storage_av2_stackbrowser.align = tom_av2_create_alignfromstack(storage_av2_stackbrowser.stackfile);
    %select all particles by default
    for i=1:size(storage_av2_stackbrowser.align,2)
        storage_av2_stackbrowser.align(i).selected = 1;
    end
end

for i=1:size(storage_av2_stackbrowser.align,2)
         storage_av2_stackbrowser.align(i).partno = i;
end

if ~isfield(storage_av2_stackbrowser.align,'ref_class')
    for i=1:size(storage_av2_stackbrowser.align,2)
        storage_av2_stackbrowser.align(i).ref_class = 0;
    end
end

loadstack(1);

% if (storage_av2_stackbrowser.incre_load==1)
%     
% else
%     
% end;


numpages = ceil(size(storage_av2_stackbrowser.align,2)./storage_av2_stackbrowser.numimages);

if numpages > 1
    %set(findobj('Tag','slider_particles'),'Min',0,'Max',numpages,'SliderStep',[1./numpages 1./numpages]);
    set(findobj('Tag','slider_particles'),'Min',0,'Max',numpages-1,'SliderStep',[1./(numpages-1) 1./(numpages-1)]);
else
    set(findobj('Tag','slider_particles'),'Enable','off');
end
set(findobj('Tag','text_pageinfo'),'String',['page 1 of ' num2str(ceil(size(storage_av2_stackbrowser.align,2)./storage_av2_stackbrowser.numimages))]);
set(hObject,'Enable','off');


render_image();


% --- Executes on button press in button_compress_stack.
function button_compress_stack_Callback(hObject, eventdata, handles)
% hObject    handle to button_compress_stack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global storage_av2_stackbrowser;

in.input_str= storage_av2_stackbrowser.stackfile;

st=tom_av2_compress_stack_gui(in);

tmp_stack=tom_emreadc(st.input_str);


[st_comp align2d]=tom_av2_compress_stack(tmp_stack.Value,st.e_num_of_classes,'not_equal',st.e_pre_alg_th,eval(st.e_eigs),st.e_bin,st.e_cluster_method,'',0);

[a_st b_st c_st]=fileparts(storage_av2_stackbrowser.stackfile);
[a_al b_al c_al]=fileparts(storage_av2_stackbrowser.alignmentfile);


tom_emwrite([a_st '/' b_st '_comp' c_st],st_comp);

save([a_al '/' b_al '_comp' c_al],'align2d');


storage_av2_stackbrowser.stackfile_org=storage_av2_stackbrowser.stackfile;
storage_av2_stackbrowser.alignmentfile_org=storage_av2_stackbrowser.alignmentfile;
storage_av2_stackbrowser.align_org=storage_av2_stackbrowser.align;

storage_av2_stackbrowser.stackfile=[a_st '/' b_st '_comp' c_st];
storage_av2_stackbrowser.alignmentfile=[a_al '/' b_al '_comp' c_al];

storage_av2_stackbrowser.align=align2d;


for i=1:size(storage_av2_stackbrowser.align,2)
    storage_av2_stackbrowser.align(i).selected = 1;
    storage_av2_stackbrowser.align(i).partno = i;
end

loadstack(1);
numpages = ceil(size(storage_av2_stackbrowser.align,2)./storage_av2_stackbrowser.numimages);

if numpages > 1
    set(findobj('Tag','slider_particles'),'Min',0,'Max',numpages-1,'SliderStep',[1./(numpages-1) 1./(numpages-1)]);
    %set(findobj('Tag','slider_particles'),'Min',0,'Max',numpages,'SliderStep',[1./numpages 1./numpages]);
else
    set(findobj('Tag','slider_particles'),'Enable','off');
end
set(findobj('Tag','text_pageinfo'),'String',['page 1 of ' num2str(ceil(size(storage_av2_stackbrowser.align,2)./storage_av2_stackbrowser.numimages))]);
set(hObject,'Enable','off');

set(handles.chk_show_parts,'visible','on');
set(handles.sh_part_filt_kern,'visible','on');
set(handles.chk_show_picks,'visible','on');

storage_av2_stackbrowser.stack_in_mem=tmp_stack;
storage_av2_stackbrowser.stack_is_compressed=1;

render_image();






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  number of columns                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_nocols_Callback(hObject, eventdata, handles)

global storage_av2_stackbrowser;

val = str2num(get(hObject,'String'));

if ~isempty(val) & val > 1
    storage_av2_stackbrowser.numcols = round(val);
    storage_av2_stackbrowser.numimages = storage_av2_stackbrowser.numcols.^2;
end
loadstack(storage_av2_stackbrowser.pos);
render_image();
numpages = floor(size(storage_av2_stackbrowser.align,2)./storage_av2_stackbrowser.numimages);
set(findobj('Tag','slider_particles'),'Max',numpages,'SliderStep',[1./numpages 1./numpages]);
set(findobj('Tag','text_pageinfo'),'String',['page 1 of ' num2str(ceil(size(storage_av2_stackbrowser.align,2)./storage_av2_stackbrowser.numimages))]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display good/bad marks                                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function display_marks_Callback(hObject, eventdata, handles)

global storage_av2_stackbrowser;

storage_av2_stackbrowser.display.marks = get(hObject,'Value');
render_image();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display CCC values                                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function display_ccc_Callback(hObject, eventdata, handles)

global storage_av2_stackbrowser;

storage_av2_stackbrowser.display.ccc = get(hObject,'Value');
render_image();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display quality                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function display_quality_Callback(hObject, eventdata, handles)

global storage_av2_stackbrowser;

storage_av2_stackbrowser.display.quality = get(hObject,'Value');
render_image();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display classno                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function display_classno_Callback(hObject, eventdata, handles)

global storage_av2_stackbrowser;

storage_av2_stackbrowser.display.ref_classno = get(hObject,'Value');
render_image();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  save to stack                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_savetostack_Callback(hObject, eventdata, handles)

global storage_av2_stackbrowser;

rest_quality_flag=get(handles.checkbox_reset_quality,'Value');

for i=1:size(storage_av2_stackbrowser.align,2)
    if storage_av2_stackbrowser.align(i).selected == 1
        storage_av2_stackbrowser.align(i).quality = storage_av2_stackbrowser.quality;
    else
        if (rest_quality_flag==1)
            storage_av2_stackbrowser.align(i).quality = 0;
        end;
    end;
end;


align2d = storage_av2_stackbrowser.align;
align2d = rmfield(align2d,'selected');
align2d = rmfield(align2d,'partno');

%outfile = get(findobj('Tag','input_alignmentfile'),'String');

outfile = storage_av2_stackbrowser.alignmentfile;

if isempty(outfile)
    [filename, pathname] = uigetfile({'*.mat'}, 'Pick an output alignment file');
    if ischar(filename)
        outfile = [pathname '/' filename];
        save(outfile,'align2d');
    end
else
    save(outfile,'align2d');
end
render_image();

disp('Alignment File saved.');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  input quality                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_quality_Callback(hObject, eventdata, handles)

global storage_av2_stackbrowser;

quality = round(str2num(get(hObject,'String')));

if quality > 10 | quality < 0
    errordlg('quality must be between 0 and 10.');
    set(hObject,'String',num2str(storage_av2_stackbrowser.quality));
    return;
end

storage_av2_stackbrowser.quality = quality;
set(hObject,'String',num2str(storage_av2_stackbrowser.quality));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  histogram reset                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_parthisto_reset_Callback(hObject, eventdata, handles)

calc_histogram();
render_image();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  histogram set                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_parthisto_set_Callback(hObject, eventdata, handles)

global storage_av2_stackbrowser;

axes(findobj('Tag','partbrowser_histogram'));
waitforbuttonpress;
point1 = get(gca,'CurrentPoint');    % button down detected
finalRect = rbbox;                   % return figure units
point2 = get(gca,'CurrentPoint');    % button up detected
point1 = point1(1,1:2);              % extract x and y
point2 = point2(1,1:2);
p1 = min(point1,point2);             % calculate locations
offset = abs(point1-point2);         % and dimensions
x = [p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)];
set(gca,'Xlim',[x(1) x(2)]);

storage_av2_stackbrowser.Scale=[x(1) x(2)];
render_image();
set(findobj('Tag','input_parthisto_low'),'String',x(1));
set(findobj('Tag','input_parthisto_high'),'String',x(2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  histogram set manually                                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_parthisto_manually_Callback(hObject, eventdata, handles)

global storage_av2_stackbrowser;

min=str2num(get(findobj('Tag','input_parthisto_low'),'String'));
max=str2num(get(findobj('Tag','input_parthisto_high'),'String'));

if max>min
    storage_av2_stackbrowser.Scale = [min max];
    set(findobj('Tag','partbrowser_histogram'),'Xlim',[min max]);
else
    set(findobj('Tag','input_parthisto_low'),'String',num2str(storage_av2_stackbrowser.Scale(1)));
    set(findobj('Tag','input_parthisto_high'),'String',num2str(storage_av2_stackbrowser.Scale(2)));
end

render_image();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  bandpass low                                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function partbandpass_low_Callback(hObject, eventdata, handles)

global storage_av2_stackbrowser;

val = str2num(get(hObject,'String'));

if ~isempty(val) & val > 0
    storage_av2_stackbrowser.bandpass.low = round(val);
end

if storage_av2_stackbrowser.bandpass.enable == 1
    render_image();
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  bandpass high                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function partbandpass_high_Callback(hObject, eventdata, handles)

global storage_av2_stackbrowser;

val = str2num(get(hObject,'String'));

if ~isempty(val) & val > 0
    storage_av2_stackbrowser.bandpass.high = round(val);
end

if storage_av2_stackbrowser.bandpass.enable == 1
    render_image();
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  bandpass enable                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkbox_partbandpass_Callback(hObject, eventdata, handles)

global storage_av2_stackbrowser;

storage_av2_stackbrowser.bandpass.enable = get(hObject,'Value');

render_image();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  browse for output filtered stack                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output_stack_browse_Callback(hObject, eventdata, handles)

[filename, pathname] = uiputfile({'*.em'}, 'Pick a stack file');
if ischar(filename)
    set(findobj('Tag','output_stack'),'String',[pathname '/' filename]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  browse for output align file                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_outputalignmentfile_browse_Callback(hObject, eventdata, handles)

[filename, pathname] = uiputfile({'*.mat'}, 'Pick a mat file');
if ischar(filename)
    set(findobj('Tag','output_alignmentfile'),'String',[pathname '/' filename]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  filter stack                                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_filterstack_Callback(hObject, eventdata, handles)

global storage_av2_stackbrowser;

%instack = get(findobj('Tag','input_particlestack'),'String');
%inalign = get(findobj('Tag','input_alignmentfile'),'String');

inalign  =  storage_av2_stackbrowser.alignmentfile;
instack = storage_av2_stackbrowser.stackfile;


outstack = get(findobj('Tag','output_stack'),'String');
outalign = get(findobj('Tag','output_alignmentfile'),'String');
value = str2num(get(findobj('Tag','filterstack_value'),'String'));

contents = get(findobj('Tag','filterstack_fieldname'),'String');
fieldname = contents{get(findobj('Tag','filterstack_fieldname'),'Value')};

contents = get(findobj('Tag','filterstack_operator'),'String');
operator = contents{get(findobj('Tag','filterstack_operator'),'Value')};

if isempty(outstack)
    errordlg('Select output stack file name!');
    return;
end

if isempty(outalign)
    errordlg('Select output alignment file name!');
    return;
end

if isempty(operator)
    errordlg('Input a relational operator!');
    return;
end

if isempty(operator)
    errordlg('Input a relational operator!');
    return;
end

if isempty(fieldname)
    errordlg('Input a fieldname!');
    return;
end

if isempty(value)
    errordlg('Input a value!');
    return;
end

try
    tom_av2_filterstack(instack,inalign,fieldname,operator,value,outstack,outalign);
catch Me
    disp(Me.message);
    errordlg('Filtering failed, please check input!');
end

if (storage_av2_stackbrowser.stack_is_compressed==1)

     [align2d stack_uncompr]=tom_av2_uncompress_stack(outalign, storage_av2_stackbrowser.alignmentfile_org,...
                                                  storage_av2_stackbrowser.stackfile_org,outalign,outstack);
end;
                                              

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  norm particles                                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dropdown_normparticles_Callback(hObject, eventdata, handles)

global storage_av2_stackbrowser;

position = round(get(findobj('Tag','slider_particles'),'Value')).*storage_av2_stackbrowser.numimages+1;
loadstack(position);
render_image();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  make reference stack                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_stackbrowse_makeref_Callback(hObject, eventdata, handles)

global storage_av2_stackbrowser;

if storage_av2_stackbrowser.makerefmode == 1

    instackfile = get(findobj('Tag','input_particlestack'),'String');
    [filename, pathname] = uiputfile({'*.em'}, 'Save reference stack as');
    if ischar(filename)
        outstackfile = [pathname '/' filename];
    else
        return;
    end

    if tom_isemfile(instackfile) ~= 1
        error('Could not load stack file.');
    end

    h = tom_reademheader(instackfile);
    newstack = [];
    c = 1;

    stack_sz=0;
    %check out size of new stack
    for i=1:size(storage_av2_stackbrowser.align,2)
        if storage_av2_stackbrowser.align(i).selected == 1;
            stack_sz=stack_sz+1;
        end
    end

    if stack_sz > 0
        %allocate some memory
        newstack=zeros(h.Header.Size(1),h.Header.Size(2),stack_sz);

        zz=1;
        for i=1:size(storage_av2_stackbrowser.align,2)
            if storage_av2_stackbrowser.align(i).selected == 1
                im = tom_emreadc(instackfile,'subregion',[1 1 i],[h.Header.Size(1)-1 h.Header.Size(2)-1 0]);
                newstack(:,:,zz)=im.Value;
                zz=zz+1;
            end
        end

        tom_emwrite(outstackfile,newstack);
        handles.output = outstackfile;
        guidata(hObject, handles);
        uiresume(handles.figure1);
    else
        error('New stack would contain no particles, please check your input.');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Sort particles                                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_sortparticles_Callback(hObject, eventdata, handles)

global storage_av2_stackbrowser;

fieldnamecontents = get(findobj('Tag','dropdown_partsort_fieldname'),'String');
fieldname = fieldnamecontents{get(findobj('Tag','dropdown_partsort_fieldname'),'Value')};

if strcmp(fieldname,'CCC')
    fieldname = 'ccc';
elseif strcmp(fieldname,'Class')
    fieldname = 'Class';
else
    fieldname = 'partno';
end

fieldordercontents = get(findobj('Tag','dropdown_partsort_order'),'String');
fieldorder = fieldordercontents{get(findobj('Tag','dropdown_partsort_order'),'Value')};

if strcmp(fieldorder,'ASC')
    fieldorder = 'ascend';
else
    fieldorder = 'descend';
end

storage_av2_stackbrowser.align = tom_sortstruct(storage_av2_stackbrowser.align,fieldname, fieldorder);

loadstack(1);
render_image();
set(findobj('Tag','slider_particles'),'Value',0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Sort particles                                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_av2stackbrowser_pca_Callback(hObject, eventdata, handles)

global storage_av2_stackbrowser;

pos = 0;

for i=1:size(storage_av2_stackbrowser.align,2)
    if storage_av2_stackbrowser.align(i).selected == 1;
        pos = i;
        break;
    end
end

if pos == 0;
    errordlg('No class selected');return;
end

pos = ceil(pos ./ 3);

pathstr = fileparts(storage_av2_stackbrowser.alignmentfile);
tom_pcagui([pathstr '/class_' num2str(pos) '/class_' num2str(pos) '.mat']);


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
%%  calculate histogram                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function calc_histogram(im)

global storage_av2_stackbrowser;

if nargin < 1
    [mean max min std] = tom_dev(storage_av2_stackbrowser.stack,'noinfo');
else
    [mean max min std] = tom_dev(im,'noinfo');
end
min = mean-3.*std;
max = mean+3.*std;

[h,n] = tom_hist3d(storage_av2_stackbrowser.stack);
h = 200 .* h ./ (100.*size(storage_av2_stackbrowser.stack,1) .* size(storage_av2_stackbrowser.stack,2) .* size(storage_av2_stackbrowser.stack,3));
axesobj = findobj('Tag','partbrowser_histogram');
axes(axesobj); 
bar(n,h);
axis auto;
set(axesobj,'Tag','partbrowser_histogram');
set(findobj('Tag','input_parthisto_low'),'String',num2str(min));
set(findobj('Tag','input_parthisto_high'),'String',num2str(max));
storage_av2_stackbrowser.Scale = [min max];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  load stack                                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loadstack(start);

global storage_av2_stackbrowser;

if ischar(start)
  if isequal(start, 'next')
    start = storage_av2_stackbrowser.pos + storage_av2_stackbrowser.numimages;
  elseif isequal(start, 'previous')
    start = storage_av2_stackbrowser.pos - storage_av2_stackbrowser.numimages;
  end
end

stop = storage_av2_stackbrowser.numimages;

header = tom_reademheader(storage_av2_stackbrowser.stackfile);
storage_av2_stackbrowser.imsize = header.Header.Size;

if start+stop > size(storage_av2_stackbrowser.align,2)
    missingpages = start+stop-size(storage_av2_stackbrowser.align,2)-1;
    stop = size(storage_av2_stackbrowser.align,2)-start+1;
end

storage_av2_stackbrowser.stack = zeros(header.Header.Size(1),header.Header.Size(2),storage_av2_stackbrowser.numimages);
lauf = 1;

for i=start:start+stop-1
     im = tom_emreadc(storage_av2_stackbrowser.stackfile,'subregion',[1 1 storage_av2_stackbrowser.align(i).partno], [header.Header.Size(1)-1 header.Header.Size(2)-1 0]);
     storage_av2_stackbrowser.stack(:,:,lauf) = im.Value;
     lauf = lauf+1;
end

%if start+stop > header.Header.Size(3)
%    missingpages = start+stop-header.Header.Size(3)-1;
%    stop = header.Header.Size(3)-start+1;
%end
%im = tom_emreadc(storage_av2_stackbrowser.stackfile,'subregion',[1 1 start],[header.Header.Size(1)-1 header.Header.Size(2)-1 stop-1]);
%storage_av2_stackbrowser.stack = im.Value;
storage_av2_stackbrowser.pos = start;

%norm stack
normcontents = get(findobj('Tag','dropdown_normparticles'),'String');
normval = normcontents{get(findobj('Tag','dropdown_normparticles'),'Value')};

if ~strcmp(normval,'none')
    for i=1:size(storage_av2_stackbrowser.stack,3)
        storage_av2_stackbrowser.stack(:,:,i) = tom_norm(storage_av2_stackbrowser.stack(:,:,i),normval);
    end
end

try;calc_histogram(storage_av2_stackbrowser.stack(:,:,lauf-1));catch;calc_histogram(storage_av2_stackbrowser.stack(:,:,end));end;

if storage_av2_stackbrowser.bandpass.high == 0
    set(findobj('Tag','partbandpass_high')','String',num2str(header.Header.Size(1)./2));
    storage_av2_stackbrowser.bandpass.high = header.Header.Size(1)./2;
end

%if exist('missingpages') & missingpages > 0
%    storage_av2_stackbrowser.stack = cat(3,storage_av2_stackbrowser.stack,zeros(header.Header.Size(1),header.Header.Size(2),missingpages));
%end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display image                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function render_image()

global storage_av2_stackbrowser;

%storage_av2_stackbrowser.fontsize = storage_av2_stackbrowser.imsize(1)./storage_av2_stackbrowser.numcols./1.3;

%bandpass filter
if storage_av2_stackbrowser.bandpass.enable == 1
    im = zeros(size(storage_av2_stackbrowser.stack));
    for i=1:size(storage_av2_stackbrowser.stack,3)
        im(:,:,i) = tom_bandpass(storage_av2_stackbrowser.stack(:,:,i),storage_av2_stackbrowser.bandpass.low,storage_av2_stackbrowser.bandpass.high);
    end
    scale = storage_av2_stackbrowser.Scale  ./(storage_av2_stackbrowser.imsize(1) .* storage_av2_stackbrowser.imsize(2));
    if  storage_av2_stackbrowser.bandpass.low > 0
            scale = scale - mean(mean(mean(storage_av2_stackbrowser.stack))) ./ (storage_av2_stackbrowser.imsize(1) .* storage_av2_stackbrowser.imsize(1));
    end
else
    im = storage_av2_stackbrowser.stack;
    scale = storage_av2_stackbrowser.Scale;
end
tmpobj = findobj('Tag','imchoosebox');
axes(tmpobj);


h = tom_dspcub(im, 0, storage_av2_stackbrowser.numcols,scale);
set(h,'buttonDownFcn',@mark_particle);

%draw ccc values
delete(findobj('Tag','cccstring'));
if storage_av2_stackbrowser.display.ccc == 1
    numrows = storage_av2_stackbrowser.numimages ./ storage_av2_stackbrowser.numcols;
    for row=1:numrows
        pos_y = (row-1).*storage_av2_stackbrowser.imsize(2)+storage_av2_stackbrowser.imsize(2).*0.9;
        for column=1:storage_av2_stackbrowser.numcols
            abspartnumber = storage_av2_stackbrowser.page.*storage_av2_stackbrowser.numimages+(row-1).*storage_av2_stackbrowser.numcols+column;
            if abspartnumber <= size(storage_av2_stackbrowser.align,2)
                cccval = storage_av2_stackbrowser.align(abspartnumber).ccc;
                pos_x = (column-1).*storage_av2_stackbrowser.imsize(1)+storage_av2_stackbrowser.imsize(1).*0.65;
                t = text(pos_x,pos_y,sprintf('%0.2f',double(cccval)));
                set(t,'Tag','cccstring','Color',[1 1 0],'FontSize',storage_av2_stackbrowser.fontsize,'FontWeight','bold');
            end
        end
    end
end

%draw quality values
delete(findobj('Tag','qualitystring'));
if storage_av2_stackbrowser.display.quality == 1
    numrows = storage_av2_stackbrowser.numimages ./ storage_av2_stackbrowser.numcols;
    for row=1:numrows
        pos_y = (row-1).*storage_av2_stackbrowser.imsize(2)+storage_av2_stackbrowser.imsize(2).*0.2;
        for column=1:storage_av2_stackbrowser.numcols
            abspartnumber = storage_av2_stackbrowser.page.*storage_av2_stackbrowser.numimages+(row-1).*storage_av2_stackbrowser.numcols+column;
            if abspartnumber <= size(storage_av2_stackbrowser.align,2)
                qualityval = storage_av2_stackbrowser.align(abspartnumber).quality;
                pos_x = (column-1).*storage_av2_stackbrowser.imsize(1)+storage_av2_stackbrowser.imsize(1).*0.8;
                t = text(pos_x,pos_y,sprintf('%i',qualityval));
                set(t,'Tag','qualitystring','Color',[1 1 0],'FontSize',storage_av2_stackbrowser.fontsize,'FontWeight','bold');
            end
        end
    end
end

%draw class numbers
delete(findobj('Tag','classstring'));
if storage_av2_stackbrowser.display.ref_classno == 1
    numrows = storage_av2_stackbrowser.numimages ./ storage_av2_stackbrowser.numcols;
    for row=1:numrows
        pos_y = (row-1).*storage_av2_stackbrowser.imsize(2)+storage_av2_stackbrowser.imsize(2).*0.9;
        for column=1:storage_av2_stackbrowser.numcols
            abspartnumber = storage_av2_stackbrowser.page.*storage_av2_stackbrowser.numimages+(row-1).*storage_av2_stackbrowser.numcols+column;
            if abspartnumber <= size(storage_av2_stackbrowser.align,2)
                classval = storage_av2_stackbrowser.align(abspartnumber).ref_class;
                pos_x = (column-1).*storage_av2_stackbrowser.imsize(1)+storage_av2_stackbrowser.imsize(1).*0.1;
                t = text(pos_x,pos_y,sprintf('%i',classval));
                set(t,'Tag','classstring','Color',[1 1 0],'FontSize',storage_av2_stackbrowser.fontsize,'FontWeight','bold');
            end
        end
    end
end


%draw marks
delete(findobj('-regexp','Tag','rect_*'));
if storage_av2_stackbrowser.display.marks == 1
    numrows = storage_av2_stackbrowser.numimages ./ storage_av2_stackbrowser.numcols;
    for row=1:numrows
        pos_y = (row-1).*storage_av2_stackbrowser.imsize(2)+storage_av2_stackbrowser.imsize(2).*0.1;
        for column=1:storage_av2_stackbrowser.numcols
            abspartnumber = storage_av2_stackbrowser.page.*storage_av2_stackbrowser.numimages+(row-1).*storage_av2_stackbrowser.numcols+column;
            if abspartnumber <= size(storage_av2_stackbrowser.align,2)
                %mark particle green if selected and red if not selected
                if storage_av2_stackbrowser.align(abspartnumber).selected == 1
                    color = [0 1 0];
                else
                    color = [1 0 0];
                end
                pos_x = (column-1).*storage_av2_stackbrowser.imsize(1)+storage_av2_stackbrowser.imsize(1).*0.1;
                r = rectangle('Position',[pos_x,pos_y,storage_av2_stackbrowser.imsize(1).*0.1,storage_av2_stackbrowser.imsize(2).*0.1],'FaceColor',color,'Tag',['rect_' abspartnumber]);
            end
        end
    end
end

set(tmpobj,'Tag','imchoosebox');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  mark particle                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mark_particle(a,b)

global storage_av2_stackbrowser;

point1=get(gca,'currentpoint');
button = get(gcf,'selectiontype');

%button values:
%normal: left mouse button
%alt: right mouse button
%extend: middle mouse buttons

pt = point1(1,1:2);
x1 = round(pt(1));
y1 = round(pt(2));
    
    %determine which particle was picked
    row = ceil(y1./storage_av2_stackbrowser.imsize(2));
    column = ceil(x1./storage_av2_stackbrowser.imsize(1));
    abspartnumber = storage_av2_stackbrowser.page.*storage_av2_stackbrowser.numimages+(row-1).*storage_av2_stackbrowser.numcols+column;

if strcmp(button,'normal') == true

    
    if abspartnumber <= size(storage_av2_stackbrowser.align,2)
        %toggle selection on particle
        if storage_av2_stackbrowser.align(abspartnumber).selected == 1
            storage_av2_stackbrowser.align(abspartnumber).selected = 0;
            color = [1 0 0];
        else
            storage_av2_stackbrowser.align(abspartnumber).selected = 1;
            color = [0 1 0];
        end

        %redraw mark
        delete(findobj('Tag',['rect_' abspartnumber]));
        pos_y = (row-1).*storage_av2_stackbrowser.imsize(2)+storage_av2_stackbrowser.imsize(2).*0.1;
        pos_x = (column-1).*storage_av2_stackbrowser.imsize(1)+storage_av2_stackbrowser.imsize(1).*0.1;
        r = rectangle('Position',[pos_x,pos_y,storage_av2_stackbrowser.imsize(1).*0.1,storage_av2_stackbrowser.imsize(2).*0.1],'FaceColor',color,'Tag',['rect_' abspartnumber]);
    end

elseif strcmp(button,'alt') == true
    
    if (storage_av2_stackbrowser.sn_show_parts==0 && storage_av2_stackbrowser.sh_picks==0)
        
        string = ['Filename: ' storage_av2_stackbrowser.align(abspartnumber).filename];
        string = strvcat(string, ['Position: X ' num2str(storage_av2_stackbrowser.align(abspartnumber).position.x) ', Y ' num2str(storage_av2_stackbrowser.align(abspartnumber).position.y)]);
        string = strvcat(string, ['Class: ' storage_av2_stackbrowser.align(abspartnumber).ref_class]);
        string = strvcat(string, ['Shift: X ' num2str(storage_av2_stackbrowser.align(abspartnumber).shift.x) ', Y ' num2str(storage_av2_stackbrowser.align(abspartnumber).shift.y)]);
        string = strvcat(string, ['Angle: ' num2str(storage_av2_stackbrowser.align(abspartnumber).angle)]);
        string = strvcat(string, ['Isaligned? ', num2str(storage_av2_stackbrowser.align(abspartnumber).isaligned)]);
        try
            string = strvcat(string, ['Dataset: ', storage_av2_stackbrowser.align(abspartnumber).dataset]);
        catch
            string = strvcat(string, ['Dataset: ']);
        end;
        
        string = strvcat(string, ['Number in stack: ', num2str(abspartnumber)]);
        msgbox(string);
    end; 
  
    if (storage_av2_stackbrowser.sn_show_parts==1)
        if (storage_av2_stackbrowser.sh_part_filt_kern~=0)
            tmp_stack=storage_av2_stackbrowser.stack_in_mem.Value(:,:,storage_av2_stackbrowser.align(abspartnumber).dataset);
            for i=1:size(tmp_stack,3)
                tmp_stack(:,:,i)=tom_filter(tmp_stack(:,:,i),storage_av2_stackbrowser.sh_part_filt_kern);
            end;
             tom_display(tmp_stack);
        else
            tom_display(storage_av2_stackbrowser.stack_in_mem.Value(:,:,storage_av2_stackbrowser.align(abspartnumber).dataset));
        end
    end;
    
    if (storage_av2_stackbrowser.sh_picks==1)
        tmp_pl=storage_av2_stackbrowser.align_org(storage_av2_stackbrowser.align(abspartnumber).dataset);
        tmp_pl(1,1).radius=storage_av2_stackbrowser.align_org(1,1).radius;
        tom_av2_particlepickergui(tmp_pl,storage_av2_stackbrowser.sh_part_filt_kern);
        
    end;
    
    
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Unused Callbacks/Create functions                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_particlestack_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_alignmentfile_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function slider_particles_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function input_quality_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_parthisto_high_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_parthisto_low_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_parthisto_low_Callback(hObject, eventdata, handles)
function input_parthisto_high_Callback(hObject, eventdata, handles)
function input_norows_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_nocols_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function partbandpass_high_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function partbandpass_low_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function output_filteredstack_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function filterstack_fieldname_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function filterstack_operator_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function filterstack_value_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function filterstack_operator_Callback(hObject, eventdata, handles)
function filterstack_fieldname_Callback(hObject, eventdata, handles)
function filterstack_value_Callback(hObject, eventdata, handles)
function output_stack_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function output_alignmentfile_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function output_alignmentfile_Callback(hObject, eventdata, handles)
function output_stack_Callback(hObject, eventdata, handles)
function dropdown_normparticles_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function dropdown_partsort_fieldname_Callback(hObject, eventdata, handles)
function dropdown_partsort_fieldname_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function dropdown_partsort_order_Callback(hObject, eventdata, handles)
function dropdown_partsort_order_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_reset_quality.
function checkbox_reset_quality_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_reset_quality (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_reset_quality


% --- Executes on button press in chk_show_parts.
function chk_show_parts_Callback(hObject, eventdata, handles)
% hObject    handle to chk_show_parts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chk_show_parts


global storage_av2_stackbrowser



storage_av2_stackbrowser.sn_show_parts=get(handles.chk_show_parts,'Value');


% --- Executes on button press in chk_incre_load.
function chk_incre_load_Callback(hObject, eventdata, handles)
% hObject    handle to chk_incre_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chk_incre_load

global storage_av2_stackbrowser;
storage_av2_stackbrowser.incre_load=get(handles.chk_show_parts,'Value');


% --- Executes on button press in chk_inv_all.
function chk_inv_all_Callback(hObject, eventdata, handles)
% hObject    handle to chk_inv_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chk_inv_all



function sh_part_filt_kern_Callback(hObject, eventdata, handles)
% hObject    handle to sh_part_filt_kern (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sh_part_filt_kern as text
%        str2double(get(hObject,'String')) returns contents of sh_part_filt_kern as a double
global storage_av2_stackbrowser;
storage_av2_stackbrowser.sh_part_filt_kern=str2double(get(handles.sh_part_filt_kern,'String'));

% --- Executes during object creation, after setting all properties.
function sh_part_filt_kern_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sh_part_filt_kern (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in chk_show_picks.
function chk_show_picks_Callback(hObject, eventdata, handles)
% hObject    handle to chk_show_picks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chk_show_picks

global storage_av2_stackbrowser;
storage_av2_stackbrowser.sh_picks=get(handles.chk_show_picks,'Value');


% --- Executes on button press in button_reset_sel.
function button_reset_sel_Callback(hObject, eventdata, handles)
% hObject    handle to button_reset_sel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA
global storage_av2_stackbrowser;


asw=questdlg('Really Reset selection ...prev. selection will be lost!');

if (strcmp(asw,'No') || strcmp(asw,'Cancel'))
    return;
end;

if (strcmp(asw,'Yes'))
    for i=1:length(storage_av2_stackbrowser.align)
        storage_av2_stackbrowser.align(i).selected = 0;
    end;
    render_image();
end;






