function varargout = tom_av3_oscargui(varargin)
%TOM_AV3_OSCARGUI is a GUI for oscar
%
%   varargout = tom_av3_oscargui(varargin)
%
%This GUI facilitates the use of oscar
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
%   tom_amira_createisosurface(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   TOM_AV3_PASTE4OSCARGUI
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
                   'gui_OpeningFcn', @tom_av3_oscargui_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_av3_oscargui_OutputFcn, ...
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
function tom_av3_oscargui_OpeningFcn(hObject, eventdata, handles, varargin)

%due to problems accessing and updating the handles structure in some
%callbacks this global is used for application data throughout this code.
global storage_oscar;

storage_oscar.volume.filename = '';
storage_oscar.volume.infostring = '';
storage_oscar.volume.size.x = 0;
storage_oscar.volume.size.y = 0;
storage_oscar.volume.size.z = 0;
storage_oscar.volume.thumbnail.image = [];
storage_oscar.volume.thumbnail.DataScale = [];
storage_oscar.volume.objpixelsize = 0;

storage_oscar.template.filename = '';
storage_oscar.template.infostring = '';
storage_oscar.template.size.x = 0;
storage_oscar.template.size.y = 0;
storage_oscar.template.size.z = 0;
storage_oscar.template.thumbnail.image = [];
storage_oscar.template.thumbnail.DataScale = [];
storage_oscar.template.objpixelsize = 0;

storage_oscar.mask.filename = '';
storage_oscar.mask.infostring = '';
storage_oscar.mask.size.x = 0;
storage_oscar.mask.size.y = 0;
storage_oscar.mask.size.z = 0;
storage_oscar.mask.thumbnail.image = [];
storage_oscar.mask.thumbnail.DataScale = [];

storage_oscar.psf.filename = '';
storage_oscar.psf.infostring = '';
storage_oscar.psf.size.x = 0;
storage_oscar.psf.size.y = 0;
storage_oscar.psf.size.z = 0;
storage_oscar.psf.thumbnail.image = [];
storage_oscar.psf.thumbnail.DataScale = [];

storage_oscar.oscar.phi.start = 0;
storage_oscar.oscar.phi.stop = 360;
storage_oscar.oscar.phi.step = 30;
storage_oscar.oscar.psi.start = 0;
storage_oscar.oscar.psi.stop = 360;
storage_oscar.oscar.psi.step = 30;
storage_oscar.oscar.theta.start = 0;
storage_oscar.oscar.theta.stop = 180;
storage_oscar.oscar.theta.step = 30;

storage_oscar.oscar.pe = 16;
storage_oscar.oscar.fftsize = 256;

storage_oscar.output.filename = '';
storage_oscar.output.logfile = '';

storage_oscar.alignlist = '';
storage_oscar.generatescript = 0;

storage_oscar.iterativemode.enable = 0;
storage_oscar.iterativemode.templaterecursive = 0;
storage_oscar.iterativemode.nrrefinements = 3;
storage_oscar.iterativemode.maxiterations = 10;
storage_oscar.iterativemode.convcrit = 0.005;
storage_oscar.iterativemode.datadir = '';
storage_oscar.iterativemode.extractionmask = '';
storage_oscar.iterativemode.pastemask = '';

%setappdata(0,'UseNativeSystemDialogs',0);

storage_oscar.lam.cpus = check_lamnodes();
set(findobj('Tag','pe'),'String',num2str(storage_oscar.lam.cpus));

handles.output = hObject;
guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Output Function                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = tom_av3_oscargui_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Exit                                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_exit_Callback(hObject, eventdata, handles)

storage_oscar = [];
delete(findobj('Tag','tom_av3_oscargui'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                      Volume Callbacks                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Show file dialog to select volume                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_browsevolume_Callback(hObject, eventdata, handles)

global storage_oscar;

[filename, pathname] = uigetfile({'*.vol';'*.em'}, 'Pick an EM-file');
if ischar(filename)
    storage_oscar.volume.filename = [pathname, filename];
    set(findobj('Tag','input_volume'),'String',storage_oscar.volume.filename);

    if tom_isemfile(storage_oscar.volume.filename) ~= 1
        errordlg('File is not an EM file.','File Error');
    else
        header = tom_reademheader(storage_oscar.volume.filename);
        storage_oscar.volume.size.x = header.Header.Size(1);
        storage_oscar.volume.size.y = header.Header.Size(2);
        storage_oscar.volume.size.z = header.Header.Size(3);
        storage_oscar.volume.objpixelsize = header.Header.Objectpixelsize;

        %Update fftsize dropdown box
        fftval = get(findobj('Tag','fftsize'),'Value');
        fftsize = 64;
        boxstring = [];
        while fftsize <= min([storage_oscar.volume.size.x,storage_oscar.volume.size.y,storage_oscar.volume.size.z,256])
            boxstring = strvcat(boxstring,num2str(fftsize));
            fftsize = fftsize*2;
        end
        if fftval > size(boxstring,1)
            set(findobj('Tag','fftsize'),'Value',1);
            fftval = 1;
        end
        set(findobj('Tag','fftsize'),'String',boxstring);
        storage_oscar.oscar.fftsize = str2num(boxstring(fftval,:));
        
        %Display volume info
        storage_oscar.volume.infostring = strvcat(['Size: ', num2str(header.Header.Size(1)), ' x ', num2str(header.Header.Size(2)), ' x ', num2str(header.Header.Size(3))],['object pixel size: ', num2str(header.Header.Objectpixelsize),' A']);
        set(findobj('Tag','volume_info'),'String',storage_oscar.volume.infostring);

        %Display volume thumbnail
        im = tom_emreadc(storage_oscar.volume.filename,'subregion',[1 1 round(storage_oscar.volume.size.z./2)],[storage_oscar.volume.size.x-1 storage_oscar.volume.size.y-1 0]);
        storage_oscar.volume.thumbnail.image = im.Value;
        [mean max min2 std] = tom_dev(storage_oscar.volume.thumbnail.image,'noinfo');
        storage_oscar.volume.thumbnail.DataScale = [mean-2.*std, mean+2.*std];
        axesobj = findobj('Tag','volume_preview');
        axes(axesobj);
        imagesc(storage_oscar.volume.thumbnail.image',storage_oscar.volume.thumbnail.DataScale);
        axis off;colormap gray;axis ij;
        set(axesobj,'Tag','volume_preview');       
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Directly edit volume                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_volume_Callback(hObject, eventdata, handles)

global storage_oscar;

filename = get(hObject, 'String');

if ischar(filename)
    if tom_isemfile(filename) ~= 1
        errordlg('File is not an EM file.','File Error');
    else
        storage_oscar.volume.filename = filename;
        header = tom_reademheader(storage_oscar.volume.filename);
        storage_oscar.volume.size.x = header.Header.Size(1);
        storage_oscar.volume.size.y = header.Header.Size(2);
        storage_oscar.volume.size.z = header.Header.Size(3);
        storage_oscar.volume.objpixelsize = header.Header.Objectpixelsize;

        %Display volume info
        storage_oscar.volume.infostring = strvcat(['Size: ', num2str(header.Header.Size(1)), ' x ', num2str(header.Header.Size(2)), ' x ', num2str(header.Header.Size(3))],['object pixel size: ', num2str(header.Header.Objectpixelsize),' A']);
        set(findobj('Tag','volume_info'),'String',storage_oscar.volume.infostring);

        %Update fftsize dropdown box
        fftval = get(findobj('Tag','fftsize'),'Value');
        fftsize = 64;
        boxstring = [];
        while fftsize <= min([storage_oscar.volume.size.x,storage_oscar.volume.size.y,storage_oscar.volume.size.z,256])
            boxstring = strvcat(boxstring,num2str(fftsize));
            fftsize = fftsize*2;
        end
        if fftval > size(boxstring,1)
            set(findobj('Tag','fftsize'),'Value',1);
            fftval = 1;
        end
        set(findobj('Tag','fftsize'),'String',boxstring);
        storage_oscar.oscar.fftsize = str2num(boxstring(fftval,:));
        
        %Display volume thumbnail
        im = tom_emreadc(storage_oscar.volume.filename,'subregion',[1 1 round(storage_oscar.volume.size.z./2)],[storage_oscar.volume.size.x-1 storage_oscar.volume.size.y-1 0]);
        %storage_oscar.volume.thumbnail.image = tom_bin(im.Value,3);
        [mean max min2 std] = tom_dev(storage_oscar.volume.thumbnail.image,'noinfo');
        storage_oscar.volume.thumbnail.DataScale = [mean-2.*std, mean+2.*std];
        axesobj = findobj('Tag','volume_preview');
        axes(axesobj);
        imagesc(storage_oscar.volume.thumbnail.image',storage_oscar.volume.thumbnail.DataScale);
        axis off;colormap gray;axis ij;
        set(axesobj,'Tag','volume_preview');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Show volume in tom_dspcub                                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_volume_dspcub_Callback(hObject, eventdata, handles)

global storage_oscar;

if ~isempty(storage_oscar.volume.filename)
    vol = tom_emreadc(storage_oscar.volume.filename);
    figure;tom_dspcub(vol.Value);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Show volume in tom_volxyz                                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_volume_xyz_Callback(hObject, eventdata, handles)

global storage_oscar;

if ~isempty(storage_oscar.volume.filename)
    vol = tom_emreadc(storage_oscar.volume.filename);
    tom_volxyz(vol);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Generate volume with tom_av3_paste4oscar                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_volumegenerate_Callback(hObject, eventdata, handles)

global storage_oscar;

[filename, alignlist] = tom_av3_paste4oscargui();

if isempty(filename)
    return;
end

storage_oscar.alignlist = alignlist;
storage_oscar.volume.filename = filename;

set(findobj('Tag','input_volume'),'String',storage_oscar.volume.filename);

if tom_isemfile(storage_oscar.volume.filename) ~= 1
    errordlg('File is not an EM file.','File Error');
else
    header = tom_reademheader(storage_oscar.volume.filename);
    storage_oscar.volume.size.x = header.Header.Size(1);
    storage_oscar.volume.size.y = header.Header.Size(2);
    storage_oscar.volume.size.z = header.Header.Size(3);
    
    %Display volume info
    storage_oscar.volume.infostring = strvcat(['Size: ', num2str(header.Header.Size(1)), ' x ', num2str(header.Header.Size(2)), ' x ', num2str(header.Header.Size(3))],['object pixel size: ', num2str(header.Header.Objectpixelsize),' A']);
    set(findobj('Tag','volume_info'),'String',storage_oscar.volume.infostring);

    %Update fftsize dropdown box
    fftval = get(findobj('Tag','fftsize'),'Value');
    fftsize = 64;
    boxstring = [];
    while fftsize <= min([storage_oscar.volume.size.x,storage_oscar.volume.size.y,storage_oscar.volume.size.z,256])
        boxstring = strvcat(boxstring,num2str(fftsize));
        fftsize = fftsize*2;
    end
    if fftval > size(boxstring,1)
        set(findobj('Tag','fftsize'),'Value',1);
        fftval = 1;
    end
    set(findobj('Tag','fftsize'),'String',boxstring);
    storage_oscar.oscar.fftsize = str2num(boxstring(fftval,:));
    
    %Display volume thumbnail
    im = tom_emreadc(storage_oscar.volume.filename,'subregion',[1 1 round(storage_oscar.volume.size.z./2)],[storage_oscar.volume.size.x-1 storage_oscar.volume.size.y-1 0]);
    storage_oscar.volume.thumbnail.image = tom_bin(im.Value,3);
    [mean max min2 std] = tom_dev(storage_oscar.volume.thumbnail.image,'noinfo');
    storage_oscar.volume.thumbnail.DataScale = [mean-2.*std, mean+2.*std];
    axesobj = findobj('Tag','volume_preview');
    axes(axesobj);
    imagesc(storage_oscar.volume.thumbnail.image',storage_oscar.volume.thumbnail.DataScale);
    axis off;colormap gray;axis ij;
    set(axesobj,'Tag','volume_preview');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                      Template Callbacks                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Show file dialog to select template                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_browsetemplate_Callback(hObject, eventdata, handles)

global storage_oscar;

[filename, pathname] = uigetfile({'*.vol';'*.em'}, 'Pick an EM-file');
if ischar(filename)
    storage_oscar.template.filename = [pathname, filename];
    set(findobj('Tag','input_template'),'String',storage_oscar.template.filename);

    if tom_isemfile(storage_oscar.template.filename) ~= 1
        errordlg('File is not an EM file.','File Error');
    else
        header = tom_reademheader(storage_oscar.template.filename);
        storage_oscar.template.size.x = header.Header.Size(1);
        storage_oscar.template.size.y = header.Header.Size(2);
        storage_oscar.template.size.z = header.Header.Size(3);
        storage_oscar.template.objpixelsize = header.Header.Objectpixelsize;

        %Display template info
        storage_oscar.template.infostring = strvcat(['Size: ', num2str(header.Header.Size(1)), ' x ', num2str(header.Header.Size(2)), ' x ', num2str(header.Header.Size(3))],['object pixel size: ', num2str(header.Header.Objectpixelsize),' A']);
        set(findobj('Tag','template_info'),'String',storage_oscar.template.infostring);

        %Display template thumbnail
        im2 = tom_emreadc(storage_oscar.template.filename,'subregion',[1 1 round(storage_oscar.template.size.z./2)-5],[storage_oscar.template.size.x-1 storage_oscar.template.size.y-1 10]);
        im = 0;
        for i=1:10
            im = im + im2.Value(:,:,i);
        end
        storage_oscar.template.thumbnail.image = im;
        %storage_oscar.template.thumbnail.image = tom_bin(im.Value,1);
        [mean max min std] = tom_dev(storage_oscar.template.thumbnail.image,'noinfo');
        storage_oscar.template.thumbnail.DataScale = [mean-2.*std, mean+2.*std];
        axesobj = findobj('Tag','template_preview');
        axes(axesobj);
        tom_imagesc(storage_oscar.template.thumbnail.image,'noinfo');
        %imagesc(storage_oscar.template.thumbnail.image',storage_oscar.template.thumbnail.DataScale);
        axis off;colormap gray;axis ij;
        set(axesobj,'Tag','template_preview');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Directly edit template                                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_template_Callback(hObject, eventdata, handles)

global storage_oscar;

filename = get(hObject, 'String');

if ischar(filename)
    if tom_isemfile(filename) ~= 1
        errordlg('File is not an EM file.','File Error');
    else
        storage_oscar.template.filename = filename;
        header = tom_reademheader(storage_oscar.template.filename);
        storage_oscar.template.size.x = header.Header.Size(1);
        storage_oscar.template.size.y = header.Header.Size(2);
        storage_oscar.template.size.z = header.Header.Size(3);
        storage_oscar.template.objpixelsize = header.Header.Objectpixelsize;

        %Display template info
        storage_oscar.template.infostring = strvcat(['Size: ', num2str(header.Header.Size(1)), ' x ', num2str(header.Header.Size(2)), ' x ', num2str(header.Header.Size(3))],['object pixel size: ', num2str(header.Header.Objectpixelsize),' A']);
        set(findobj('Tag','template_info'),'String',storage_oscar.template.infostring);

        %Display template thumbnail
        %im = tom_emreadc(storage_oscar.template.filename,'subregion',[1 1 round(storage_oscar.template.size.z./2)],[storage_oscar.template.size.x-1 storage_oscar.template.size.y-1 0]);
        im2 = tom_emreadc(storage_oscar.template.filename,'subregion',[1 1 round(storage_oscar.template.size.z./2)-5],[storage_oscar.template.size.x-1 storage_oscar.template.size.y-1 10]);
        im = 0;
        for i=1:10
            im = im + im2.Value(:,:,i);
        end
        storage_oscar.template.thumbnail.image = im;
        %storage_oscar.template.thumbnail.image = tom_bin(im.Value,1);
        [mean max min std] = tom_dev(storage_oscar.template.thumbnail.image,'noinfo');
        storage_oscar.template.thumbnail.DataScale = [mean-2.*std, mean+2.*std];
        axesobj = findobj('Tag','template_preview');
        axes(axesobj);
        imagesc(storage_oscar.template.thumbnail.image',storage_oscar.template.thumbnail.DataScale);
        axis off;colormap gray;axis ij;
        set(axesobj,'Tag','template_preview');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  template resize                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_templaterescale_Callback(hObject, eventdata, handles)

global storage_oscar;

if isempty(storage_oscar.template.filename)
    errordlg('Select template first!','Error');
    return;
end

prompt = {'Enter old object pixel size:','Enter new object pixel size:'};
dlg_title = 'Rescale template';
num_lines = 1;


def = {num2str(storage_oscar.template.objpixelsize),num2str(storage_oscar.volume.objpixelsize)};
answer = inputdlg(prompt,dlg_title,num_lines,def,'on');

if ~isempty(answer)
    [filename, pathname] = uiputfile({'*.vol','*.em'}, 'Save template as');
    if ischar(filename)
        scalefactor = str2num(answer{1})./str2num(answer{2});
        newsize = [round(storage_oscar.template.size.x.*scalefactor), round(storage_oscar.template.size.y.*scalefactor), round(storage_oscar.template.size.z.*scalefactor)];
        storage_oscar.template.size.x = newsize(1);
        storage_oscar.template.size.y = newsize(2);
        storage_oscar.template.size.z = newsize(3);
        
        temp = tom_emreadc(storage_oscar.template.filename);
        temp=tom_rescale3d(temp.Value,newsize,'bicubic',1);
        
        storage_oscar.template.objpixelsize = str2num(answer{2});
        storage_oscar.template.filename = [pathname,filename];
        tom_emwrite(storage_oscar.template.filename, temp);
        
        header = tom_reademheader(storage_oscar.template.filename);
        header.Header.Objectpixelsize = storage_oscar.template.objpixelsize;
        tom_writeemheader(storage_oscar.template.filename, header.Header);
        
        set(findobj('Tag','input_template'),'String',storage_oscar.template.filename);
        
        %Display template info
        storage_oscar.template.infostring = strvcat(['Size: ', num2str(storage_oscar.template.size.x), ' x ', num2str(storage_oscar.template.size.y), ' x ', num2str(storage_oscar.template.size.z)],['object pixel size: ', num2str(storage_oscar.template.objpixelsize),' A']);
        set(findobj('Tag','template_info'),'String',storage_oscar.template.infostring);

        %Display template thumbnail
        im = tom_emreadc(storage_oscar.template.filename,'subregion',[1 1 round(storage_oscar.template.size.z./2)],[storage_oscar.template.size.x-1 storage_oscar.template.size.y-1 0]);
        storage_oscar.template.thumbnail.image = tom_bin(im.Value,1);
        [mean max min std] = tom_dev(storage_oscar.template.thumbnail.image,'noinfo');
        storage_oscar.template.thumbnail.DataScale = [mean-2.*std, mean+2.*std];
        axesobj = findobj('Tag','template_preview');
        axes(axesobj);
        imagesc(storage_oscar.template.thumbnail.image',storage_oscar.template.thumbnail.DataScale);
        axis off;colormap gray;axis ij;
        set(axesobj,'Tag','template_preview');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  template phase normalization                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_normphase_Callback(hObject, eventdata, handles)

global storage_oscar;

if isempty(storage_oscar.template.filename)
    errordlg('Select template first!','Error');
    return;
end

[filename, pathname] = uiputfile({'*.vol','*.em'}, 'Save template as');
if ischar(filename)
    temp = tom_emreadc(storage_oscar.template.filename);
    temp = tom_norm(temp.Value,'phase');
    storage_oscar.template.filename = [pathname,filename];
    tom_emwrite(storage_oscar.template.filename, temp);
    header = tom_reademheader(storage_oscar.template.filename);
    header.Header.Objectpixelsize = storage_oscar.template.objpixelsize;
    tom_writeemheader(storage_oscar.template.filename, header.Header);
    set(findobj('Tag','input_template'),'String',storage_oscar.template.filename);
    
    %Display template info
    storage_oscar.template.infostring = strvcat(['Size: ', num2str(storage_oscar.template.size.x), ' x ', num2str(storage_oscar.template.size.y), ' x ', num2str(storage_oscar.template.size.z)],['object pixel size: ', num2str(storage_oscar.template.objpixelsize),' A']);
    set(findobj('Tag','template_info'),'String',storage_oscar.template.infostring);

    %Display template thumbnail
    im = tom_emreadc(storage_oscar.template.filename,'subregion',[1 1 round(storage_oscar.template.size.z./2)],[storage_oscar.template.size.x-1 storage_oscar.template.size.y-1 0]);
    storage_oscar.template.thumbnail.image = tom_bin(im.Value,1);
    [mean max min std] = tom_dev(storage_oscar.template.thumbnail.image,'noinfo');
    storage_oscar.template.thumbnail.DataScale = [mean-2.*std, mean+2.*std];
    axesobj = findobj('Tag','template_preview');
    axes(axesobj);
    imagesc(storage_oscar.template.thumbnail.image',storage_oscar.template.thumbnail.DataScale);
    axis off;colormap gray;axis ij;
    set(axesobj,'Tag','template_preview');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Show template in tom_dspcub                                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_template_dspcub_Callback(hObject, eventdata, handles)

global storage_oscar;

if ~isempty(storage_oscar.template.filename)
    vol = tom_emreadc(storage_oscar.template.filename);
    figure;tom_dspcub(vol.Value);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Show template in tom_volxyz                                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_template_xyz_Callback(hObject, eventdata, handles)

global storage_oscar;

if ~isempty(storage_oscar.template.filename)
    vol = tom_emreadc(storage_oscar.template.filename);
    tom_volxyz(vol);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                      Mask Callbacks                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Show file dialog to select mask                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_browsemask_Callback(hObject, eventdata, handles)

global storage_oscar;

[filename, pathname] = uigetfile({'*.vol';'*.em'}, 'Pick an EM-file');
if ischar(filename)
    storage_oscar.mask.filename = [pathname, filename];
    set(findobj('Tag','input_mask'),'String',storage_oscar.mask.filename);

    if tom_isemfile(storage_oscar.mask.filename) ~= 1
        errordlg('File is not an EM file.','File Error');
    else
        header = tom_reademheader(storage_oscar.mask.filename);
        storage_oscar.mask.size.x = header.Header.Size(1);
        storage_oscar.mask.size.y = header.Header.Size(2);
        storage_oscar.mask.size.z = header.Header.Size(3);

        %Display mask info
        storage_oscar.mask.infostring = ['Size: ', num2str(header.Header.Size(1)), ' x ', num2str(header.Header.Size(2)), ' x ', num2str(header.Header.Size(3))];
        set(findobj('Tag','mask_info'),'String',storage_oscar.mask.infostring);

        %Display mask thumbnail
        %im = tom_emreadc(storage_oscar.mask.filename,'subregion',[1 1 round(storage_oscar.mask.size.z./2)],[storage_oscar.mask.size.x-1 storage_oscar.mask.size.y-1 0]);
        im2 = tom_emreadc(storage_oscar.mask.filename,'subregion',[1 1 round(storage_oscar.mask.size.z./2)-5],[storage_oscar.mask.size.x-1 storage_oscar.mask.size.y-1 10]);
        im = 0;
        for i=1:10
            im = im + im2.Value(:,:,i);
        end
        storage_oscar.mask.thumbnail.image = im;
        %storage_oscar.mask.thumbnail.image = tom_bin(im.Value,1);
        %a = tom_norm(storage_oscar.template.thumbnail.image,'phase');
        [mean max min std] = tom_dev(storage_oscar.template.thumbnail.image,'noinfo');

        storage_oscar.mask.thumbnail.DataScale = [mean-3*std, mean+3*std];

        try
            storage_oscar.mask.thumbnail.image = storage_oscar.mask.thumbnail.image.*a;
        end


        axesobj = findobj('Tag','mask_preview');
        axes(axesobj);
        tom_imagesc(storage_oscar.mask.thumbnail.image,'noinfo');
        %imagesc(storage_oscar.mask.thumbnail.image',storage_oscar.mask.thumbnail.DataScale);
        axis off;colormap gray;axis ij;
        set(axesobj,'Tag','mask_preview');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Directly edit maskfile                                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_mask_Callback(hObject, eventdata, handles)

global storage_oscar;

filename = get(hObject, 'String');

if ischar(filename)
    if tom_isemfile(filename) ~= 1
        errordlg('File is not an EM file.','File Error');
    else
        storage_oscar.mask.filename = filename;
        set(findobj('Tag','input_mask'),'String',storage_oscar.mask.filename);
        header = tom_reademheader(storage_oscar.mask.filename);
        storage_oscar.mask.size.x = header.Header.Size(1);
        storage_oscar.mask.size.y = header.Header.Size(2);
        storage_oscar.mask.size.z = header.Header.Size(3);

        %Display mask info
        storage_oscar.mask.infostring = ['Size: ', num2str(header.Header.Size(1)), ' x ', num2str(header.Header.Size(2)), ' x ', num2str(header.Header.Size(3))];
        set(findobj('Tag','mask_info'),'String',storage_oscar.mask.infostring);

        %Display mask thumbnail
        %im = tom_emreadc(storage_oscar.mask.filename,'subregion',[1 1 round(storage_oscar.mask.size.z./2)],[storage_oscar.mask.size.x-1 storage_oscar.mask.size.y-1 0]);
        im2 = tom_emreadc(storage_oscar.mask.filename,'subregion',[1 1 round(storage_oscar.mask.size.z./2)-5],[storage_oscar.mask.size.x-1 storage_oscar.mask.size.y-1 10]);
        im = 0;
        for i=1:10
            im = im + im2.Value(:,:,i);
        end
        storage_oscar.mask.thumbnail.image = im;
%        storage_oscar.mask.thumbnail.image = tom_bin(im.Value,1);
        try
            storage_oscar.mask.thumbnail.image = storage_oscar.mask.thumbnail.image.*storage_oscar.template.thumbnail.image;
        end

        [mean max min std] = tom_dev(storage_oscar.mask.thumbnail.image,'noinfo');
        storage_oscar.mask.thumbnail.DataScale = [mean-2.*std, mean+2.*std];
        axesobj = findobj('Tag','mask_preview');
        axes(axesobj);
        imagesc(storage_oscar.mask.thumbnail.image',storage_oscar.mask.thumbnail.DataScale);
        axis off;colormap gray;axis ij;
        set(axesobj,'Tag','mask_preview');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Generate mask                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_generatemask_Callback(hObject, eventdata, handles)

global storage_oscar;

if isempty(storage_oscar.template.filename)
    errordlg('Select template first!','Error');
    return;
end

prompt = {'Enter radius in pixel:','Enter smoothing in pixel (outside of radius)'};
dlg_title = 'Generate mask';
num_lines = 1;
def = {num2str(round(storage_oscar.template.size.x./2)),'1'};
answer = inputdlg(prompt,dlg_title,num_lines,def,'on');

if ~isempty(answer)
    [filename, pathname] = uiputfile({'*.vol','*.em'}, 'Save mask as');
    if ischar(filename)

        storage_oscar.mask.size.x = storage_oscar.template.size.x;
        storage_oscar.mask.size.y = storage_oscar.template.size.y;
        storage_oscar.mask.size.z = storage_oscar.template.size.z;

        mask=ones(storage_oscar.mask.size.x,storage_oscar.mask.size.y,storage_oscar.mask.size.z);
        mask=tom_spheremask(mask,str2num(answer{1}),str2num(answer{2}),[storage_oscar.mask.size.x./2+1,storage_oscar.mask.size.y./2+1,storage_oscar.mask.size.z./2+1]);
        
        storage_oscar.mask.filename = [pathname,filename];
        tom_emwrite(storage_oscar.mask.filename, mask);
        set(findobj('Tag','input_mask'),'String',storage_oscar.mask.filename);
        
        %Display mask info
        storage_oscar.mask.infostring = ['Size: ', num2str(storage_oscar.mask.size.x), ' x ', num2str(storage_oscar.mask.size.y), ' x ', num2str(storage_oscar.mask.size.x)];
        set(findobj('Tag','mask_info'),'String',storage_oscar.mask.infostring);

        %Display mask thumbnail
        im = tom_emreadc(storage_oscar.mask.filename,'subregion',[1 1 round(storage_oscar.mask.size.z./2)],[storage_oscar.mask.size.x-1 storage_oscar.mask.size.y-1 0]);
        storage_oscar.mask.thumbnail.image = tom_bin(im.Value,1);
        [mean max min std] = tom_dev(storage_oscar.mask.thumbnail.image,'noinfo');
        
        try
            storage_oscar.mask.thumbnail.image = storage_oscar.mask.thumbnail.image.*storage_oscar.template.thumbnail.image;
        end
        storage_oscar.mask.thumbnail.DataScale = [mean-2.*std, mean+2.*std];
        axesobj = findobj('Tag','mask_preview');
        axes(axesobj);
        imagesc(storage_oscar.mask.thumbnail.image',storage_oscar.mask.thumbnail.DataScale);
        axis off;colormap gray;axis ij;
        set(axesobj,'Tag','mask_preview');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Show mask in tom_dspcub                                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_mask_dspcub_Callback(hObject, eventdata, handles)

global storage_oscar;

if ~isempty(storage_oscar.mask.filename)
    vol = tom_emreadc(storage_oscar.mask.filename);
    template = tom_emreadc(storage_oscar.template.filename);
    try
        vol.Value = vol.Value*template.Value;
    end

    figure;tom_dspcub(vol.Value);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Show mask in tom_volxyz                                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_mask_xyz_Callback(hObject, eventdata, handles)

global storage_oscar;

if ~isempty(storage_oscar.mask.filename)
    vol = tom_emreadc(storage_oscar.mask.filename);
    template = tom_emreadc(storage_oscar.template.filename);
    try
        vol.Value = vol.Value*template.Value;
    end
    tom_volxyz(vol);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                        PSF Callbacks                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Show file dialog to select PSF                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_browsepsf_Callback(hObject, eventdata, handles)

global storage_oscar;

[filename, pathname] = uigetfile({'*.vol';'*.em'}, 'Pick an EM-file');
if ischar(filename)
    storage_oscar.psf.filename = [pathname, filename];
    set(findobj('Tag','input_psf'),'String',storage_oscar.psf.filename);

    if tom_isemfile(storage_oscar.psf.filename) ~= 1
        errordlg('File is not an EM file.','File Error');
    else
        header = tom_reademheader(storage_oscar.psf.filename);
        storage_oscar.psf.size.x = header.Header.Size(1);
        storage_oscar.psf.size.y = header.Header.Size(2);
        storage_oscar.psf.size.z = header.Header.Size(3);

        %Display psf info
        storage_oscar.psf.infostring = ['Size: ', num2str(header.Header.Size(1)), ' x ', num2str(header.Header.Size(2)), ' x ', num2str(header.Header.Size(3))];
        set(findobj('Tag','psf_info'),'String',storage_oscar.psf.infostring);

        %Display psf thumbnail
        im = tom_emreadc(storage_oscar.psf.filename,'subregion',[1 1 round(storage_oscar.psf.size.z./2)],[storage_oscar.psf.size.x-1 storage_oscar.psf.size.y-1 0]);
        storage_oscar.psf.thumbnail.image = tom_bin(im.Value,1);
        [mean max min std] = tom_dev(storage_oscar.psf.thumbnail.image,'noinfo');
        storage_oscar.psf.thumbnail.DataScale = [mean-2.*std, mean+2.*std];
        axesobj = findobj('Tag','psf_preview');
        axes(axesobj);
        imagesc(storage_oscar.psf.thumbnail.image',storage_oscar.psf.thumbnail.DataScale);
        axis off;colormap gray;axis ij;
        set(axesobj,'Tag','psf_preview');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Directly edit psffile                                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_psf_Callback(hObject, eventdata, handles)

global storage_oscar;

filename = get(hObject, 'String');

if ischar(filename)
   if tom_isemfile(filename) ~= 1
        errordlg('File is not an EM file.','File Error');
   else
        storage_oscar.psf.filename = filename;
        header = tom_reademheader(storage_oscar.psf.filename);
        storage_oscar.psf.size.x = header.Header.Size(1);
        storage_oscar.psf.size.y = header.Header.Size(2);
        storage_oscar.psf.size.z = header.Header.Size(3);

        %Display psf info
        storage_oscar.psf.infostring = ['Size: ', num2str(header.Header.Size(1)), ' x ', num2str(header.Header.Size(2)), ' x ', num2str(header.Header.Size(3))];
        set(findobj('Tag','psf_info'),'String',storage_oscar.psf.infostring);

        %Display psf thumbnail
        im = tom_emreadc(storage_oscar.psf.filename,'subregion',[1 1 round(storage_oscar.psf.size.z./2)],[storage_oscar.psf.size.x-1 storage_oscar.psf.size.y-1 0]);
        storage_oscar.psf.thumbnail.image = tom_bin(im.Value,1);
        [mean max min std] = tom_dev(storage_oscar.psf.thumbnail.image,'noinfo');
        storage_oscar.psf.thumbnail.DataScale = [mean-2.*std, mean+2.*std];
        axesobj = findobj('Tag','psf_preview');
        axes(axesobj);
        imagesc(storage_oscar.psf.thumbnail.image',storage_oscar.psf.thumbnail.DataScale);
        axis off;colormap gray;axis ij;
        set(axesobj,'Tag','psf_preview');
   end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Generate psffile                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_generatepsf_Callback(hObject, eventdata, handles)

global storage_oscar;

if isempty(storage_oscar.template.filename)
    errordlg('Select template first!','Error');
    return;
end

prompt = {'Enter missing wedge in degrees:'};
dlg_title = 'Generate point spread function';
num_lines = 1;
def = {'30'};
answer = inputdlg(prompt,dlg_title,num_lines,def,'on');

if ~isempty(answer)
    [filename, pathname] = uiputfile({'*.vol','*.em'}, 'Save PSF as');
    if ischar(filename)
        yyy = zeros(storage_oscar.template.size.x, storage_oscar.template.size.y, storage_oscar.template.size.z);
        wedge=tom_wedge(yyy,str2num(answer{1}));
        yyy(1,1,1)=1;
        psf = real(tom_ifourier(ifftshift(fftshift(tom_fourier(yyy)).*wedge)));
        storage_oscar.psf.filename = [pathname,filename];
        tom_emwrite(storage_oscar.psf.filename, psf);

        set(findobj('Tag','input_psf'),'String',storage_oscar.psf.filename);
        
        storage_oscar.psf.size.x = storage_oscar.template.size.x;
        storage_oscar.psf.size.y = storage_oscar.template.size.y;
        storage_oscar.psf.size.z = storage_oscar.template.size.z;
        
        %Display psf info
        storage_oscar.psf.infostring = ['Size: ', num2str(storage_oscar.psf.size.x), ' x ', num2str(storage_oscar.psf.size.y), ' x ', num2str(storage_oscar.psf.size.x)];
        set(findobj('Tag','psf_info'),'String',storage_oscar.psf.infostring);

        %Display psf thumbnail
        im = tom_emreadc(storage_oscar.psf.filename,'subregion',[1 1 round(storage_oscar.psf.size.z./2)],[storage_oscar.psf.size.x-1 storage_oscar.psf.size.y-1 0]);
        storage_oscar.psf.thumbnail.image = tom_bin(im.Value,1);
        [mean max min std] = tom_dev(storage_oscar.psf.thumbnail.image,'noinfo');
        storage_oscar.psf.thumbnail.DataScale = [mean-2.*std, mean+2.*std];
        axesobj = findobj('Tag','psf_preview');
        axes(axesobj);
        imagesc(storage_oscar.psf.thumbnail.image',storage_oscar.psf.thumbnail.DataScale);
        axis off;colormap gray;axis ij;
        set(axesobj,'Tag','psf_preview');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Show psf in tom_dspcub                                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_psf_dspcub_Callback(hObject, eventdata, handles)

global storage_oscar;

if ~isempty(storage_oscar.psf.filename)
    vol = (tom_emreadc(storage_oscar.psf.filename));
    figure;tom_dspcub(tom_ps(vol.Value));title('Powerspectrum of PSF')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Show psf in tom_volxyz                                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_psf_xyz_Callback(hObject, eventdata, handles)

global storage_oscar;

if ~isempty(storage_oscar.psf.filename)
    vol = (tom_emreadc(storage_oscar.psf.filename));
    tom_volxyz(tom_ps(vol.Value));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  Oscar Options Callbacks                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Phi start                                                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function phi_start_Callback(hObject, eventdata, handles)

global storage_oscar;

storage_oscar.oscar.phi.start = str2num(get(hObject,'String'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Phi stop                                                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function phi_stop_Callback(hObject, eventdata, handles)

global storage_oscar;

storage_oscar.oscar.phi.stop = str2num(get(hObject,'String'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Phi step                                                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function phi_step_Callback(hObject, eventdata, handles)

global storage_oscar;

storage_oscar.oscar.phi.step = str2num(get(hObject,'String'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Psi start                                                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function psi_start_Callback(hObject, eventdata, handles)

global storage_oscar;

storage_oscar.oscar.psi.start = str2num(get(hObject,'String'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Psi stop                                                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function psi_stop_Callback(hObject, eventdata, handles)

global storage_oscar;

storage_oscar.oscar.psi.stop = str2num(get(hObject,'String'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Psi step                                                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function psi_step_Callback(hObject, eventdata, handles)

global storage_oscar;

storage_oscar.oscar.psi.step = str2num(get(hObject,'String'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Theta start                                                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function theta_start_Callback(hObject, eventdata, handles)

global storage_oscar;

storage_oscar.oscar.theta.start = str2num(get(hObject,'String'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Theta stop                                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function theta_stop_Callback(hObject, eventdata, handles)

global storage_oscar;

storage_oscar.oscar.theta.stop = str2num(get(hObject,'String'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Theta step                                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function theta_step_Callback(hObject, eventdata, handles)

global storage_oscar;

storage_oscar.oscar.theta.step = str2num(get(hObject,'String'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Processing elements                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pe_Callback(hObject, eventdata, handles)

global storage_oscar;

storage_oscar.oscar.pe = round(str2num(get(hObject,'String')));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Set Processing elements to max cpus available from lamnodes        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_peset_Callback(hObject, eventdata, handles)

global storage_oscar;

storage_oscar.lam.cpus = check_lamnodes();
set(findobj('Tag','pe'),'String',num2str(storage_oscar.lam.cpus));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  fftsize                                                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fftsize_Callback(hObject, eventdata, handles)

global storage_oscar;

fftval = get(hObject,'Value');
fftstring = get(hObject,'String');
storage_oscar.oscar.fftsize = str2num(fftstring(fftval,:));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                        Output Callbacks                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Show file dialog to select Outfiles                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function browse_outfiles_Callback(hObject, eventdata, handles)

global storage_oscar;

[filename, pathname] = uiputfile({'*.*'}, 'Pick an Outfile');
if ischar(filename)
    storage_oscar.output.filename = [pathname, filename];
    %storage_oscar.output.filename = storage_oscar.output.filename(1:max(strfind(storage_oscar.output.filename,'.'))-1);
    set(findobj('Tag','output_files'),'String',storage_oscar.output.filename);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Directly edit Outfile                                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output_files_Callback(hObject, eventdata, handles)

global storage_oscar;

filename = get(hObject,'String');

if ischar(filename)
    storage_oscar.output.filename = filename;
    %storage_oscar.output.filename = storage_oscar.output.filename(1:max(strfind(storage_oscar.output.filename,'.'))-1);
    set(findobj('Tag','output_files'),'String',storage_oscar.output.filename);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Show file dialog to select logfile                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function browse_log_Callback(hObject, eventdata, handles)

global storage_oscar;

[filename, pathname] = uiputfile({'*.*'}, 'Pick an Outfile');
if ischar(filename)
    storage_oscar.output.logfile = [pathname, filename];
    set(findobj('Tag','output_log'),'String',storage_oscar.output.logfile);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Directly edit logfile                                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output_log_Callback(hObject, eventdata, handles)

global storage_oscar;

filename = get(hObject,'String');

if ischar(filename)
    storage_oscar.output.logfile = filename;
    set(findobj('Tag','output_log'),'String',storage_oscar.output.logfile);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Show file dialog to load settings                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function buttons_settingsload_Callback(hObject, eventdata, handles)

global storage_oscar;

[filename, pathname] = uigetfile('*.mat', 'Load Current Settings');
if ischar(filename)
    s = load([pathname '/' filename]);

    if isfield(s,'storagesave_oscar') ~= 1
        errordlg('This is not a settings file.','File Error');
    end

    storage_oscar = s.storagesave_oscar;
    
    %Restore Oscar options
    set(findobj('Tag','phi_start'),'String',num2str(storage_oscar.oscar.phi.start));
    set(findobj('Tag','phi_stop'),'String',num2str(storage_oscar.oscar.phi.stop));
    set(findobj('Tag','phi_step'),'String',num2str(storage_oscar.oscar.phi.step));
    set(findobj('Tag','psi_start'),'String',num2str(storage_oscar.oscar.psi.start));
    set(findobj('Tag','psi_stop'),'String',num2str(storage_oscar.oscar.psi.stop));
    set(findobj('Tag','psi_step'),'String',num2str(storage_oscar.oscar.psi.step));
    set(findobj('Tag','theta_start'),'String',num2str(storage_oscar.oscar.theta.start));
    set(findobj('Tag','theta_stop'),'String',num2str(storage_oscar.oscar.theta.stop));
    set(findobj('Tag','theta_step'),'String',num2str(storage_oscar.oscar.theta.step));
    
    fftsize = 64;
    boxstring = [];
    storedfftsize = storage_oscar.oscar.fftsize;
    storedfftcounter = 1;
    while fftsize <= min([storage_oscar.volume.size.x,storage_oscar.volume.size.y,storage_oscar.volume.size.z,256])
        if storedfftsize == fftsize
            selectedfftsize = storedfftcounter;
        end
        storedfftcounter = storedfftcounter + 1;
        boxstring = strvcat(boxstring,num2str(fftsize));
        fftsize = fftsize*2;
    end
    set(findobj('Tag','fftsize'),'String',boxstring,'Value',selectedfftsize);
    

    set(findobj('Tag','pe'),'String',num2str(storage_oscar.oscar.pe));
    
    set(findobj('Tag','checkbox_generatescript'),'Value',storage_oscar.generatescript);
    
    %Restore output files
    set(findobj('Tag','output_files'),'String',storage_oscar.output.filename);
    set(findobj('Tag','output_log'),'String',storage_oscar.output.logfile);
    
    %Restore volume
    set(findobj('Tag','input_volume'),'String',storage_oscar.volume.filename);
    set(findobj('Tag','volume_info'),'String',storage_oscar.volume.infostring);
    axesobj = findobj('Tag','volume_preview');
    axes(axesobj);
    imagesc(storage_oscar.volume.thumbnail.image',storage_oscar.volume.thumbnail.DataScale);
    axis off;colormap gray;axis ij;
    set(axesobj,'Tag','volume_preview');
    
    %Restore template
    set(findobj('Tag','input_template'),'String',storage_oscar.template.filename);
    set(findobj('Tag','template_info'),'String',storage_oscar.template.infostring);
    axesobj = findobj('Tag','template_preview');
    axes(axesobj);
    imagesc(storage_oscar.template.thumbnail.image',storage_oscar.template.thumbnail.DataScale);
    axis off;colormap gray;axis ij;
    set(axesobj,'Tag','template_preview');
    
    %Restore mask
    set(findobj('Tag','input_mask'),'String',storage_oscar.mask.filename);
    set(findobj('Tag','mask_info'),'String',storage_oscar.mask.infostring);
    axesobj = findobj('Tag','mask_preview');
    axes(axesobj);
    imagesc(storage_oscar.mask.thumbnail.image',storage_oscar.mask.thumbnail.DataScale);
    axis off;colormap gray;axis ij;
    set(axesobj,'Tag','mask_preview');

    %Restore psf
    set(findobj('Tag','input_psf'),'String',storage_oscar.psf.filename);
    set(findobj('Tag','psf_info'),'String',storage_oscar.psf.infostring);
    axesobj = findobj('Tag','psf_preview');
    axes(axesobj);
    imagesc(storage_oscar.psf.thumbnail.image',storage_oscar.psf.thumbnail.DataScale);
    axis off;colormap gray;axis ij;
    set(axesobj,'Tag','psf_preview');
    
    %Restore iterative mode box
    set(findobj('Tag','checkbox_runiterativemode'),'Value',storage_oscar.iterativemode.enable);
    set(findobj('Tag','checkbox_recursivetemplate'),'Value',storage_oscar.iterativemode.templaterecursive);
    set(findobj('Tag','output_runtimedata'),'String',storage_oscar.iterativemode.datadir);
    set(findobj('Tag','input_extractionmask'),'String',storage_oscar.iterativemode.extractionmask);
    set(findobj('Tag','input_pastemask'),'String',storage_oscar.iterativemode.pastemask);
    set(findobj('Tag','iterative_number'),'String',num2str(storage_oscar.iterativemode.nrrefinements));
    set(findobj('Tag','iterative_maxiterations'),'String',num2str(storage_oscar.iterativemode.maxiterations));
    set(findobj('Tag','iterative_convergence'),'String',num2str(storage_oscar.iterativemode.convcrit));
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Show file dialog to save settings                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_settingssave_Callback(hObject, eventdata, handles)

global storage_oscar;

[filename, pathname] = uiputfile('*.mat', 'Save Current Settings as');
if ischar(filename)
    storagesave_oscar = storage_oscar;
    save([pathname '/' filename], 'storagesave_oscar');
    disp('Settings Saved');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  iterative mode                                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkbox_runiterativemode_Callback(hObject, eventdata, handles)

global storage_oscar;

storage_oscar.iterativemode.enable = get(hObject,'Value');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  number of iterations                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function iterative_number_Callback(hObject, eventdata, handles)

global storage_oscar;

storage_oscar.iterativemode.nrrefinements = round(str2num(get(hObject,'String')));
set(hObject,'String',num2str(storage_oscar.iterativemode.nrrefinements));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  recursive template                                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkbox_recursivetemplate_Callback(hObject, eventdata, handles)

global storage_oscar;

storage_oscar.iterativemode.templaterecursive = get(hObject,'Value');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  convergence criteria                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function iterative_convergence_Callback(hObject, eventdata, handles)

global storage_oscar;

storage_oscar.iterativemode.convcrit = str2num(get(hObject,'String'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  maxiterations                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function iterative_maxiterations_Callback(hObject, eventdata, handles)

global storage_oscar;

storage_oscar.iterativemode.maxiterations = str2num(get(hObject,'String'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  recursive runtime data directory browse                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function browse_runtimedata_Callback(hObject, eventdata, handles)

global storage_oscar;

dirname = uigetdir('Pick a directory');
if ischar(dirname)
    storage_oscar.iterativemode.datadir = dirname;
    set(findobj('Tag','output_runtimedata'),'String',storage_oscar.iterativemode.datadir);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  recursive runtime data directory edit                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output_runtimedata_Callback(hObject, eventdata, handles)

global storage_oscar;

storage_oscar.iterativemode.datadir = get(hObject,'String');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  extraction mask browse                                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_browseextractionmask_Callback(hObject, eventdata, handles)

global storage_oscar;

[filename, pathname] = uigetfile({'*.em','*.vol'}, 'Pick a extraction mask file');
if ischar(filename)
    storage_oscar.iterativemode.extractionmask = [pathname, filename];
    set(findobj('Tag','input_extractionmask'),'String',storage_oscar.iterativemode.extractionmask);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  directy edit extraction mask                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_extractionmask_Callback(hObject, eventdata, handles)

global storage_oscar;

storage_oscar.iterativemode.extractionmask = get(hObject,'String');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  paste mask browse                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_browsepastemask_Callback(hObject, eventdata, handles)

global storage_oscar;

[filename, pathname] = uigetfile({'*.em','*.vol'}, 'Pick a paste mask file');
if ischar(filename)
    storage_oscar.iterativemode.pastemask = [pathname, filename];
    set(findobj('Tag','input_pastemask'),'String',storage_oscar.iterativemode.pastemask);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  directy edit paste mask                                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_pastemask_Callback(hObject, eventdata, handles)

global storage_oscar;

storage_oscar.iterativemode.pastemask = get(hObject,'String');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  lamboot                                                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_lamboot_Callback(hObject, eventdata, handles)

tom_lamboot;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Run Oscar                                                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_runoscar_Callback(hObject, eventdata, handles)

global storage_oscar;

%%%%%%%%%%%%%%%%%%%%%%%
%Check all the settings
%%%%%%%%%%%%%%%%%%%%%%%

%Check for files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(storage_oscar.volume.filename)
    errordlg('No volume file selected!');
    return;
end 
if isempty(storage_oscar.template.filename)
    errordlg('No template file selected!');
    return;
end 
if isempty(storage_oscar.mask.filename)
    errordlg('No mask file selected!');
    return;
end 
if isempty(storage_oscar.psf.filename)
    errordlg('No psf file selected!');
    return;
end 
if isempty(storage_oscar.output.filename) & storage_oscar.iterativemode.enable == 0
    errordlg('No output file selected!');
    return;
end 
if isempty(storage_oscar.output.logfile) & storage_oscar.iterativemode.enable == 0
    errordlg('No log file selected!');
    return;
end 


%Check for template, mask and psf size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if storage_oscar.mask.size.x ~= storage_oscar.template.size.x | storage_oscar.mask.size.y ~= storage_oscar.template.size.y | storage_oscar.mask.size.z ~= storage_oscar.template.size.z
    errordlg('Template and Mask must have the same dimensions!');
    return;
end

if storage_oscar.psf.size.x ~= storage_oscar.template.size.x | storage_oscar.psf.size.y ~= storage_oscar.template.size.y | storage_oscar.psf.size.z ~= storage_oscar.template.size.z
    errordlg('Template and Point Spread Function must have the same dimensions!');
    return;
end

%Check for angles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isnumeric(storage_oscar.oscar.phi.start) 
    errordlg('Angle Phi start must be a number!','Input Error');
    return;
end
if ~isnumeric(storage_oscar.oscar.phi.stop) 
    errordlg('Angle Phi stop must be a number!','Input Error'); 
    return;
end
if ~isnumeric(storage_oscar.oscar.phi.step) 
    errordlg('Angle Phi step must be a number!','Input Error'); 
    return;
end

if ~isnumeric(storage_oscar.oscar.psi.start) 
    errordlg('Angle Psi start must be a number!','Input Error'); 
    return;
end
if ~isnumeric(storage_oscar.oscar.psi.stop) 
    errordlg('Angle Psi stop must be a number!','Input Error'); 
    return;
end
if ~isnumeric(storage_oscar.oscar.psi.step) 
    errordlg('Angle Psi step must be a number!','Input Error'); 
    return;
end

if ~isnumeric(storage_oscar.oscar.theta.start) 
    errordlg('Angle Theta start must be a number!','Input Error'); 
    return;
end
if ~isnumeric(storage_oscar.oscar.theta.stop)
    errordlg('Angle Theta stop must be a number!','Input Error'); 
    return;
end
if ~isnumeric(storage_oscar.oscar.theta.step) 
    errordlg('Angle Theta step must be a number!','Input Error'); 
    return;
end

%Check for oscar options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isnumeric(storage_oscar.oscar.fftsize)
    errordlg('FFTsize is invalid!','Input Error'); 
    return;
end

if ~isnumeric(storage_oscar.oscar.pe)
    errordlg('number of cpus is invalid!','Input Error'); 
    return;
end

if storage_oscar.oscar.pe == 0
    errordlg('Number of cpus = 0, try lamboot first.','No cpus available');
end


%% All checks passed, run oscar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%noniterative mode
%%%%%%%%%%%%%%%%%%%%
if get(findobj('Tag','checkbox_runiterativemode'),'Value') == 0
    angle_start = [storage_oscar.oscar.phi.start, storage_oscar.oscar.psi.start, storage_oscar.oscar.theta.start];
    angle_end = [storage_oscar.oscar.phi.stop, storage_oscar.oscar.psi.stop, storage_oscar.oscar.theta.stop];
    angle_incr = [storage_oscar.oscar.phi.step, storage_oscar.oscar.psi.step, storage_oscar.oscar.theta.step];

    tom_av3_oscar_feed(storage_oscar.volume.filename, storage_oscar.template.filename, storage_oscar.output.filename, storage_oscar.psf.filename, storage_oscar.mask.filename, storage_oscar.output.logfile, storage_oscar.oscar.fftsize, storage_oscar.oscar.pe, angle_start, angle_incr, angle_end, storage_oscar.generatescript);

%iterative mode
%%%%%%%%%%%%%%%%%%%
else
    %check for mask file
    if isempty(storage_oscar.iterativemode.extractionmask)
        errordlg('Select extraction mask file first!');
        return;
    end
    
    %check for align file
    if isempty(storage_oscar.alignlist)
        [filename, pathname] = uigetfile('*.mat', 'Load Particle Alignment File');
        if ischar(filename)
            storage_oscar.alignlist = [pathname, filename];
        else
            errordlg('A Particle Alignment File must be provided for iterative mode!','Missing Alignment File');
            return;
        end
    end
   
    %execute directly
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    if get(findobj('Tag','checkbox_generatescript'),'Value') == 0
        
        disp(strvcat([''],['**************************'],['Iterative mode activated.'],['**************************'],['']));
        if storage_oscar.iterativemode.templaterecursive == 1
            disp(strvcat(['Will run for ' num2str(storage_oscar.iterativemode.nrrefinements) ' refinement steps, convergence criteria is ' num2str(storage_oscar.iterativemode.convcrit) '.'],[' ']));
        else
            disp(strvcat(['Will run for ' num2str(storage_oscar.iterativemode.nrrefinements) ' refinement steps.'],[' ']));
        end
        %copy files to datadir
        disp(['Copying data to datadir ' storage_oscar.iterativemode.datadir ' ...']);
        unix(['cp ' storage_oscar.volume.filename ' ' storage_oscar.iterativemode.datadir '/run1_volume.vol']);
        unix(['cp ' storage_oscar.psf.filename ' ' storage_oscar.iterativemode.datadir '/psf.vol']);
        unix(['cp ' storage_oscar.mask.filename ' ' storage_oscar.iterativemode.datadir '/mask.vol']);
        unix(['cp ' storage_oscar.iterativemode.extractionmask ' ' storage_oscar.iterativemode.datadir '/extract_mask.vol']);
        unix(['cp ' storage_oscar.iterativemode.pastemask ' ' storage_oscar.iterativemode.datadir '/paste_mask.vol']);
        unix(['cp ' storage_oscar.alignlist ' ' storage_oscar.iterativemode.datadir '/align.mat']);

        %create filenames of subvolumes
        a = load(storage_oscar.alignlist);
        filenames = '';
        for k = 1:size(a.Align,2)
            filenames = strvcat(filenames,a.Align(1,k).Filename);
        end

        %create angles
        angle_start = [storage_oscar.oscar.phi.start, storage_oscar.oscar.psi.start, storage_oscar.oscar.theta.start];
        angle_end = [storage_oscar.oscar.phi.stop, storage_oscar.oscar.psi.stop, storage_oscar.oscar.theta.stop];
        angle_incr = [storage_oscar.oscar.phi.step, storage_oscar.oscar.psi.step, storage_oscar.oscar.theta.step];

        if storage_oscar.iterativemode.templaterecursive == 1
            unix(['cp ' storage_oscar.template.filename ' ' storage_oscar.iterativemode.datadir '/run1_template_1.vol']);
            tmplatestring = [storage_oscar.iterativemode.datadir '/run1_template_1.vol'];
            tmpheader = tom_reademheader([storage_oscar.iterativemode.datadir '/run1_template_1.vol']);
        else
            unix(['cp ' storage_oscar.template.filename ' ' storage_oscar.iterativemode.datadir '/template.vol']);
            tmplatestring = [storage_oscar.iterativemode.datadir '/template.vol'];
        end

        mask = tom_emreadc([storage_oscar.iterativemode.datadir '/extract_mask.vol']);
        mask = mask.Value;
        avgrefinementlist = {};

        %template recursive mode
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if storage_oscar.iterativemode.templaterecursive == 1

            for n = 1:storage_oscar.iterativemode.nrrefinements

                if mod(angle_incr.*1000,10) ~= 0
                    disp('Angular Increment reached precision limit, terminating.');
                    break;
                end

                for it=1:storage_oscar.iterativemode.maxiterations

                    disp(strvcat([' '],['==========================================='],[num2str(n) '. refinement step, ' num2str(it) '. iteration'],['==========================================='],[' ']));
                    disp(['Angles: ' num2str(angle_start(1)) ' ' num2str(angle_end(1)) ' ' num2str(angle_incr(1)) ' ' num2str(angle_start(2)) ' ' num2str(angle_end(2)) ' ' num2str(angle_incr(2)) ' ' num2str(angle_start(3)) ' ' num2str(angle_end(3)) ' ' num2str(angle_incr(3))]);

                    %run oscar
                    status = tom_av3_oscar_feed([storage_oscar.iterativemode.datadir '/run' num2str(n) '_volume.vol'],...
                        tmplatestring,...
                        [storage_oscar.iterativemode.datadir '/run' num2str(n) '_Out'],...
                        [storage_oscar.iterativemode.datadir '/psf.vol'],...
                        [storage_oscar.iterativemode.datadir '/mask.vol'],...
                        [storage_oscar.iterativemode.datadir '/run' num2str(n) '_oscar_' num2str(it) '.log'],...
                        storage_oscar.oscar.fftsize, storage_oscar.oscar.pe, angle_start, angle_incr, angle_end, 0);
                    if status == 0
                        disp('Oscar terminated with problems, exiting.');
                        return;
                    end

                    %extract
                    disp(strvcat([' '],['Extracting data ...']));
                    a = load([storage_oscar.iterativemode.datadir '/align.mat']);
                    Align = tom_av3_extract_anglesshifts([storage_oscar.iterativemode.datadir '/run' num2str(n) '_Out'], a.Align, mask, 2);
                    save([storage_oscar.iterativemode.datadir '/align.mat'], 'Align');

                    %generate average
                    disp(strvcat([' '],['Generating average ...']));
                    al = tom_av3_align_sum(Align);
                    average=tom_av3_average(al,'sum',0,0,2);
                    average = tom_emheader(average);
                    tempheader = tom_reademheader(al(1,1).Filename);
                    average.Header.Objectpixelsize = tempheader.Header.Objectpixelsize;
                    tom_emwrite([storage_oscar.iterativemode.datadir '/run' num2str(n) '_avg_' num2str(it) '.vol'],average);
                    avgrefinementlist{n} = [storage_oscar.iterativemode.datadir '/run' num2str(n) '_avg_' num2str(it) '.vol'];

                    %calculate derivative of templates from this run and the
                    %last two runs, if it is smaller than the convergence
                    %criteria, refine the angular range and start again.
                    moveon = 0;
                    if it > 2
                        last_1=tom_emreadc(['run' num2str(n) '_avg_' num2str(it-1) '.vol']);
                        last_2=tom_emreadc(['run' num2str(n) '_avg_' num2str(it-2) '.vol']);
                        %average_ccc(1)=tom_ccc(last_1.Value .* mask,last_2.Value .* mask,'norm');
                        %average_ccc(2)=tom_ccc(average.Value .* mask,last_1.Value .* mask,'norm');
                        average_ccc(1)=tom_ccc(last_1.Value,last_2.Value,'norm');
                        average_ccc(2)=tom_ccc(average.Value,last_1.Value,'norm');
                        disp(strvcat([' '],['absolute ccc values of last three averages: ', num2str(average_ccc(1)), ', ' num2str(average_ccc(2))]));
                        diff_average_ccc=diff(average_ccc);
                        disp(strvcat([' '],['*** differential average ccc: ' num2str(diff_average_ccc), ' ***'],[' ']));
                        if abs(diff_average_ccc) < storage_oscar.iterativemode.convcrit | average_ccc(2) > 1-storage_oscar.iterativemode.convcrit
                            moveon = 1;
                        end
                    end
                    
                    if moveon == 0 & it == storage_oscar.iterativemode.maxiterations
                        disp(['WARNING: Could not stabilize template after ' num2str(storage_oscar.iterativemode.maxiterations) ' iterations, giving up.']);
                        moveon = 1;
                    end
                    
                    if moveon == 0 %next iteration with same angular range

                        % crop average
                        average=average.Value;
                        average = average((size(average,1)./2+1)-(tmpheader.Header.Size(1)./2):(size(average,1)./2+1)+(tmpheader.Header.Size(1)./2-1), ...
                            (size(average,2)./2+1)-(tmpheader.Header.Size(2)./2):(size(average,2)./2+1)+(tmpheader.Header.Size(2)./2-1), ...
                            (size(average,3)./2+1)-(tmpheader.Header.Size(3)./2):(size(average,3)./2+1)+(tmpheader.Header.Size(3)./2-1));
                        average = tom_norm(average,'phase');
                        average = tom_emheader(average);
                        average.Header.Objectpixelsize = tempheader.Header.Objectpixelsize;
                        tom_emwrite([storage_oscar.iterativemode.datadir '/run' num2str(n) '_template_' num2str(it+1) '.vol'],average);
                        tmplatestring = [storage_oscar.iterativemode.datadir '/run' num2str(n) '_template_' num2str(it+1) '.vol'];

                    else %next refinement step with new angles
                        
                        if n < storage_oscar.iterativemode.nrrefinements
                        
                            % crop average
                            average=average.Value;
                            average = average((size(average,1)./2+1)-(tmpheader.Header.Size(1)./2):(size(average,1)./2+1)+(tmpheader.Header.Size(1)./2-1), ...
                                (size(average,2)./2+1)-(tmpheader.Header.Size(2)./2):(size(average,2)./2+1)+(tmpheader.Header.Size(2)./2-1), ...
                                (size(average,3)./2+1)-(tmpheader.Header.Size(3)./2):(size(average,3)./2+1)+(tmpheader.Header.Size(3)./2-1));
                            average = tom_norm(average,'phase');
                            average = tom_emheader(average);
                            average.Header.Objectpixelsize = tempheader.Header.Objectpixelsize;
                            tom_emwrite([storage_oscar.iterativemode.datadir '/run' num2str(n+1) '_template_1.vol'],average);
                            tmplatestring = [storage_oscar.iterativemode.datadir '/run' num2str(n+1) '_template_1.vol'];

                            %paste new volume
                            disp(strvcat([' '], ['Creating new input volume ...']));
                            Align = tom_av3_create_alignlist(filenames, Align, 0);
                            Align = tom_av3_paste4oscar(filenames, [storage_oscar.iterativemode.datadir '/run' num2str(n+1) '_volume.vol'], [storage_oscar.iterativemode.datadir '/paste_mask.vol'], Align, Align(1,1).Filter, Align(1,1).NormFlag, 2);
                            save([storage_oscar.iterativemode.datadir '/align.mat'], 'Align');

                            disp('template stabilized, moving to next refinement step ... ' );
                        
                        end
                        %create new angles
                        angle_incr = ceil(angle_incr./2);
                        angle_start = [-angle_incr(1).*2, -angle_incr(2).*2, -angle_incr(3).*2];
                        angle_end = [angle_incr(1).*2, angle_incr(2).*2, angle_incr(3).*2];

                        moveon = 0;
                        break;
                    end

                end

            end

            %refinement overview picture
            save([storage_oscar.iterativemode.datadir '/align.mat'], 'Align');
            disp('Generating refinement overview ...');
            pic=tom_show_refinement(avgrefinementlist);
            tom_emwrite([storage_oscar.iterativemode.datadir '/avgrefinement.em'],pic');
            
            avgvol = tom_emreadc([storage_oscar.iterativemode.datadir '/run' num2str(storage_oscar.iterativemode.nrrefinements) '_avg_' num2str(it) '.vol']);
            tom_volxyz(avgvol.Value);
            imtool(pic');

        
        %simple mode
        %%%%%%%%%%%%%%%%%%%%%%%
        else

            for n = 1:storage_oscar.iterativemode.nrrefinements

                if mod(angle_incr.*1000,10) ~= 0
                    disp('Angular Increment reached precision limit, terminating.');
                    break;
                end
                
                disp(strvcat([' '],[num2str(n) '. refinement step'],['********************************'],[' ']));
                disp(['Angles: ' num2str(angle_start(1)) ' ' num2str(angle_end(1)) ' ' num2str(angle_incr(1)) ' ' num2str(angle_start(2)) ' ' num2str(angle_end(2)) ' ' num2str(angle_incr(2)) ' ' num2str(angle_start(3)) ' ' num2str(angle_end(3)) ' ' num2str(angle_incr(3))]);

                %run oscar
                status = tom_av3_oscar_feed([storage_oscar.iterativemode.datadir '/run' num2str(n) '_volume.vol'],...
                    tmplatestring,...
                    [storage_oscar.iterativemode.datadir '/run' num2str(n) '_Out'],...
                    [storage_oscar.iterativemode.datadir '/psf.vol'],...
                    [storage_oscar.iterativemode.datadir '/mask.vol'],...
                    [storage_oscar.iterativemode.datadir '/run' num2str(n) '_oscar.log'],...
                    storage_oscar.oscar.fftsize, storage_oscar.oscar.pe, angle_start, angle_incr, angle_end, 0);
                if status == 0
                    disp('Oscar terminated with problems, exiting.');
                    return;
                end

                %extract
                disp(strvcat([' '],['Extracting data ...']));
                a = load([storage_oscar.iterativemode.datadir '/align.mat']);
                Align = tom_av3_extract_anglesshifts([storage_oscar.iterativemode.datadir '/run' num2str(n) '_Out'], a.Align, mask, 2);
                save([storage_oscar.iterativemode.datadir '/align.mat'], 'Align');

                %generate average
                disp(strvcat([' '],['Generating average ...']));
                al = tom_av3_align_sum(Align);
                average=tom_av3_average(al,'sum',0,0,2);
                average = tom_emheader(average);
                tempheader = tom_reademheader(al(1,1).Filename);
                average.Header.Objectpixelsize = tempheader.Header.Objectpixelsize;
                tom_emwrite([storage_oscar.iterativemode.datadir '/run' num2str(n) '_avg.vol'],average);
                avgrefinementlist{n} = [storage_oscar.iterativemode.datadir '/run' num2str(n) '_avg.vol'];

                %paste new volume
                if n < storage_oscar.iterativemode.nrrefinements
                    disp(strvcat([' '], ['Creating new input volume ...']));
                    Align = tom_av3_create_alignlist(filenames, Align, 0);
                    Align = tom_av3_paste4oscar(filenames, [storage_oscar.iterativemode.datadir '/run' num2str(n+1) '_volume.vol'], [storage_oscar.iterativemode.datadir '/paste_mask.vol'], Align, Align(1,1).Filter, Align(1,1).NormFlag, 2);
                    save([storage_oscar.iterativemode.datadir '/align.mat'], 'Align');
                end

                %create new angles
                angle_incr = ceil(angle_incr./2);
                angle_start = [-angle_incr(1).*2, -angle_incr(2).*2, -angle_incr(3).*2];
                angle_end = [angle_incr(1).*2, angle_incr(2).*2, angle_incr(3).*2];


            end

            %refinement overview picture
            %disp('Generating refinement overview ...');
            %pic=tom_show_refinement(avgrefinementlist);
            %tom_emwrite([storage_oscar.iterativemode.datadir '/avgrefinement.em'],pic');

            avgvol = tom_emreadc([storage_oscar.iterativemode.datadir '/run' num2str(storage_oscar.iterativemode.nrrefinements) '_avg.vol']);
            tom_volxyz(avgvol.Value);
            %imtool(pic');
            
        end
        
    
    %generate script
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    else

        filename = '/particle_auto_refinement.m';
        pathname = storage_oscar.iterativemode.datadir;

        %open script file for writing
        fid = fopen([pathname filename],'wt');
        fprintf(fid,'%%This is an automatically generated file from tom_av3_oscargui.\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n');
        fprintf(fid,'function particle_auto_refinement()\n\n');
        [r,s]=unix('whoami');
        fprintf(fid,'disp(''********************************************************************************************'');\n');
        fprintf(fid,'disp(''Automatic 3D particle alignment script, V1.1, generated at %s by %s'');\n', datestr(now), s(1:size(s,2)-1));
        fprintf(fid,'disp(''********************************************************************************************'');\n');
        disp(strvcat([''],['**************************'],['Generating script for iterative run, please wait ...'],['**************************'],['']));
        disp(strvcat(['Generating script for ' num2str(storage_oscar.iterativemode.nrrefinements) ' Iterations.'],[' ']));

        %fftsize and nrcpus
        fprintf(fid,'fftsize = %i;\n',storage_oscar.oscar.fftsize);
        fprintf(fid,'nrcpus = %i;\n\n',storage_oscar.oscar.pe);

        
        
        %template recursive mode
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if storage_oscar.iterativemode.templaterecursive == 1
        
            %copy files to datadir
            disp(['Copying data to datadir ' storage_oscar.iterativemode.datadir ' ...']);

            unix(['cp ' storage_oscar.template.filename ' ' storage_oscar.iterativemode.datadir '/run1_template_1.vol']);
            fprintf(fid,'tmpheader = tom_reademheader(''run1_template_1.vol'');\n\n');
            
            unix(['cp ' storage_oscar.volume.filename ' ' storage_oscar.iterativemode.datadir '/run1_volume.vol']);
            unix(['cp ' storage_oscar.psf.filename ' ' storage_oscar.iterativemode.datadir '/psf.vol']);
            unix(['cp ' storage_oscar.mask.filename ' ' storage_oscar.iterativemode.datadir '/mask.vol']);
            unix(['cp ' storage_oscar.iterativemode.extractionmask ' ' storage_oscar.iterativemode.datadir '/extract_mask.vol']);
            unix(['cp ' storage_oscar.iterativemode.pastemask ' ' storage_oscar.iterativemode.datadir '/paste_mask.vol']);
            unix(['cp ' storage_oscar.alignlist ' ' storage_oscar.iterativemode.datadir '/align.mat']);

            
            fprintf(fid,'%%Generic file definitions, used in all runs\n');
            fprintf(fid,'psffile = ''psf.vol'';\n');
            fprintf(fid,'maskfile = ''mask.vol'';\n');
            fprintf(fid,'extractionmaskfile = ''extraction_mask.vol'';\n');
            fprintf(fid,'pastemaskfile = ''paste_mask.vol'';\n');
            fprintf(fid,'alignmentfile = ''align.mat'';\n');
            
            %prepare subvolumes and adjust alignlist to the new filenames
            mkdir([storage_oscar.iterativemode.datadir '/particles']);
            disp(['Copying data to datadir ' storage_oscar.iterativemode.datadir '/particles ...']);
            fprintf(fid,'\n%%Filenames of subvolumes\n');
            a = load([storage_oscar.iterativemode.datadir '/align.mat']);
            Align = a.Align;
            fprintf(fid,'filenames = strvcat({ ');
            j = size(Align,1);
            for k = 1:size(Align,2)
                [pathstr, name, ext] = fileparts(Align(j,k).Filename);
                unix(['cp ' Align(j,k).Filename ' ' storage_oscar.iterativemode.datadir '/particles/particle_' num2str(k) ext]);
                Align(j,k).Tempfilename = Align(j,k).Filename;
                Align(j,k).Filename = ['particles/particle_' num2str(k) ext];
                fprintf(fid,'''%s'' ', Align(j,k).Filename);
                if mod(k,10) == 0
                    fprintf(fid,' ... \n');
                end
            end
            fprintf(fid,'});\n');
            save([storage_oscar.iterativemode.datadir '/align.mat'],'Align');
            
            %create angles
            angle_start = [storage_oscar.oscar.phi.start, storage_oscar.oscar.psi.start, storage_oscar.oscar.theta.start];
            angle_end = [storage_oscar.oscar.phi.stop, storage_oscar.oscar.psi.stop, storage_oscar.oscar.theta.stop];
            angle_incr = [storage_oscar.oscar.phi.step, storage_oscar.oscar.psi.step, storage_oscar.oscar.theta.step];

            %load mask and initialize refinement list
            fprintf(fid, '\n%%load mask and initialize refinement list\n');
            fprintf(fid, 'mask = tom_emreadc(''extract_mask.vol'');\nmask = mask.Value;\n');
            fprintf(fid, 'avgrefinementlist = {};\navgcounter = 1;\n');

            %build loop for refinement steps
            for n = 1:storage_oscar.iterativemode.nrrefinements
            
                if mod(angle_incr.*1000,10) ~= 0
                    disp('Angular Increment reached precision limit, terminating.');
                    break;
                end
                fprintf(fid,'\n\n%%refinement step %i\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n',n);
                fprintf(fid, 'for it=1:%i\n\n',storage_oscar.iterativemode.maxiterations);
                fprintf(fid, '   disp(strvcat(['' ''],[''===========================================''],[''%i. refinement step, '' num2str(it) ''. iteration''],[''===========================================''],['' '']));\n',n);
                
                %oscar run
                fprintf(fid,'   angle_start = [%0.2f %0.2f %0.2f];\n',angle_start(1),angle_start(2),angle_start(3));
                fprintf(fid,'   angle_incr = [%0.2f %0.2f %0.2f];\n',angle_incr(1),angle_incr(2),angle_incr(3));
                fprintf(fid,'   angle_end = [%0.2f %0.2f %0.2f];\n',angle_end(1),angle_end(2),angle_end(3));
                fprintf(fid,'   disp(''Angles: %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f'');\n\n', angle_start(1), angle_end(1), angle_incr(1), angle_start(2), angle_end(2), angle_incr(2), angle_start(3), angle_end(3), angle_incr(3));
                fprintf(fid,'   %%run oscar\n');
                fprintf(fid,'   status = tom_av3_oscar_feed(''run%i_volume.vol'', [''run%i_template_'' num2str(it) ''.vol''], ''run%i_Out'', ''psf.vol'', ''mask.vol'', ''run%i_oscar.log'', fftsize, nrcpus, angle_start, angle_incr, angle_end, 0);\n',n,n,n,n);
                fprintf(fid,'   if status == 0\n      disp(''Oscar terminated with problems, exiting.'');\n      return;\n   end\n\n');

                %extract
                fprintf(fid,'   %%extract values from oscar result\n');
                fprintf(fid,'   disp(strvcat(['' ''],[''Extracting data ...'']));\n');
                fprintf(fid,'   a = load(''align.mat'');\n');
                fprintf(fid,'   Align = tom_av3_extract_anglesshifts(''run%i_Out'', a.Align, mask, 2);\n',n);
                
                %generate average
                fprintf(fid,'\n   %%generate average\n');
                fprintf(fid,'   disp(strvcat(['' ''],[''Generating average ...'']));\n');
                fprintf(fid,'   al = tom_av3_align_sum(Align);\n');
                fprintf(fid,'   average = tom_av3_average(al,''sum'',0,0,2);\n');
                fprintf(fid,'   average = tom_emheader(average);\n');
                fprintf(fid,'   tempheader = tom_reademheader(al(1,1).Filename);\n');
                fprintf(fid,'   average.Header.Objectpixelsize = tmpheader.Header.Objectpixelsize;\n');
                fprintf(fid,'   tom_emwrite([''run%i_avg_'' num2str(it) ''.vol''], average);\n',n);
                fprintf(fid,'   avgrefinementlist{avgcounter} =  [''run%i_avg_'' num2str(it) ''.vol''];\n',n);
                fprintf(fid,'   avgcounter = avgcounter + 1;\n');

                %calculate derivative of templates from this run and the
                %last two runs, if it is smaller than the convergence
                %criteria, refine the angular range and start again.
                fprintf(fid,'   \n   %%calculate derivative of templates from this run and the\n   %%last two runs, if it is smaller than the convergence\n   %%criteria, refine the angular range and start again.\n');
                fprintf(fid,'   moveon = 0;\n');
                fprintf(fid,'   if it > 2\n');
                fprintf(fid,'      last_1=tom_emreadc([''run%i_avg_'' num2str(it-1) ''.vol'']);\n',n);
                fprintf(fid,'      last_2=tom_emreadc([''run%i_avg_'' num2str(it-2) ''.vol'']);\n',n);
                fprintf(fid,'      average_ccc(1)=tom_ccc(last_1.Value,last_2.Value,''norm'');\n');
                fprintf(fid,'      average_ccc(2)=tom_ccc(average.Value,last_1.Value,''norm'');\n');
                fprintf(fid,'      disp(strvcat(['' ''],[''absolute ccc values of last three averages: '', num2str(average_ccc(1)), '', '' num2str(average_ccc(2))]));\n');
                fprintf(fid,'      diff_average_ccc=diff(average_ccc);\n');
                fprintf(fid,'      disp(strvcat(['' ''],[''*** differential average ccc: '' num2str(diff_average_ccc), '' ***''],['' '']));\n');
                fprintf(fid,'      if abs(diff_average_ccc) < %0.4f | average_ccc(2) > %0.4f\n',storage_oscar.iterativemode.convcrit,1-storage_oscar.iterativemode.convcrit);
                fprintf(fid,'         moveon = 1;\n');
                fprintf(fid,'      end\n');
                fprintf(fid,'   end\n\n');
                fprintf(fid,'   if moveon == 0 & it == %i\n',storage_oscar.iterativemode.maxiterations);
                fprintf(fid,'      disp([''WARNING: Could not stabilize template after %i iterations, giving up and continuing.'']);\n',storage_oscar.iterativemode.maxiterations);
                fprintf(fid,'      moveon = 1;\n');
                fprintf(fid,'   end\n\n');
                
                fprintf(fid,'   if moveon == 0 %%next iteration with same angular range\n\n');
                fprintf(fid,'      %% crop average\n');
                fprintf(fid,'      average=average.Value;\n');
                fprintf(fid,'      average = average((size(average,1)./2+1)-(tmpheader.Header.Size(1)./2):(size(average,1)./2+1)+(tmpheader.Header.Size(1)./2-1), ...\n');
                fprintf(fid,'                        (size(average,2)./2+1)-(tmpheader.Header.Size(2)./2):(size(average,2)./2+1)+(tmpheader.Header.Size(2)./2-1), ...\n');
                fprintf(fid,'                        (size(average,3)./2+1)-(tmpheader.Header.Size(3)./2):(size(average,3)./2+1)+(tmpheader.Header.Size(3)./2-1));\n');
                fprintf(fid,'      average = tom_norm(average,''phase'');\n');
                fprintf(fid,'      average = tom_emheader(average);\n');
                fprintf(fid,'      average.Header.Objectpixelsize = tempheader.Header.Objectpixelsize;\n');
                fprintf(fid,'      tom_emwrite([''run%i_template_'' num2str(it+1) ''.vol''],average);\n\n',n);
                fprintf(fid,'   else %%next refinement step with new angles\n\n');
                
                if n < storage_oscar.iterativemode.nrrefinements
                   
                    fprintf(fid,'      %% crop average\n');
                    fprintf(fid,'      average=average.Value;\n');
                    fprintf(fid,'      average = average((size(average,1)./2+1)-(tmpheader.Header.Size(1)./2):(size(average,1)./2+1)+(tmpheader.Header.Size(1)./2-1), ...\n');
                    fprintf(fid,'                        (size(average,2)./2+1)-(tmpheader.Header.Size(2)./2):(size(average,2)./2+1)+(tmpheader.Header.Size(2)./2-1), ...\n');
                    fprintf(fid,'                        (size(average,3)./2+1)-(tmpheader.Header.Size(3)./2):(size(average,3)./2+1)+(tmpheader.Header.Size(3)./2-1));\n');
                    fprintf(fid,'      average = tom_norm(average,''phase'');\n');
                    fprintf(fid,'      average = tom_emheader(average);\n');
                    fprintf(fid,'      average.Header.Objectpixelsize = tempheader.Header.Objectpixelsize;\n');
                    fprintf(fid,'      tom_emwrite([''run%i_template_1.vol''], average);\n\n',n+1);

                    fprintf(fid,'      %%paste new volume\n');
                    fprintf(fid,'      disp(strvcat(['' ''], [''Creating new input volume ...'']));\n');
                    fprintf(fid,'      Align = tom_av3_create_alignlist(filenames, Align, 0);\n');
                    fprintf(fid,'      Align = tom_av3_paste4oscar(filenames, ''run%i_volume.vol'', ''paste_mask.vol'', Align, Align(1,1).Filter, Align(1,1).NormFlag, 2);\n',n+1);
                    fprintf(fid,'      save(''align.mat'', ''Align'');\n');
                    fprintf(fid,'      disp(''template stabilized, moving to next refinement step ... '');\n');
                
                end

                fprintf(fid,'      moveon = 0;\n');
                fprintf(fid,'      break;\n\n');
                fprintf(fid,'   end\n');
                
                fprintf(fid, 'end');

                angle_incr = ceil(angle_incr./2);
                angle_start = [-angle_incr(1).*2, -angle_incr(2).*2, -angle_incr(3).*2];
                angle_end = [angle_incr(1).*2, angle_incr(2).*2, angle_incr(3).*2];
                

            end

            %refinement overview picture
            fprintf(fid,'\n\nsave(''align.mat'', ''Align'');\n');
            %fprintf(fid,'\n\n%%refinement overview picture\n');
            %fprintf(fid,'disp(''Generating refinement overview ...'');\n');
            %fprintf(fid,'pic=tom_show_refinement(avgrefinementlist);\n');
            %fprintf(fid,'tom_emwrite(''avgrefinement.em'',pic'');\n');
            fprintf(fid,'\n\ndisp(''FINISHED.'');\n');
            fclose(fid);
            disp('Finished.');
            
            
            
        %simple mode
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else
            
            %copy files to datadir
            disp(['Copying data to datadir ' storage_oscar.iterativemode.datadir ' ...']);

            unix(['cp ' storage_oscar.template.filename ' ' storage_oscar.iterativemode.datadir '/template.vol']);
            templatestring = 'template.vol';
            fprintf(fid,'tmpheader = tom_reademheader(''template.vol'');\n\n');
            
            unix(['cp ' storage_oscar.volume.filename ' ' storage_oscar.iterativemode.datadir '/run1_volume.vol']);
            unix(['cp ' storage_oscar.psf.filename ' ' storage_oscar.iterativemode.datadir '/psf.vol']);
            unix(['cp ' storage_oscar.mask.filename ' ' storage_oscar.iterativemode.datadir '/mask.vol']);
            unix(['cp ' storage_oscar.iterativemode.extractionmask ' ' storage_oscar.iterativemode.datadir '/extract_mask.vol']);
            unix(['cp ' storage_oscar.iterativemode.pastemask ' ' storage_oscar.iterativemode.datadir '/paste_mask.vol']);
            unix(['cp ' storage_oscar.alignlist ' ' storage_oscar.iterativemode.datadir '/align.mat']);

            fprintf(fid,'%%Generic file definitions, used in all runs\n');
            fprintf(fid,'psffile = ''psf.vol'';\n');
            fprintf(fid,'maskfile = ''mask.vol'';\n');
            fprintf(fid,'extractionmaskfile = ''extraction_mask.vol'';\n');
            fprintf(fid,'pastemaskfile = ''paste_mask.vol'';\n');
            fprintf(fid,'alignmentfile = ''align.mat'';\n');
            fprintf(fid,'templatefile = ''template.vol'';\n');
            
            %prepare subvolumes and adjust alignlist to the new filenames
            mkdir([storage_oscar.iterativemode.datadir '/particles']);
            disp(['Copying data to datadir ' storage_oscar.iterativemode.datadir '/particles ...']);
            fprintf(fid,'\n%%Filenames of subvolumes\n');
            a = load([storage_oscar.iterativemode.datadir '/align.mat']);
            Align = a.Align;
            fprintf(fid,'filenames = strvcat({ ');
            j = size(Align,1);
            for k = 1:size(Align,2)
                [pathstr, name, ext] = fileparts(Align(j,k).Filename);
                unix(['cp ' Align(j,k).Filename ' ' storage_oscar.iterativemode.datadir '/particles/particle_' num2str(k) ext]);
                Align(j,k).Tempfilename = Align(j,k).Filename;
                Align(j,k).Filename = ['particles/particle_' num2str(k) ext];
                fprintf(fid,'''%s'' ', Align(j,k).Filename);
                if mod(k,10) == 0
                    fprintf(fid,' ... \n');
                end
            end
            fprintf(fid,'});\n');
            save([storage_oscar.iterativemode.datadir '/align.mat'],'Align');

            %create angles
            angle_start = [storage_oscar.oscar.phi.start, storage_oscar.oscar.psi.start, storage_oscar.oscar.theta.start];
            angle_end = [storage_oscar.oscar.phi.stop, storage_oscar.oscar.psi.stop, storage_oscar.oscar.theta.stop];
            angle_incr = [storage_oscar.oscar.phi.step, storage_oscar.oscar.psi.step, storage_oscar.oscar.theta.step];

            %load mask and initialize refinement list
            fprintf(fid, 'mask = tom_emreadc(''extract_mask.vol'');\nmask = mask.Value;\n');
            avgrefinementlist = '{';                

            %Loop over iterations
            for n = 1:storage_oscar.iterativemode.nrrefinements

                if mod(angle_incr.*1000,10) ~= 0
                    disp('Angular Increment reached precision limit, terminating.');
                    break;
                end;

                fprintf(fid,'\n\n%%refinement step %i\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n',n);
                fprintf(fid,'disp(strvcat(['' ''],['' ''],[''%i. refinement step''],[''********************************''],['' '']));\n',n);

                %oscar run
                fprintf(fid,'run%i_volumefile = ''%s'';\n',n,['run' num2str(n) '_volume.vol']);
                fprintf(fid,'angle_start = [%0.2f %0.2f %0.2f];\n',angle_start(1),angle_start(2),angle_start(3));
                fprintf(fid,'angle_incr = [%0.2f %0.2f %0.2f];\n',angle_incr(1),angle_incr(2),angle_incr(3));
                fprintf(fid,'angle_end = [%0.2f %0.2f %0.2f];\n',angle_end(1),angle_end(2),angle_end(3));
                fprintf(fid,'disp(''Angles: %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f'');\n\n', angle_start(1), angle_end(1), angle_incr(1), angle_start(2), angle_end(2), angle_incr(2), angle_start(3), angle_end(3), angle_incr(3));
                fprintf(fid,'%%run oscar\n');
                fprintf(fid,'status = tom_av3_oscar_feed(''run%i_volume.vol'', ''template.vol'', ''run%i_Out'', ''psf.vol'', ''mask.vol'', ''run%i_oscar.log'', fftsize, nrcpus, angle_start, angle_incr, angle_end, 0);\n',n,n,n);
                fprintf(fid,'if status == 0\n   disp(''Oscar terminated with problems, exiting.'');\n   return;\nend\n\n');

                %extract
                fprintf(fid,'%%extract values from oscar result\n');
                fprintf(fid,'disp(strvcat(['' ''],[''Extracting data ...'']));\n');
                fprintf(fid,'a = load(''align.mat'');\n');
                fprintf(fid,'Align = tom_av3_extract_anglesshifts(''run%i_Out'', a.Align, mask, 2);\n',n);

                %generate average
                fprintf(fid,'\n%%generate average\n');
                fprintf(fid,'disp(strvcat(['' ''],[''Generating average ...'']));\n');
                fprintf(fid,'al = tom_av3_align_sum(Align);\n');
                fprintf(fid,'average = tom_av3_average(al,''sum'',0,0,2);\n');
                fprintf(fid,'average = tom_emheader(average);\n');
                fprintf(fid,'tempheader = tom_reademheader(al(1,1).Filename);\n');
                fprintf(fid,'average.Header.Objectpixelsize = tmpheader.Header.Objectpixelsize;\n');
                fprintf(fid,'tom_emwrite(''run%i_avg.vol'',average);\n',n);
                avgrefinementlist = [avgrefinementlist ' ''run' num2str(n) '_avg.vol'' '];

                %paste new volume
                if n < storage_oscar.iterativemode.nrrefinements
                    fprintf(fid,'\n%%paste new volume\n');
                    fprintf(fid,'disp(strvcat(['' ''], [''Creating new input volume ...'']));\n');
                    fprintf(fid,'Align = tom_av3_create_alignlist(filenames, Align, 0);\n');
                    fprintf(fid,'Align = tom_av3_paste4oscar(filenames, ''run%i_volume.vol'', ''paste_mask.vol'', Align, [%i %i], %i, 2);\n', n+1, Align(1,1).Filter(1), Align(1,1).Filter(2), Align(1,1).NormFlag);
                    fprintf(fid,'save(''align.mat'', ''Align'');\n');

                    %create new angles
                    angle_incr = ceil(angle_incr./2);
                    angle_start = [-angle_incr(1).*2, -angle_incr(2).*2, -angle_incr(3).*2];
                    angle_end = [angle_incr(1).*2, angle_incr(2).*2, angle_incr(3).*2];

                end

            end

            %refinement overview picture
            fprintf(fid,'\nsave(''align.mat'', ''Align'');\n');
            avgrefinementlist = [avgrefinementlist '}'];
            %fprintf(fid,'\n%%refinement overview picture\n');
            %fprintf(fid,'disp(''Generating refinement overview ...'');\n');
            %fprintf(fid,'pic=tom_show_refinement(%s);\n',avgrefinementlist);
            %fprintf(fid,'tom_emwrite(''avgrefinement.em'',pic'');\n');
            fprintf(fid,'\n\ndisp(''FINISHED.'');\n');
            fclose(fid);
            disp('Finished.');
        end
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  check how many lamnodes are available                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function numnodes = check_lamnodes()

global storage_oscar;

[status,result] = unix('lamnodes');

numnodes = 0;

if strcmp(result(1:2),'n0') == 1
    match = regexp(result,':[0-9]:','match');
    
    for i = 1:size(match,2)
        string = match(i);
        string = string{1}(2:size(string{1},2)-1);
        numnodes = numnodes + str2num(string);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Proceed to next gui tom_av3_extractanglesshiftsgui                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_extract_Callback(hObject, eventdata, handles)

global storage_oscar;

tom_av3_extractanglesshiftsgui(storage_oscar.output.filename,storage_oscar.alignlist,storage_oscar.mask.filename);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Checkbox generate oscar script                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkbox_generatescript_Callback(hObject, eventdata, handles)

global storage_oscar;

storage_oscar.generatescript = get(hObject,'Value');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Unused Callbacks/Create functions                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function input_volume_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_template_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_mask_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_psf_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function phi_start_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function phi_stop_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function phi_step_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function psi_start_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function psi_stop_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function psi_step_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function theta_start_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function theta_stop_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function theta_step_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function fftsize_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function pe_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function output_files_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function output_log_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_logfile_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function output_runtimedata_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_extractionmask_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function iterative_number_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function iterative_maxiterations_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function iterative_convergence_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_pastemask_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end