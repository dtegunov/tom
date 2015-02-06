function varargout = tom_processdirgui(varargin)
%TOM_PROCESSDIRGUI processes files from input directory to output directory
%
%   varargout = tom_processdirgui(varargin)
%
%   tom_processdirgui processes all files in the input directory and writes
%   them to the output directory. A variety of processing steps can be
%   selected by the user.
%   
%PARAMETERS
%
%  INPUT
%  
%  OUTPUT
%
%EXAMPLE
%   tom_amira_createisosurface(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   tom_processdir
%
%   created by AK 08/16/06
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
                   'gui_OpeningFcn', @tom_processdirgui_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_processdirgui_OutputFcn, ...
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
%%  opening function                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tom_processdirgui_OpeningFcn(hObject, eventdata, handles, varargin)


handles.actionstruct = [];

handles.output = hObject;



% Update handles structure
guidata(hObject, handles);

% UIWAIT makes tom_processdirgui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  output function                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = tom_processdirgui_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  browse input dir                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_processdir_inputdirbrowse_Callback(hObject, eventdata, handles)

pathname = uigetdir(pwd,'Input directory');
if ischar(pathname)
    set(handles.input_processdir_inputdir,'String',pathname);
end

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  browse output dir                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_processdir_browseoutdir_Callback(hObject, eventdata, handles)

pathname = uigetdir(pwd,'Output directory');
if ischar(pathname)
    set(handles.input_processdir_outputdir,'String',pathname);
end

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  browse button for field 5                                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_processdir_browsefield5_Callback(hObject, eventdata, handles)

contents = get(handles.dropdown_processdir_actions,'String');
selectedaction  = contents{get(handles.dropdown_processdir_actions,'Value')};

if isequal(selectedaction,'norm')
    [filename,pathname] = uigetfile('Mask file');
    if ischar(pathname)
        set(handles.edit_processdir_field5,'String',[pathname '/' filename]);
    end
end

guidata(hObject, handles);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  dropdown menu processing actions                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dropdown_processdir_actions_Callback(hObject, eventdata, handles)

handles = make_all_fields_invisible(handles);


contents = get(hObject,'String');
selectedaction  = contents{get(hObject,'Value')};

switch lower(selectedaction)
    case 'binning'
        set(handles.text_processdir_1,'Visible','on','String','number of times');
        set(handles.edit_processdir_field1,'Visible','on','String','1');

    case 'norm'
        set(handles.text_processdir_6,'Visible','on','String','norm method');
        set(handles.popupmenu_processdir_1,'Visible','on','String',{'phase','3std','2std','scaling factor','oscar'},'Value',1);
        set(handles.button_processdir_browsefield5,'Visible','on');
        set(handles.edit_processdir_field5,'Visible','on','String','');
        set(handles.text_processdir_5,'Visible','on','String','mask file (optional)');        
        
    case 'filter'
        set(handles.text_processdir_6,'Visible','on','String','method');
        set(handles.popupmenu_processdir_1,'Visible','on','String',{'fft','real space','anisotropic diffusion'},'Value',1);
        set(handles.edit_processdir_field1,'Visible','on','String','3');
        set(handles.text_processdir_1,'Visible','on','String','radius');
        set(handles.popupmenu_processdir_2,'Visible','on','String',{'circular kernel','quadratic kernel'},'Value',1);
        set(handles.text_processdir_7,'Visible','on','String','kernel');
        
    case 'taper'
        set(handles.text_processdir_1,'Visible','on','String','new size ([x y] or [x y z])');
        set(handles.edit_processdir_field1,'Visible','on','String','[1024,1024]');

    case 'mirror'
        set(handles.text_processdir_6,'Visible','on','String','axis');
        set(handles.popupmenu_processdir_1,'Visible','on','String',{'x','y','z'});
        
    case 'shift'
        set(handles.text_processdir_1,'Visible','on','String','shift vector ([x y] or [x y z])');
        set(handles.edit_processdir_field1,'Visible','on','String','[10,10]');
        
    case 'x-ray correction'
        
    case 'sort tilt series'
        set(handles.text_processdir_1,'Visible','on','String','new file name');
        set(handles.edit_processdir_field1,'Visible','on','String','sorted_');
                
    case 'convert file format'
        set(handles.text_processdir_6,'Visible','on','String','axis');
        set(handles.popupmenu_processdir_1,'Visible','on','String',{'mrc -> em','spi -> em','tif -> em'},'Value',1);

    case 'edit header'
        set(handles.popupmenu_processdir_1,'Visible','on','String',{'Voltage','Cs','Aperture','Magnification','Postmagnification','Exposuretime','Objectpixelsize','Pixelsize','CCDArea','Defocus','Astigmatism','Astigmatismangle','FocusIncrement','CountsPerElectron','Intensity','EnergySlitwidth','EnergyOffset','Tiltangle','Tiltaxis','Marker_X','Marker_Y','Microscope'},'Value',1);
        set(handles.edit_processdir_field1,'Visible','on','String','');
        set(handles.text_processdir_1,'Visible','on','String','value');

end

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  dropdown menu 1                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function popupmenu_processdir_1_Callback(hObject, eventdata, handles)

contents = get(hObject,'String');
selectedaction  = contents{get(hObject,'Value')};

switch lower(selectedaction)
    case 'scaling factor'
        set(handles.edit_processdir_field1,'Visible','on','String','1');
        set(handles.text_processdir_1,'Visible','on','String','scaling factor');
        set(handles.button_processdir_browsefield5,'Visible','on');
        set(handles.edit_processdir_field5,'Visible','on','String','');
        set(handles.text_processdir_5,'Visible','on','String','mask file (optional)');
        
    case 'phase'
        set(handles.edit_processdir_field1,'Visible','off');
        set(handles.text_processdir_1,'Visible','off');
        set(handles.button_processdir_browsefield5,'Visible','on');
        set(handles.edit_processdir_field5,'Visible','on','String','');
        set(handles.text_processdir_5,'Visible','on','String','mask file (optional)');
        
    case '3std'
        set(handles.edit_processdir_field1,'Visible','off');
        set(handles.text_processdir_1,'Visible','off');
        set(handles.button_processdir_browsefield5,'Visible','on');
        set(handles.edit_processdir_field5,'Visible','on','String','');
        set(handles.text_processdir_5,'Visible','on','String','mask file (optional)');        
        
    case '2std'
        set(handles.edit_processdir_field1,'Visible','off');
        set(handles.text_processdir_1,'Visible','off');
        set(handles.button_processdir_browsefield5,'Visible','on');
        set(handles.edit_processdir_field5,'Visible','on','String','');
        set(handles.text_processdir_5,'Visible','on','String','mask file (optional)');
        
    case 'oscar'
        set(handles.edit_processdir_field1,'Visible','off');
        set(handles.text_processdir_1,'Visible','off');
        set(handles.edit_processdir_field5,'Visible','on');
        set(handles.text_processdir_5,'Visible','on','String','mask file (optional)');
        set(handles.button_processdir_browsefield5,'Visible','on');
        
    case 'fft'
        set(handles.edit_processdir_field1,'Visible','on','String','3');
        set(handles.text_processdir_1,'Visible','on','String','radius');
        set(handles.popupmenu_processdir_2,'Visible','on','String',{'circular kernel','quadratic kernel'},'Value',1);
        set(handles.text_processdir_7,'Visible','on','String','kernel');
        set(handles.text_processdir_2,'Visible','off');
        set(handles.text_processdir_3,'Visible','off');
        set(handles.text_processdir_4,'Visible','off');
        set(handles.text_processdir_5,'Visible','off');
        set(handles.edit_processdir_field2,'Visible','off');
        set(handles.edit_processdir_field3,'Visible','off');
        set(handles.edit_processdir_field4,'Visible','off');
        set(handles.edit_processdir_field5,'Visible','off');

        
    case 'real space'
        set(handles.edit_processdir_field1,'Visible','on','String','3');
        set(handles.text_processdir_1,'Visible','on','String','radius');
        set(handles.popupmenu_processdir_2,'Visible','on','String',{'circular kernel','quadratic kernel'},'Value',1);
        set(handles.text_processdir_7,'Visible','on','String','kernel');
        set(handles.text_processdir_2,'Visible','off');
        set(handles.text_processdir_3,'Visible','off');
        set(handles.text_processdir_4,'Visible','off');
        set(handles.text_processdir_5,'Visible','off');
        set(handles.edit_processdir_field2,'Visible','off');
        set(handles.edit_processdir_field3,'Visible','off');
        set(handles.edit_processdir_field4,'Visible','off');
        set(handles.edit_processdir_field5,'Visible','off');
        
    case 'anisotropic diffusion'
        set(handles.popupmenu_processdir_2,'Visible','on','String',{'linearPsi','lorentzianPsi','tukeyPsi'},'Value',1);
        set(handles.text_processdir_7,'Visible','on','String','function');
        set(handles.edit_processdir_field1,'Visible','on','String','1');
        set(handles.text_processdir_1,'Visible','on','String','sigma');
        set(handles.edit_processdir_field2,'Visible','on','String','10');
        set(handles.text_processdir_2,'Visible','on','String','iterations');
        set(handles.edit_processdir_field3,'Visible','on','String','0.25');
        set(handles.text_processdir_3,'Visible','on','String','lambda'); 

end

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  add processing step to list                                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_processdir_addaction_Callback(hObject, eventdata, handles)

contents = get(handles.dropdown_processdir_actions,'String');
selectedaction  = contents{get(handles.dropdown_processdir_actions,'Value')};

switch lower(selectedaction)
    case '--- Select an action ---'
        errordlg('Select an action first!');
        return;
    
    case 'binning'
        timesbinning = round(str2double(get(handles.edit_processdir_field1,'String')));
        if ~isnumeric(timesbinning) || isempty(timesbinning)
            errordlg('number of times binning must be numeric!');
            return;
        end
        actioncmd = ['''bin'',' num2str(timesbinning)];
        actiondescription = ['bin ' num2str(timesbinning) ' times.'];
        
    case 'norm'
        contents = get(handles.popupmenu_processdir_1,'String');
        normmethod = contents{get(handles.popupmenu_processdir_1,'Value')};
        maskfile = get(handles.edit_processdir_field5,'String');

        if ~isempty(maskfile)
            try; 
                maskfile = tom_emreadc(maskfile);
                maskfile = maskfile.Value;
            catch
                errordlg('Mask file could not be opened.');
                return;
            end
        end
        
        if isequal(normmethod,'scaling factor')
            scalingfactor = str2double(get(handles.edit_processdir_field1,'String'));
            if ~isnumeric(scalingfactor) || isempty(scalingfactor)
                errordlg('scaling factor must be numeric!');
                return;
            end
            actioncmd = ['''norm'',' num2str(scalingfactor)];
            actiondescription = ['norm with scaling factor of ' num2str(scalingfactor)];

        else
            actioncmd = ['''norm'',''' normmethod ''''];
            actiondescription = [normmethod ' norm.'];
        end
        
    case 'filter'
        contents = get(handles.popupmenu_processdir_1,'String');
        filtermethod = contents{get(handles.popupmenu_processdir_1,'Value')};
        if isequal(filtermethod,'fft') || isequal(filtermethod,'real space')
            radius = round(str2double(get(handles.edit_processdir_field1,'String')));
            contents = get(handles.popupmenu_processdir_2,'String');
            kerneltype = contents{get(handles.popupmenu_processdir_2,'Value')};
            if isequal(kerneltype,'circular kernel')
                kerneltype2 = 'circ';
            else
                kerneltype2 = 'quadr';
            end
            if isequal(filtermethod,'fft')
                space = 'fourier';
            else
                space = 'real';
            end
            if ~isnumeric(radius) || isempty(radius)
                errordlg('radius must be numeric!');
                return;
            end
            actioncmd = ['''filter'',' num2str(radius) ',''' kerneltype2 ''',''' space ''''];
            actiondescription = ['filter in ' space ' space with a ' kerneltype ' of radius ' num2str(radius) ' pixels.'];
            
        else
            sigma = str2double(get(handles.edit_processdir_field1,'String'));
            iterations = round(str2double(get(handles.edit_processdir_field2,'String')));
            lambda = str2double(get(handles.edit_processdir_field3,'String'));
            if ~isnumeric(sigma) || isempty(sigma)
                errordlg('sigma must be numeric!');
                return;
            end
            if ~isnumeric(iterations) || isempty(iterations)
                errordlg('number of iterations must be numeric!');
                return;
            end
            if ~isnumeric(lambda) || isempty(lambda)
                errordlg('lambda must be numeric!');
                return;
            end
            contents = get(handles.popupmenu_processdir_2,'String');
            aniso_function = ['aniso_' contents{get(handles.popupmenu_processdir_2,'Value')}];
            actioncmd = ['''filter'',''' aniso_function ''',' num2str(sigma) ',' num2str(iterations) ',' num2str(lambda)];
            actiondescription = ['filter using anisotropic diffusion function ' aniso_function ' using sigma = ' num2str(sigma) ' for ' num2str(iterations) ' iterations, lamba = ' num2str(lambda)];
        end
        
        
    case 'taper'
        newsize = round(str2num(get(handles.edit_processdir_field1,'String')));
        if ~isnumeric(newsize) || isempty(newsize) || (length(newsize) ~=2 && length(newsize) ~=3)
            errordlg('number of times binning must be a 2 or 3 element vector!');
            return;
        end
        actioncmd = ['''taper'',' num2str(newsize)];
        actiondescription = ['taper to new size of ' num2str(newsize) ' pixels.'];
        
    case 'mirror'
        contents = get(handles.popupmenu_processdir_1,'String');
        axis = contents{get(handles.popupmenu_processdir_1,'Value')};
        actioncmd = ['''mirror'',''' axis ''''];
        actiondescription = ['mirror around ' axis ' axis.'];

    case 'shift'
        shiftvec = round(str2num(get(handles.edit_processdir_field1,'String')));
        if ~isnumeric(shiftvec) || isempty(shiftvec) || (length(shiftvec) ~=2 && length(shiftvec) ~=3)
            errordlg('number of times binning must be a 2 or 3 element vector!');
            return;
        end
        actioncmd = ['''shift'',' num2str(shiftvec)];
        actiondescription = ['shift by ' num2str(shiftvec) ' pixels.'];
        
    case 'x-ray correction'
        actioncmd = '''xraycorrect''';
        actiondescription = 'correct for x rays';
        
    case 'convert file format'
        contents = get(handles.popupmenu_processdir_1,'String');
        filetype = contents{get(handles.popupmenu_processdir_1,'Value')};
        if isequal(filetype,'mrc -> em')
            filetype = 'mrc';
        elseif isequal(filetype,'spi -> em');
            filetype = 'spi';
        elseif isequal(filetype,'tif -> em');
            filetype = 'tif';
        end
        actioncmd = ['convert_' filetype];
        actiondescription = ['convert files from ' filetype ' to em'];
        
    case 'sort tilt series'
        newfilename = get(handles.edit_processdir_field1,'String');
        if isempty(newfilename)
            errordlg('No filename given!');
            return;
        end

        actioncmd = ['sort ' newfilename];
        actiondescription = ['sort tilt series'];
        
    case 'edit header'
        contents = get(handles.popupmenu_processdir_1,'String');
        headerfieldname = contents{get(handles.popupmenu_processdir_1,'Value')};
        val = get(handles.edit_processdir_field1,'String');
        if ~isequal(headerfieldname,'Microscope')
            val = str2num(val);
        else
            val = ['''' val ''''];
        end
        if isempty(val)
            errordlg('Enter a value!');
        end
        actioncmd = ['''editheader'',''' headerfieldname ''',' num2str(val)];
        actiondescription = ['edit file header, set  ' headerfieldname ' to ' num2str(val)];
        
end

pos = length(handles.actionstruct) + 1;
handles.actionstruct(pos).cmds = actioncmd;
handles.actionstruct(pos).descriptions = actiondescription;
if exist('maskfile','var') && ~isempty(maskfile)
    handles.actionstruct(pos).maskfiles = maskfile;
else
    handles.actionstruct(pos).maskfiles = [];
end

%automatically select x ray correction if tiltseries is selected
if isequal(lower(selectedaction), 'sort tilt series')
    pos = length(handles.actionstruct) + 1;
    handles.actionstruct(pos).cmds = '''xraycorrect''';
    handles.actionstruct(pos).descriptions = 'correct for x rays';
end

handles = display_actions(handles);

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  remove action from list                                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_processdir_removeaction_Callback(hObject, eventdata, handles)

val = get(handles.listbox_processdir_processingsteps,'Value');

lauf = 1;
for i=1:length(handles.actionstruct)
    if i ~= val
        tmpstruct(lauf) = handles.actionstruct(i);
        lauf = lauf + 1;
    end
end

if exist('tmpstruct','var')
    handles.actionstruct = tmpstruct;
else
    handles.actionstruct = [];
end

set(handles.listbox_processdir_processingsteps,'Value',1);

handles = display_actions(handles);

guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  start processing                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function button_processdir_go_Callback(hObject, eventdata, handles)

indir = get(handles.input_processdir_inputdir,'String');
outdir = get(handles.input_processdir_outputdir,'String');

if isempty(indir)
    errordlg('No input directory specified.');
    return;
end
if isempty(outdir)
    errordlg('No output directory specified.');
    return;
end

filetype = 'em';
outstring = '';
sortfiles = 0;
for i=1:length(handles.actionstruct)
    if isequal(handles.actionstruct(i).cmds,'convert_mrc') 
        filetype = 'mrc';
    elseif isequal(handles.actionstruct(i).cmds,'convert_spi') 
        filetype = 'spi';
    elseif isequal(handles.actionstruct(i).cmds,'convert_tif') 
        filetype = 'tif';
    elseif findstr(handles.actionstruct(i).cmds,'sort ')
        sortfiles = 1;
        newname = handles.actionstruct(i).cmds(6:end);
    elseif ~isempty(findstr(handles.actionstruct(i).cmds,'norm')) & ~isempty(handles.actionstruct(i).maskfiles)
        outstring = [outstring ',' handles.actionstruct(i).cmds ', handles.actionstruct(' num2str(i) ').maskfiles'];
    else
        outstring = [outstring ',' handles.actionstruct(i).cmds];
    end
end

if sortfiles == 0
    if ~isempty(outstring)
        eval(['tom_processdir(''' indir ''',''' outdir ''',''' filetype ''',' outstring(2:end) ');']);
    else
        errordlg('Nothing to do!');
        return;
    end
    
else
    handles = sortseries(handles, indir, outdir, newname);
    if ~isempty(outstring)
        eval(['tom_processdir(''' outdir ''',''' outdir ''',''' filetype ''',' outstring(2:end) ');']);
    end
end

try
    h = fopen([indir '/protocol.txt'],'w');
    fprintf(h,['tom_processdirgui log file, generated at ' datestr(now) '\r\n\r\n']);
    fprintf(h,['input directory: ' indir '\r\n']);
    fprintf(h,['output directory: ' indir '\r\n']);
    
    for i=1:length(handles.actionstruct)
        fprintf(h,[handles.actionstruct(i).descriptions '\r\n']);
    end
    fclose(h);
catch
    disp('Warning: could not write protocol.txt in input dir.');
end
try
    h = fopen([outdir '/protocol.txt'],'w');
    fprintf(h,['tom_processdirgui log file, generated at ' datestr(now) '\r\n\r\n']);
    fprintf(h,['input directory: ' indir '\r\n']);
    fprintf(h,['output directory: ' indir '\r\n']);
    
    for i=1:length(handles.actionstruct)
        fprintf(h,[handles.actionstruct(i).descriptions '\r\n']);
    end
    fclose(h);
catch
    disp('Warning: could not write protocol.txt in output dir.');
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
%%  Make all fields invisible                                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = make_all_fields_invisible(handles)

set(handles.edit_processdir_field1,'Visible','off');
set(handles.edit_processdir_field2,'Visible','off');
set(handles.edit_processdir_field3,'Visible','off');
set(handles.edit_processdir_field4,'Visible','off');
set(handles.edit_processdir_field5,'Visible','off');
set(handles.popupmenu_processdir_1,'Visible','off');
set(handles.popupmenu_processdir_2,'Visible','off');

set(handles.text_processdir_1,'Visible','off');
set(handles.text_processdir_2,'Visible','off');
set(handles.text_processdir_3,'Visible','off');
set(handles.text_processdir_4,'Visible','off');
set(handles.text_processdir_5,'Visible','off');
set(handles.text_processdir_6,'Visible','off');
set(handles.text_processdir_7,'Visible','off');

set(handles.button_processdir_browsefield5,'Visible','off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display actions in listbox                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = display_actions(handles)

string = '';
for i=1:length(handles.actionstruct)
    string = strvcat(string, handles.actionstruct(i).descriptions);
end

set(handles.listbox_processdir_processingsteps,'String',string);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  sort tilt series                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = sortseries(handles, indir, outdir, newname)

h = waitbar(0,'Sorting tilt series...');
disp('Sorting tilt series...');
if ~isequal(newname(end),'_')
    newname = [newname '_'];
end

if exist(outdir,'dir') ~= 7
    disp('Output directory does not exist, creating new directory...');
    try
        mkdir(outdir);
    catch
        error('Could not create new directory!');
    end
end

d=dir(indir);
angles=zeros((size(d,1)-2),1);
laufx=0;
for lauf=3:(size(d,1))
    if (tom_isemfile([indir '/' d(lauf).name])==1)
        laufx=laufx+1;
        em=tom_reademheader([indir '/' d(lauf).name]);
        angles(laufx)=em.Header.Tiltangle;
    else    
        laufx=laufx+1;
        angles(laufx)=-9999; %just a dummy Value to keep the index in the right order :-)
    end;    

end
[y,Index]=sort(angles);

lauf2=1;
for lauf=1:(size(angles,1))
    if (tom_isemfile([indir '/' d(Index(lauf)+2).name]) == 1)
        [y, n, newext] = fileparts(d(Index(lauf)+2).name);
        [s,mess]=copyfile([indir '/' d(Index(lauf)+2).name],[outdir '/'  newname num2str(lauf2) newext],'f');
        if s == 0
            close(h);
            error(['Copying failed: ' mess]);
            
        end
        lauf2=lauf2+1;    
    end 
end
close(h);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Unused Callbacks/Create functions                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function listbox_processdir_processingsteps_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_processdir_outputdir_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_processdir_inputdir_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function input_processdir_inputdir_Callback(hObject, eventdata, handles)
function input_processdir_outputdir_Callback(hObject, eventdata, handles)
function dropdown_processdir_actions_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function listbox_processdir_processingsteps_Callback(hObject, eventdata, handles)


function edit_processdir_field1_Callback(hObject, eventdata, handles)
function edit_processdir_field1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_processdir_field2_Callback(hObject, eventdata, handles)
function edit_processdir_field2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_processdir_field3_Callback(hObject, eventdata, handles)
function edit_processdir_field3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_processdir_field4_Callback(hObject, eventdata, handles)
function edit_processdir_field4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_processdir_field5_Callback(hObject, eventdata, handles)
function edit_processdir_field5_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupmenu_processdir_1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function popupmenu_processdir_2_Callback(hObject, eventdata, handles)
function popupmenu_processdir_2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
