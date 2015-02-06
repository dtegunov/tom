function varargout = tom_recparticles(varargin)
%TOM_RECPARTICLES reconstructs a particle with the hight resolution
%
%   varargout = tom_recparticles(varargin)
%
%  tom_recparticles is used to reconstruct a particle with the hight resolution.
%   It's very similar to tom_rec3d but it 's just for reconstructing a 3D
%   particles instead of the full volume. But for that, the function needs
%   to have a particle list file bilt just with tom_particles.
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
%   ... = tom_recparticles(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   TOM_PARTICLES ,TOM_RECONSTRUCTION3D, TOM_REC3D
%
%   created by WDN 05/20/2003
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

global mypathname myfilename myfirstnb mylastnb myext;


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_recparticles_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_recparticles_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


function tom_recparticles_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;

handles.WeightedProjection.Proj='';
handles.WeightedProjection.Path='';
handles.WeightedProjection.Name='';
handles.WeightedProjection.Ext='';
handles.WeightedProjection.First='';
handles.WeightedProjection.Last='';
handles.ParticlesList='';
handles.Sizex=128;
handles.Sizey=128;
handles.Sizez=128;
handles.Bin=0;
handles.Myfile='';
handles.Path='';
handles.File='';
handles.Ext='';
handles.ParaRecNum=2;
if nargin>3
    handles.ListMarkerPoint=varargin{1};  
    set(handles.particles_list,'String','-- To see the particles list --> Click Show');
    set(handles.list_button,'String','Show');
else
    handles.ListMarkerPoint='';
end;
% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = tom_recparticles_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --------------------------------------------------------------------
% -------   BUTTON BROWSE WEIGHTED PROJECTION   ----------------------
% --------------------------------------------------------------------
function weighted_button_Callback(hObject, eventdata, handles)
org_path=pwd;
if ~isempty(handles.WeightedProjection.Path)
    cd (handles.WeightedProjection.Path)
end
[filename, pathname]=uigetfile('*.*', 'Load an EM-file');
if isequal(filename,0) | isequal(pathname,0) 
    %error('No data loaded.'); 
    return;%nothing because Cancel is ckicked
else    
    a=findstr(filename,'.');
    if ~isempty(a)        
        b=size(a,2);
        v=[];w=[];
        for i=a(b):size(filename,2)
            v=strcat(v,filename(i));
        end
        handles.WeightedProjection.Ext=v;
        set(handles.weighted_ext,'String',v);
        for i=1:a(b)-1
            w=strcat(w,filename(i));
        end
        filename=w;
    else
        handles.WeightedProjection.Ext='';
    end
    a=findstr(filename,'_');
    b=size(a,2);
    v=[];w=[];
    for i=a(b)+1:size(filename,2)
        v=strcat(v,filename(i));
    end
    set(handles.weighted_last,'String',v);
    handles.WeightedProjection.Last=v;
    set(handles.weighted_first,'String','1');
    handles.WeightedProjection.First='1';    
    for i=1:a(b)
        w=strcat(w,filename(i));
    end
    set(handles.weighted_file,'String',w);
    handles.WeightedProjection.Name=w;
    set(handles.weighted_path,'String',pathname);
    handles.WeightedProjection.Path=pathname;
    handles.WeightedProjection.Proj=[pathname w];
end
cd (org_path);
guidata(hObject,handles);

% --------------------------------------------------------------------
% -------   EDITBOX PATH WEIGHTED PROJECTION   -----------------------
% --------------------------------------------------------------------
function weighted_path_Callback(hObject, eventdata, handles)
handles.WeightedProjection.Path=get(hObject,'String');
guidata(hObject,handles);
% --------------------------------------------------------------------
% -------   EDITBOX NAME WEIGHTED PROJECTION   -----------------------
% --------------------------------------------------------------------
function weighted_file_Callback(hObject, eventdata, handles)
nf=get(hObject,'String');
handles.WeightedProjection.Name=nf;
handles.WeightedProjection.Proj=[handles.WeightedProjection.Path nf];
guidata(hObject,handles);
% --------------------------------------------------------------------
% -------   EDITBOX EXTENSION WEIGHTED PROJECTION   ------------------
% --------------------------------------------------------------------
function weighted_ext_Callback(hObject, eventdata, handles)
e=get(hObject,'String');
if isempty(findstr(e,'.'))
    e=strcat('.',e);
end
handles.WeightedProjection.Ext=e;
guidata(hObject,handles);
% --------------------------------------------------------------------
% -------   EDITBOX FIRST WEIGHTED PROJECTION   ----------------------
% --------------------------------------------------------------------
function weighted_first_Callback(hObject, eventdata, handles)
handles.WeightedProjection.First=get(hObject,'String');
guidata(hObject,handles);
% --------------------------------------------------------------------
% -------   EDITBOX EDITBOX LAST WEIGHTED PROJECTION   ---------------
% --------------------------------------------------------------------
function weighted_last_Callback(hObject, eventdata, handles)
handles.WeightedProjection.Last=get(hObject,'String');
guidata(hObject,handles);
% --------------------------------------------------------------------
% -------   BUTTON BROWSE PARTICLES LIST   ---------------------------
% --------------------------------------------------------------------
function list_button_Callback(hObject, eventdata, handles)
a=get(hObject,'String'),
switch a
    case 'Browse'
        [filename, pathname]=uigetfile({'*.em';'*.*'}, 'Load an EM-file');
        if isequal(filename,0) | isequal(pathname,0)     
            return;%nothing because Cancel is ckicked
        else
            temp2=tom_emread([pathname filename]);
            handles.ListMarkerPoint=temp2.Value';    
            set(handles.particles_list,'String',[pathname filename]);
            handles.ParticlesList=[pathname filename];
            guidata(hObject,handles);
        end
    case 'Show'
        ListMarkerPoint=handles.ListMarkerPoint
end


% --------------------------------------------------------------------
% -------   EDIT TEXT SIZE X   ---------------------------------------
% --------------------------------------------------------------------
function sizex_Callback(hObject, eventdata, handles)
handles.Sizex=str2num(get(handles.sizex,'String'));
guidata(hObject,handles);
% --------------------------------------------------------------------
% -------   EDIT TEXT SIZE Y   ---------------------------------------
% --------------------------------------------------------------------
function sizey_Callback(hObject, eventdata, handles)
handles.Sizey=str2num(get(handles.sizey,'String'));
guidata(hObject,handles);
% --------------------------------------------------------------------
% -------   EDIT TEXT SIZE Z   ---------------------------------------
% --------------------------------------------------------------------
function sizez_Callback(hObject, eventdata, handles)
handles.Sizez=str2num(get(handles.sizez,'String'));
guidata(hObject,handles);
% --------------------------------------------------------------------
% -------   POPUP MENU BINNING   -------------------------------------
% --------------------------------------------------------------------
function binning_Callback(hObject, eventdata, handles)
a=get(handles.binning,'String');
b=get(handles.binning,'Value');
handles.Bin=str2num(a{b});
guidata(hObject,handles);
% --------------------------------------------------------------------
% -------   BUTTON BROWSE RECONSTRUCTION FILE   -------------------------------
% --------------------------------------------------------------------
function file_button_Callback(hObject, eventdata, handles)
org_path=pwd;
if ~isempty(handles.WeightedProjection.Path)
    cd (handles.WeightedProjection.Path)
end
[filename, pathname]=uiputfile('*.em', 'Save as');
if isequal(filename,0) | isequal(pathname,0) 
    %error('No data loaded.'); 
    return;%nothing because Cancel is ckicked
else
    set(handles.path,'String',pathname);
    if isempty(findstr(filename,'.em'))
        if ~isempty(findstr(filename,'.'))
            filename=strrep(filename,'.','_');
        end
        set(handles.name,'String',filename);
        extname='.em';
        set(handles.ext,'String',extname);
    else
        extname='.em';
        filename=strrep(filename,'.em','');
        set(handles.name,'String',filename);
        set(handles.ext,'String',extname);        
    end
    myfile=[pathname filename extname]
end
handles.Myfile=myfile;
handles.Path=pathname;
handles.File=filename;
handles.Ext=extname;
cd (org_path);
guidata(hObject,handles);
% --------------------------------------------------------------------
% -------   EDIT NUMBER OF RECONSTRUCTION   --------------------------
% --------------------------------------------------------------------
function parallel_Callback(hObject, eventdata, handles)
handles.ParaRecNum=str2num(get(handles.parallel,'String'));
guidata(hObject,handles);
% --------------------------------------------------------------------
% -------   BUTTON DO RECONSTRUCTION   -------------------------------
% --------------------------------------------------------------------
function Do_reconstruction_Callback(hObject, eventdata, handles)
if isempty(handles.WeightedProjection.Proj)
    message=['ERROR!!!   Enter the  name of a weighted projections.'];
    uiwait(msgbox(message,'Error','error'));
    return;
end
if isempty(handles.ListMarkerPoint)
    message=['ERROR!!!   There is no particles list.'];
    uiwait(msgbox(message,'Error','error'));
    return;
end
if isempty(handles.sizex)
    message=['ERROR!!!   Enter a number for the size X of the volume.'];
    uiwait(msgbox(message,'Error','error'));
    return;
end
if isempty(handles.sizey)
    message=['ERROR!!!   Enter a number for the size Y of the volume.'];
    uiwait(msgbox(message,'Error','error'));
    return;
end
if isempty(handles.sizez)
    message=['ERROR!!!   Enter a number for the size Z of the volume.'];
    uiwait(msgbox(message,'Error','error'));
    return;
end
if isempty(handles.Myfile)
    message=['ERROR!!!   Enter the file''s name of the reconstruction volume.'];
    uiwait(msgbox(message,'Error','error'));
    return;
end
%Volume=single(zeros(128,128,128));
%for lauf=1:87
%Projection=tom_emreadc(['TEMPV_BPP2048_' num2str(lauf)]);
%   disp(['Read in TEMPV_BPP2048_' num2str(lauf) ' tiltangle of ' num2str(Projection.Header.Tiltangle)]);
%tom_backproj3d(Volume,Projection.Value,0,Projection.Header.Tiltangle,[208 -904 -116]);
%end



laufx=0;
angles=0;
pathname=handles.WeightedProjection.Proj;
name=strrep(handles.WeightedProjection.Ext,'.','');


for lauf=1:str2num(handles.WeightedProjection.Last)
    if (tom_isemfile([pathname num2str(lauf) '.' name])==1)
        laufx=laufx+1;
        em=tom_reademheader([pathname num2str(lauf) '.' name]);
        angles(laufx)=em.Header.Tiltangle;
    end;

end;



Number_of_parallel_particles=(floor(size(handles.ListMarkerPoint,1)./handles.ParaRecNum)).*handles.ParaRecNum;

Number_of_rest_particles=size(handles.ListMarkerPoint,1)-Number_of_parallel_particles;
for i=1:handles.ParaRecNum:Number_of_parallel_particles
%for i=1:Number_of_parallel_particles

    for lauf=1:handles.ParaRecNum
        Volume(:,:,:,lauf)=single(zeros(handles.Sizex,handles.Sizey,handles.Sizez,'single'));
    end

    % read all the projections
    for lauf=str2num(handles.WeightedProjection.First):str2num(handles.WeightedProjection.Last)

        Projection=tom_emreadc([handles.WeightedProjection.Proj num2str(lauf) handles.WeightedProjection.Ext]);
        b=handles.Bin;

        if b>0
            Projection.Value=tom_bin(double(Projection.Value),b);
            Projection.Header.Size(1)=size(Projection.Value,1);
            Projection.Header.Size(2)=size(Projection.Value,2);
        end
        %disp(['Read file ' handles.WeightedProjection.Proj num2str(lauf) handles.WeightedProjection.Ext ' tiltangle of ' num2str(Projection.Header.Tiltangle)]);

        for lauf2=1:handles.ParaRecNum
            p1=handles.ListMarkerPoint(lauf2-1+i,10);
            p2=handles.ListMarkerPoint(lauf2-1+i,11);
            p3=handles.ListMarkerPoint(lauf2-1+i,12);
            if b>0
                for n=1:b
                    p1=p1./2;
                    p2=p2./2;
                    p3=p3./2;
                end
            end
            v=squeeze(Volume(:,:,:,lauf2));
            tom_backproj3d(single(v),single(Projection.Value),Projection.Header.Tiltaxis,Projection.Header.Tiltangle,[p1 p2 p3]);
            Volume(:,:,:,lauf2)=v;
        end
    end

    for lauf3=1:handles.ParaRecNum
        temp=[handles.Path handles.File '_' num2str(lauf3-1+i) handles.Ext];
        Vt=tom_emheader(squeeze(Volume(:,:,:,lauf3)));
        sz_tmp=Vt.Header.Size;
        Vt.Header=Projection.Header;
        Vt.Header.Size=sz_tmp;
        if (isunix==1)
            ind=(strfind(handles.WeightedProjection.Path,'/'));
        else
            ind=(strfind(handles.WeightedProjection.Path,'\'));
        end;
        ind1=ind(size(ind,2));
        ind2=ind(size(ind,2)-2);
        Vt.Header.Tiltaxis=min(angles);
        Vt.Header.Tiltangle=max(angles);
        part_path=handles.WeightedProjection.Path(ind2:ind1);
        try
            z_offset=num2str(handles.ListMarkerPoint(lauf3-1+i,14));
        catch
            z_offset=num2str(0);
        end;
        Vt.Header.Comment=[num2str(handles.ListMarkerPoint(lauf3-1+i,10)) ' ' num2str(handles.ListMarkerPoint(lauf3-1+i,11)) ' ' num2str(handles.ListMarkerPoint(lauf3-1+i,12)) ' ' z_offset ' ' num2str(b) ' ' part_path];
        Vt.Value=squeeze(Volume(:,:,:,lauf3));
        tom_emwrite(temp,Vt);
    end

end

    % reconstruction of the rest
        for i=Number_of_parallel_particles+1:Number_of_parallel_particles+Number_of_rest_particles

            Volume=single(zeros(handles.Sizex,handles.Sizey,handles.Sizez));
            % read all the projections
            for lauf=str2num(handles.WeightedProjection.First):str2num(handles.WeightedProjection.Last)
            
                Projection=tom_emreadc([handles.WeightedProjection.Proj num2str(lauf) handles.WeightedProjection.Ext]);
                b=handles.Bin;
        
                if b>0
                    Projection.Value=tom_bin(double(Projection.Value),b);
                    Projection.Header.Size(1)=size(Projection.Value,1);
                    Projection.Header.Size(2)=size(Projection.Value,2);
                end
                p1=handles.ListMarkerPoint(i,10);
                p2=handles.ListMarkerPoint(i,11);
                p3=handles.ListMarkerPoint(i,12);
                if b>0
                    for n=1:b
                        p1=p1./2;                    
                        p2=p2./2;                    
                        p3=p3./2;
                    end
                end                
                tom_backproj3d(single(Volume),single(Projection.Value),Projection.Header.Tiltaxis,Projection.Header.Tiltangle,[p1 p2 p3]);
            end
            temp=[handles.Path handles.File '_' num2str(i) handles.Ext];
            
            Vt=tom_emheader(squeeze(Volume));
            sz_tmp=Vt.Header.Size;
            Vt.Header=Projection.Header;
            Vt.Header.Size=sz_tmp;
            if (isunix==1)
               ind=(strfind(handles.WeightedProjection.Path,'/')); 
            else
               ind=(strfind(handles.WeightedProjection.Path,'\'));
            end;
            ind1=ind(size(ind,2));
            ind2=ind(size(ind,2)-2);
            Vt.Header.Tiltaxis=min(angles);
            Vt.Header.Tiltangle=max(angles);
            part_path=handles.WeightedProjection.Path(ind2:ind1);
            try
                z_offset=num2str(handles.ListMarkerPoint(lauf3-1+i,14));
            catch
                z_offset=num2str(0);
            end;

            Vt.Header.Comment=[num2str(handles.ListMarkerPoint(i,10)) ' ' num2str(handles.ListMarkerPoint(i,11)) ' ' num2str(handles.ListMarkerPoint(i,12)) ' ' z_offset ' ' num2str(b) ' ' part_path];
            Vt.Value=squeeze(Volume);
            tom_emwrite(temp,Vt);
        end
        %end

    
    
% --------------------------------------------------------------------
% -------   BUTTON EXIT   --------------------------------------------
% --------------------------------------------------------------------
function varargout = exit_Callback(h, eventdata, handles, varargin)
message=['Do you really want to quit 3D Reconstruction Particles?'];
Question=questdlg(message,'QUIT','Yes','No','Yes');
if strcmp(Question,'Yes')    
    if ~isempty(handles.rec_particles)
        delete(handles.rec_particles);
    end
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


