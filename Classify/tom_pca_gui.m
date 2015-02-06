function varargout = tom_pca_gui(varargin)
%TOM_PCA_GUI creates ...
%
%   varargout = tom_pca_gui(varargin)
%
% TOM_PCA_GUI M-file for tom_pca_gui.fig
%      TOM_PCA_GUI, by itself, creates a new TOM_PCA_GUI or raises the existing
%      singleton*.
%
%      H = TOM_PCA_GUI returns the handle to a new TOM_PCA_GUI or the handle to
%      the existing singleton*.
%
%      TOM_PCA_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TOM_PCA_GUI.M with the given input arguments.
%
%      TOM_PCA_GUI('Property','Value',...) creates a new TOM_PCA_GUI or
%      raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before tom_pca_gui_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to tom_pca_gui_OpeningFcn via varargin.
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
%   ... = tom_pca_gui(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   GUIDE, GUIDATA, GUIHANDLES
%
%   created by ... (author date)
%   updated by GUIDE v2.5 02-Jun-2006 14:21:10
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

error(nargchk(0, 1, nargin, 'struct'))


gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_pca_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_pca_gui_OutputFcn, ...
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


% --- Executes just before tom_pca_gui is made visible.
function tom_pca_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tom_pca_gui (see VARARGIN)

% Choose default command line output for tom_pca_gui
handles.output = hObject;
set(handles.basename,'Visible','off');
set(handles.extension,'Visible','off');
set(handles.last_part,'Visible','off');
set(handles.label_basename,'Visible','off');
set(handles.label_extension,'Visible','off');
set(handles.label_last_particle,'Visible','off');



% Update handles structure
guidata(hObject, handles);

% UIWAIT makes tom_pca_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = tom_pca_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;






% --- Executes on button press in browse_input.
function browse_input_Callback(hObject, eventdata, handles)
% hObject    handle to browse_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file path]=uigetfile('*.*');


set(handles.path_input,'String',[path file]);

st=get_gui_values(handles);

[path,basename,ext] = fileparts(st.in_path);
tmp_basename=basename;
tmp_extension=ext;
[tmp_basename rest]=strtok(basename,'_');

tmp_last_num=str2num(strrep(rest,'_',''));

set(handles.basename,'String',[tmp_basename '_']);
set(handles.extension,'String',tmp_extension);
set(handles.last_part,'String',tmp_last_num);









% --- Executes during object creation, after setting all properties.
function eigenvects_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eigenvects (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in disp_data.
function disp_data_Callback(hObject, eventdata, handles)
% hObject    handle to disp_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

st=get_gui_values(handles);
axes(handles.img);


% index1=size(handles.coefs,1)-(st.eigvects(1)-1);
% index2=size(handles.coefs,1)-(st.eigvects(2)-1);

index1=st.eigvects(1);
index2=st.eigvects(2);


if (size(st.eigvects,1)>2)
    %index3=size(handles.coefs,1)-(st.eigvects(3)-1);
    index3=st.eigvects(3);
    plot3(handles.scores(index1,:),handles.scores(index2,:),handles.scores(index3,:),'.');
    figure; plot3(handles.scores(index1,:),handles.scores(index2,:),handles.scores(index3,:),'.');
    set(gcf,'Name','pca_3d');
else
    plot(handles.scores(index1,:),handles.scores(index2,:),'r+');
end;



disp('mumu');

% --- Executes on button press in show_images.
function show_images_Callback(hObject, eventdata, handles)
% hObject    handle to show_images (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
st=get_gui_values(handles);

axes(handles.img);
[X,Y]=ginput(1);

v_tmp=[X Y];

% for i=1:size(st.eigvects,1)
%     index=size(handles.coefs,1)-(st.eigvects(i)-1);
%     tmp(:,i)=handles.coefs(:,index);
% end;

 for i=1:size(st.eigvects,1)
     index=st.eigvects(i);
     tmp(:,i)=handles.coefs(:,index);
 end;


img_tmp=tmp*v_tmp';



if (st.flag_3d==1)
    h=tom_reademheader(st.in_path);
    sz=h.Header.Size;
    sz=round(sz./(2^st.input_binning));   
    %figure; tom_dspcub(tom_rm_mean(reshape(img_tmp,sz(1),sz(2),sz(3)),handles.mean));
    
    figure; tom_dspcub(reshape(tom_rm_mean(img_tmp,handles.mean),sz(1),sz(2),sz(3)));
    
    %figure; tom_dspcub(reshape(img_tmp,sz(1),sz(2),sz(3)));
    
else
    h=tom_reademheader(st.in_path);
    sz=round(h.Header.Size');
    sz=round(sz./(2^st.input_binning));
    
    figure; tom_imagesc(reshape(tom_rm_mean(img_tmp,handles.mean),sz(1),sz(2) ));
end;



function st=get_gui_values(handles)


st.in_path=get(handles.path_input,'String');

tmp=get(handles.eigenvects,'String');
st.eigvects=sscanf(strrep(tmp,';',' '),'%d %d');
st.input_binning=str2num(get(handles.input_binning,'String'));


tmp=get(handles.edit_show_eig_gal,'String');
st.eigvects_gal=sscanf(strrep(tmp,'-',' '),'%d %d');

st.num_of_classes=str2num(get(handles.cluster_edit,'String'));
st.class=str2num(get(handles.edit_show_class,'String'));

st.basename=get(handles.basename,'String');
st.extension=get(handles.extension,'String');
st.last_particle=str2num(get(handles.last_part,'String'));
st.flag_3d=get(handles.check_3d,'Value');

st.num_of_eigs=str2num(get(handles.num_of_eigs,'String'));


tmp=get(handles.opt_num_classes,'String');
st.opt_num_of_classes=sscanf(strrep(tmp,'-',' '),'%d %d');

tmp=get(handles.opt_num_eig,'String');
st.opt_num_eig=sscanf(strrep(tmp,'-',' '),'%d %d');



% --- Executes on button press in show_orig_image.
function show_orig_image_Callback(hObject, eventdata, handles)
% hObject    handle to show_orig_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



st=get_gui_values(handles);
axes(handles.img);
[X,Y]=ginput(1);

v_tmp=[X Y];

for i=1:size(st.eigvects,1)
    %index=size(handles.scores,1)-(st.eigvects(i)-1);
    index=st.eigvects(i);
    tmp(:,i)=handles.scores(index,:);
end;

[pointidx, pointcoords,d] = tom_nearestpoint(v_tmp,tmp);
hold on; plot(pointcoords(1),pointcoords(2),'go'); hold off;


if (st.flag_3d==1)
    filenames=tom_path2cell(fileparts(st.in_path),st.basename,st.extension,[1 st.last_particle]);
    part=tom_emreadc(filenames{pointidx});
    figure; tom_dspcub(part.Value);
else
    h=tom_reademheader(st.in_path);
    part=tom_emreadc([st.in_path],'subregion',[1 1 pointidx],[h.Header.Size(1)-1 h.Header.Size(1)-1 0]);
    figure; tom_imagesc(part);
end;



% --- Executes on button press in show_eigvect_gallery.
function show_eigvect_gallery_Callback(hObject, eventdata, handles)
% hObject    handle to show_eigvect_gallery (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

st=get_gui_values(handles);


if (st.flag_3d==1)
   
    h=tom_reademheader(st.in_path);
    sz=h.Header.Size;
    sz=round(sz./(2^st.input_binning)); 
    
    
    start_y=1;
    for i=st.eigvects_gal(1):st.eigvects_gal(2)
        %index=size(handles.coefs,1)-(i-1);
        index=i;
        %st_tmp(:,start:(start+sz(1)-1))=reshape(handles.coefs(:,index),sz(1).*sz(2),sz(3));
        vol_tmp=reshape(handles.coefs(:,index),sz(1),sz(2),sz(3));
        start_x=1;
        for ii=1:size(vol_tmp,3)
           st_tmp(start_x:start_x+sz(1)-1,start_y:(start_y+sz(2)-1))=vol_tmp(:,:,ii);
           start_x=start_x+sz(1);
        end;
        
        start_y=start_y+sz(2);    
    end;
    figure; tom_imagesc(st_tmp);
else
    sz=sqrt(size(handles.coefs,1));
    for i=st.eigvects_gal(1):st.eigvects_gal(2)
        %index=size(handles.coefs,1)-(i-1);
        index=i;
        st_tmp(:,:,i)=reshape(handles.coefs(:,index),sz,sz);
    end;
    figure; tom_dspcub(st_tmp);
end;






% --- Executes on button press in calc_pca2.
function calc_pca2_Callback(hObject, eventdata, handles)
% hObject    handle to calc_pca2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

st=get_gui_values(handles);

disp('reading stack...');

if (st.flag_3d==1)
    filenames=tom_path2cell(fileparts(st.in_path),st.basename,st.extension,[1 st.last_particle]);
%    st_new=tom_reshape_stack(filenames,'',st.input_binning);
    st_new=1;
else
    hist_flag=0;
    if (hist_flag==0)
        st_new=tom_reshape_stack(st.in_path,'',st.input_binning);
    else
        st_new=tom_calc_hist_stack(st.in_path,20);
    end;
    
 end;


disp('removing mean...');
[st_new handles.mean]=tom_rm_mean(st_new);


disp('calculating pca...');
%load AAA;
%[handles.scores,handles.coefs]=tom_calc_pca(st_new,st.num_of_eigs,'cpca',Align);
[handles.scores,handles.coefs]=tom_calc_pca(st_new,st.num_of_eigs,'cpca');


guidata(hObject, handles);



% --- Executes on button press in define_class.
function define_class_Callback(hObject, eventdata, handles)
% hObject    handle to define_class (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


st=get_gui_values(handles);
axes(handles.img);

[xi,yi]=getline(handles.img,'closed');

% get coordinates
for i=1:size(st.eigvects,1)
    %index=size(handles.scores,1)-(st.eigvects(i)-1);
    index=st.eigvects(i);
    tmp(:,i)=handles.scores(index,:);
end;


in = inpolygon(tmp(:,1),tmp(:,2),xi,yi);
indd=find(in);
hold on; plot(tmp(indd,1),tmp(indd,2),'go'); hold off;


disp('averaging...');

if (st.flag_3d==1)
    h=tom_reademheader(st.in_path);
    avg=zeros(h.Header.Size');
    filenames=tom_path2cell(fileparts(st.in_path),st.basename,st.extension,[1 st.last_particle]);
else
    h=tom_reademheader(st.in_path);
    avg=zeros(h.Header.Size(1),h.Header.Size(2));
end;


for i=1:size(indd,1)
    if (st.flag_3d==1)
        part=tom_emread([filenames{indd(i)}]);
    else
        part=tom_emreadc([st.in_path],'subregion',[1 1 indd(i)],[h.Header.Size(1)-1 h.Header.Size(1)-1 0]);
    end;
    part=part.Value;
    avg=avg+part;
end;


if (st.flag_3d==1)
    figure; tom_dspcub(avg);
else
    figure; tom_imagesc(avg);
end;




% --- Executes on button press in cluster.
function cluster_Callback(hObject, eventdata, handles)
% hObject    handle to cluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

st=get_gui_values(handles);
st.eigvects=(4:40)';

for i=1:size(st.eigvects,1)
    %index=size(handles.scores,1)-(st.eigvects(i)-1);
    index=st.eigvects(i);
    tmp(:,i)=handles.scores(index,:);
end;


flag='k-means';


if (strcmp(flag,'k-means'))
    [classes centriod sum_dist dists] = kmeans(tmp,st.num_of_classes);
end;

if (strcmp(flag,'Kohonen_SOFM'));
    t2=ones(size(tmp,1),1);
    [a b classes]=tom_kohonen_sofm(tmp',t2',[st.num_of_classes 7],0);
end;

if (strcmp(flag,'AGHC'));
    t2=ones(size(tmp,1),1);
    [a b classes]=tom_AGHC(tmp',t2,[st.num_of_classes 'min'],0);
end;

if (strcmp(flag,'leader_follower'));
    t2=ones(size(tmp,1),1);
    [a b classes]=tom_leader_follower(tmp',t2,[0.4 0.8],0);
end;


if (strcmp(flag,'competetive_learning'));
    t2=ones(size(tmp,1),1);
    [a b classes]=tom_leader_follower(tmp',t2',[st.num_of_classes 0.01],0);
end;

if (strcmp(flag,'k-means2'));
    t2=ones(size(tmp,1),1);
    t2(1)=0; t2(10)=0; t(20)=0;
    [a b classes]=k_means(tmp',t2',st.num_of_classes,0);
end;




%if isempty(findobj('Name','pca_3d')) & size(st.eigvects,1) > 2;
%    figure; set(gcf,'Name','pca_3d')
%    hold on; plot3(tmp(:,1),tmp(:,2),tmp(:,3),'.'); hold off;
%end


handles.classes=classes;

%disp class
for ii=1:size(st.eigvects,1)

%    if (ii==2 & size(st.eigvects,1) > 2 )
%        figure(findobj('Name','pca_3d')) ;
%    end;

    for i=1:size(tmp,1)
        if (classes(i)==1)
            col=['go'];
        end;
        if (classes(i)==2)
            col=['ro'];
        end;
        if (classes(i)==3)
            col=['mo'];
        end;
        if (classes(i)==4)
            col=['co'];
        end;
        if (classes(i)==5)
            col=['yo'];
        end;

        if (classes(i)==6)
            col=['yo'];
        end;

        if (classes(i)==7)
            col=['ko'];
        end;

        if (classes(i)>=8)
            col=['ro'];
        end;


        if (size(st.eigvects,1)==3)

            if (ii==2)
                figure(findobj('Name','pca_3d'));
                hold on; plot3(tmp(i,1),tmp(i,2),tmp(i,3),col); hold off;
            else
                axes(handles.img);
                hold on; plot3(tmp(i,1),tmp(i,2),tmp(i,3),col); hold off;
            end;

        else
            hold on; plot(tmp(i,1),tmp(i,2),col); hold off;

        end;

    end;

end;





guidata(hObject, handles);





% --- Executes on button press in opt_clustering.
function opt_clustering_Callback(hObject, eventdata, handles)
% hObject    handle to opt_clustering (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

st=get_gui_values(handles);


for i=1:size(st.eigvects,1)
    %index=size(handles.scores,1)-(st.eigvects(i)-1);
    index=st.eigvects(i);
    tmp(:,i)=handles.scores(index,:);
end;



[dists opt_num_of_classes]=tom_find_opt_num_of_classes(tmp,st.opt_num_of_classes,'k-means');

figure; plot(dists(:,2),dists(:,1));

%set(



% --- Executes on button press in show_class.
function show_class_Callback(hObject, eventdata, handles)
% hObject    handle to show_class (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%disp class




st=get_gui_values(handles);



for i=1:size(st.eigvects,1)
    %index=size(handles.scores,1)-(st.eigvects(i)-1);
    index=st.eigvects(i);
    tmp(:,i)=handles.scores(index,:);
end;

% if isempty(findobj('Name','pca_3d')) & size(st.eigvects,1) > 2;
%     figure; set(gcf,'Name','pca_3d')
%     hold on; plot3(tmp(:,1),tmp(:,2),tmp(:,3),'.'); hold off;
% end


axes(handles.img);

% index1=size(handles.coefs,1)-(st.eigvects(1)-1);
% index2=size(handles.coefs,1)-(st.eigvects(2)-1);

index1=st.eigvects(1);
index2=st.eigvects(2);


if (size(st.eigvects,1)>2)
    %axes(handles.img); set(handles.img,'MenuBar','figure');

    %index3=size(handles.coefs,1)-(st.eigvects(3)-1);
    index3=st.eigvects(3);
    plot3(handles.scores(index1,:),handles.scores(index2,:),handles.scores(index3,:),'r+');
    figure; plot3(handles.scores(index1,:),handles.scores(index2,:),handles.scores(index3,:),'r+');
    set(gcf,'Name','pca_3d');
else
    plot(handles.scores(index1,:),handles.scores(index2,:),'r+');
end;


if (st.flag_3d==1)
    filenames=tom_path2cell(fileparts(st.in_path),st.basename,st.extension,[1 st.last_particle]);
    h=tom_reademheader([st.in_path]);
    avg=zeros(h.Header.Size');
else
    h=tom_reademheader(st.in_path);
    avg=zeros(h.Header.Size(1),h.Header.Size(2));
end;


avg2=zeros(h.Header.Size(1),h.Header.Size(2),st.num_of_classes);
avg2_num=zeros(st.num_of_classes,1);
col='bo';
for ii=1:size(st.eigvects,1)
    for i=1:size(tmp,1)
    %  if (handles.classes(i)==st.class)

      %      if (ii==2  & size(st.eigvects,1) > 2)
      %          figure(findobj('Name','pca_3d'));
      %      end;

            if (size(st.eigvects,1)==3)
                if (ii==2)
                    figure(findobj('Name','pca_3d'));
                    hold on; plot3(tmp(i,1),tmp(i,2),tmp(i,3),col); hold off;
                else
                    axes(handles.img);
                    hold on; plot3(tmp(i,1),tmp(i,2),tmp(i,3),col); hold off;
                end;
            else
                hold on; plot(tmp(i,1),tmp(i,2),col); hold off;
            end;

            if (st.flag_3d==1)
                part=tom_emread(filenames{i});
            else
                part=tom_emreadc([st.in_path],'subregion',[1 1 i],[h.Header.Size(1)-1 h.Header.Size(1)-1 0]);
            end;
            part=part.Value;
            avg=avg+part;
            avg2(:,:,handles.classes(i))=avg2(:,:,handles.classes(i))+part;
            avg2_num(handles.classes(i))=avg2_num(handles.classes(i))+1;
    %    end; %num of classes
    end;
end;

for i=1:st.num_of_classes
    avg2(:,:,i)=avg2(:,:,i)./avg2_num(i);
end;


if  (st.flag_3d==1)
    figure; tom_dspcub(avg);
else
    figure; tom_imagesc(avg);
end;

figure; tom_dspcub(avg2);


% --- Executes on button press in check_3d.
function check_3d_Callback(hObject, eventdata, handles)
% hObject    handle to check_3d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_3d

if (get(handles.check_3d,'Value')==1)
    set(handles.basename,'Visible','on');
    set(handles.extension,'Visible','on');
    set(handles.last_part,'Visible','on');
    set(handles.label_basename,'Visible','on');
    set(handles.label_extension,'Visible','on');
    set(handles.label_last_particle,'Visible','on');
else
    set(handles.basename,'Visible','off');
    set(handles.extension,'Visible','off');
    set(handles.last_part,'Visible','off');
    set(handles.label_basename,'Visible','off');
    set(handles.label_extension,'Visible','off');
    set(handles.label_last_particle,'Visible','off');
end;


function path_input_Callback(hObject, eventdata, handles)
% hObject    handle to path_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of path_input as text
%        str2double(get(hObject,'String')) returns contents of path_input as a double


function cluster_edit_Callback(hObject, eventdata, handles)
% hObject    handle to cluster_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cluster_edit as text
%        str2double(get(hObject,'String')) returns contents of cluster_edit as a double


function edit_show_class_Callback(hObject, eventdata, handles)
% hObject    handle to edit_show_class (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_show_class as text
%        str2double(get(hObject,'String')) returns contents of edit_show_class as a double

function input_binning_Callback(hObject, eventdata, handles)
% hObject    handle to input_binning (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of input_binning as text
%        str2double(get(hObject,'String')) returns contents of input_binning as a double


function last_part_Callback(hObject, eventdata, handles)
% hObject    handle to last_part (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of last_part as text
%        str2double(get(hObject,'String')) returns contents of last_part as a double


function extension_Callback(hObject, eventdata, handles)
% hObject    handle to extension (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of extension as text
%        str2double(get(hObject,'String')) returns contents of extension as a double



function basename_Callback(hObject, eventdata, handles)
% hObject    handle to basename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with % --- Executes during object creation, after setting all properties.


function edit_show_eig_gal_Callback(hObject, eventdata, handles)
% hObject    handle to edit_show_eig_gal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_show_eig_gal as text
%        str2double(get(hObject,'String')) returns contents of edit_show_eig_gal as a double



% function edit_show_class_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to edit_show_class (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: edit controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'Backgrounfunction extension_Callback(hObject, eventdata, handles)
% % hObject    handle to extension (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: get(hObject,'String') returns contents of extension as text
% %        str2double(get(hObject,'String')) returns contents of extension as a double
% dColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end
% 
% handles and user data (see GUIDATA)
% 
% % Hints: get(hObject,'String') returns contents of basename as text
% %        str2double(get(hObject,'String')) returns contents of basename as a double





% --- Executes during object creation, after setting all properties.
function extension_CreateFcn(hObject, eventdata, handles)
% hObject    handle to extension (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% --- Executes during object creation, after setting all properties.
function cluster_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cluster_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function last_part_CreateFcn(hObject, eventdata, handles)
% hObject    handle to last_part (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edit_show_eig_gal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_show_eig_gal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function basename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to basename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function edit_show_class_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_show_class (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
% function cluster_edit_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to cluster_edit (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: edit controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end


% --- Executes during object creation, after setting all properties.
function input_binning_CreateFcn(hObject, eventdata, handles)
% hObject    handle to input_binning (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function path_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to path_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function num_of_eigs_Callback(hObject, eventdata, handles)
% hObject    handle to num_of_eigs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_of_eigs as text
%        str2double(get(hObject,'String')) returns c% --- Executes on button press in opt_clustering.
%function opt_clustering_Callback(hObject, eventdata, handles)
% hObject    handle to opt_clustering (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)ontents of num_of_eigs as a double


% --- Executes during object creation, after setting all properties.
function num_of_eigs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_of_eigs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end








function opt_num_eig_Callback(hObject, eventdata, handles)
% hObject    handle to opt_num_eig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of opt_num_eig as text
%        str2double(get(hObject,'String')) returns contents of opt_num_eig as a double


% --- Executes during object creation, after setting all properties.
function opt_num_eig_CreateFcn(hObject, eventdata, handles)
% hObject    handle to opt_num_eig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function opt_num_classes_Callback(hObject, eventdata, handles)
% hObject    handle to opt_num_classes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of opt_num_classes as text
%        str2double(get(hObject,'String')) returns contents of opt_num_classes as a double


% --- Executes during object creation, after setting all properties.
function opt_num_classes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to opt_num_classes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


