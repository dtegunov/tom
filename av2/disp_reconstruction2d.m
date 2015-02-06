function varargout = disp_reconstruction2d(varargin)
%TOM_DISP_RECONSTRUCTION2D is a object for alignment demo mode
%
%   varargout = disp_reconstruction2d(varargin)
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
%   ... = disp_reconstruction2d(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by AK 12/01/06
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
                   'gui_OpeningFcn', @disp_reconstruction2d_OpeningFcn, ...
                   'gui_OutputFcn',  @disp_reconstruction2d_OutputFcn, ...
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
function disp_reconstruction2d_OpeningFcn(hObject, eventdata, handles, varargin)

global storage_rec2d;

func = varargin{1};


switch func
    
    case 'init'
        numimages = varargin{2};
        numclasses = varargin{3};
        imsize = varargin{4};
        columns = varargin{5};
        storage_rec2d.numclasses = numclasses;
        storage_rec2d.histogram.data = [];
        storage_rec2d.projections.data = zeros(imsize(1), imsize(2), numimages);
        storage_rec2d.classprojections.data = zeros(imsize(1), imsize(2), numclasses);
        storage_rec2d.projections.numcolumns = columns;
        storage_rec2d.projections.imsize = imsize;
        storage_rec2d.classprojections.imsize = imsize;
        storage_rec2d.classprojections.numcolumns = columns;
        storage_rec2d.projections.maximage = 0;
        display_projections();
        display_clsprojections();
        
        axes(findobj('Tag','reconst'));axis off;
        axes(findobj('Tag','matching_particle'));axis off;
        axes(findobj('Tag','matching_projection'));axis off;
        axes(findobj('Tag','matching_xcorrelation'));axis off;
        axes(findobj('Tag','matchresult_particle'));axis off;
        axes(findobj('Tag','matchresult_projection'));axis off;
        axes(findobj('Tag','im_template'));axis off;
        axes(findobj('Tag','im_isosurf'));axis off;
        set(findobj('Tag','histogram'),'XLim',[1 numclasses]);
        
    case 'addprojection'
        image = varargin{2};
        addprojection(image);
    
    case 'clsadd'
        image = varargin{2};
        number = varargin{3};
        clsadd(image, number);     

    case 'mark_projection'
        number = varargin{2};
        mark_projection(number);
        
    case 'mark_class'
        number = varargin{2};
        mark_class(number);
    
    case 'display_match'
        im1 = varargin{2};
        im2 = varargin{3};
        im3 = varargin{4};
        title_in = varargin{5};
        peak = varargin{6};
        pos_flag = varargin{7};
        display_match(im1,im2,im3,title_in,peak,pos_flag);
    
    case 'update_hist'
        clsidx = varargin{2};
        update_hist(clsidx);
    
    case 'disp_reconst'
        vol = varargin{2};
        disp_reconst(vol);
    
    case 'init_surface'
        vol = varargin{2};
        vol = tom_rotate(single(vol),[90 0 0]);
        storage_rec2d.thresh = varargin{3};
        storage_rec2d.vol = vol;
        tmp_obj=findobj('Tag','im_isosurf');
        axes(tmp_obj);cla;axis vis3d;axis off;
        t = hgtransform('Parent',tmp_obj,'Tag','hgt');
        p = patch(isosurface(storage_rec2d.vol, storage_rec2d.thresh),'FaceColor', [230./255 213./255 17./255], 'EdgeColor', 'none','Tag','iso','Parent',t);
        isonormals(storage_rec2d.vol,p);
        
    case 'update_surface'
        
        euler_in = varargin{2};
        
        update_surface(euler_in);
        
    case 'destroy'
        clear storage_rec2d;
        tmpobj = findobj('Tag','disp_rec2d');
        if tmpobj ~= 0
            delete(tmpobj);
        end
        return;
end

% Choose default command line output for disp_reconstruction2d
handles.output = hObject;
guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Output Function                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = disp_reconstruction2d_OutputFcn(hObject, eventdata, handles) 
%varargout{1} = handles.output;


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
%%  Add image to projection viewer                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function addprojection(image)

global storage_rec2d;

storage_rec2d.projections.maximage = storage_rec2d.projections.maximage + 1;
storage_rec2d.projections.data(:,:,storage_rec2d.projections.maximage) = image;

tmpobj = findobj('Tag','im_template');
axes(tmpobj);
imagesc(squeeze(storage_rec2d.projections.data(:,:,storage_rec2d.projections.maximage))');colormap gray;axis ij; axis off;
set(tmpobj,'Tag','im_template');

display_projections();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display projections                                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function display_projections()

global storage_rec2d;

tmpobj = findobj('Tag','image_projclasses');
axes(tmpobj);
tom_dspcub(storage_rec2d.projections.data,0, storage_rec2d.projections.numcolumns);
set(tmpobj,'Tag','image_projclasses');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  mark projection                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mark_projection(number)

global storage_rec2d;

tmpobj = findobj('Tag','proj_rect');
if tmpobj ~= 0
    delete(tmpobj);
end

tmpobj = findobj('Tag','im_template');
axes(tmpobj);
imagesc(squeeze(storage_rec2d.projections.data(:,:,number))');colormap gray;axis ij; axis off;
set(tmpobj,'Tag','im_template');

tmpobj = findobj('Tag','image_projclasses');
axes(tmpobj);
set(tmpobj,'Units','Pixel');
number = number - 1;
pos_y = floor(number ./ storage_rec2d.projections.numcolumns) * storage_rec2d.projections.imsize(2);
pos_x = (mod(number, storage_rec2d.projections.numcolumns)) * storage_rec2d.projections.imsize(1);
width = storage_rec2d.projections.imsize(1);
height = storage_rec2d.projections.imsize(2);
rectangle('Position',[pos_x, pos_y, width, height],'Tag','proj_rect','EdgeColor',[1 0 0],'LineStyle','-','LineWidth',2);
set(tmpobj,'Tag','image_projclasses');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Add image to class viewer                                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function clsadd(image, number)

global storage_rec2d;

storage_rec2d.classprojections.data(:,:,number) = image;
display_clsprojections();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display class projections                                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function display_clsprojections()

global storage_rec2d;

tmpobj = findobj('Tag','image_classes');
axes(tmpobj);
tom_dspcub(storage_rec2d.classprojections.data,0, storage_rec2d.classprojections.numcolumns);
set(tmpobj,'Tag','image_classes');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  mark class                                                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mark_class(number)

global storage_rec2d;

tmpobj = findobj('Tag','classproj_rect');
if tmpobj ~= 0
    delete(tmpobj);
end

update_hist(number);

tmpobj = findobj('Tag','image_classes');
axes(tmpobj);
set(tmpobj,'Units','Pixel');
number = number - 1;
pos_y = floor(number ./ storage_rec2d.classprojections.numcolumns) * storage_rec2d.classprojections.imsize(2);
pos_x = (mod(number, storage_rec2d.classprojections.numcolumns)) * storage_rec2d.classprojections.imsize(1);
width = storage_rec2d.projections.imsize(1);
height = storage_rec2d.projections.imsize(2);
rectangle('Position',[pos_x, pos_y, width, height],'Tag','classproj_rect','EdgeColor',[0 1 0],'LineStyle','-','LineWidth',2);
set(tmpobj,'Tag','image_classes');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display                                                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function display_match(im1,im2,im3,title,peak,pos_flag)

if strcmp(pos_flag,'upper')
    
    set(findobj('Tag','upper_1'),'String',title{1});
    set(findobj('Tag','upper_2'),'String',title{2});
    set(findobj('Tag','upper_3'),'String',title{3});
    
    tmpobj = findobj('Tag','matching_particle');
    axes(tmpobj);
    imagesc(tom_norm(im1',1),[0,1]);colormap(gray);axis ij;axis off;
    set(tmpobj,'Tag','matching_particle');

    tmpobj = findobj('Tag','matching_projection');
    axes(tmpobj);
    imagesc(tom_norm(im2',1),[0,1]);colormap(gray);axis ij;axis off;
    set(tmpobj,'Tag','matching_projection');

    tmpobj = findobj('Tag','matching_xcorrelation');
    axes(tmpobj);
    imagesc(tom_norm(im3',1),[0,1]);colormap(gray);axis ij;axis off;
    hold on; plot(peak(1),peak(2),'ro'); hold off; drawnow;
    set(tmpobj,'Tag','matching_xcorrelation');
    
else
    
    set(findobj('Tag','lower_1'),'String',title{1});
    set(findobj('Tag','lower_2'),'String',title{2});
    set(findobj('Tag','lower_3'),'String',title{3});
    
    tmpobj = findobj('Tag','matchresult_particle');
    axes(tmpobj);
    imagesc(tom_norm(im1',1),[0,1]);colormap(gray);axis ij;axis off;
    set(tmpobj,'Tag','matchresult_particle');

    tmpobj = findobj('Tag','matchresult_projection');
    axes(tmpobj);
    imagesc(tom_norm(im2',1),[0,1]);colormap(gray);axis ij;axis off;
    set(tmpobj,'Tag','matchresult_projection');

    tmpobj = findobj('Tag','matchresult_xcorrfunction');
    axes(tmpobj);
    plot(im3);
    hold on; plot(peak(1),peak(2),'ro'); hold off; drawnow;
    set(tmpobj,'Tag','matchresult_xcorrfunction');
        
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  update histogram                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function update_hist(idx)

global storage_rec2d;

storage_rec2d.histogram.data = [storage_rec2d.histogram.data, idx];

tmpobj = findobj('Tag','histogram');
axes(tmpobj);
hist(storage_rec2d.histogram.data);
set(tmpobj,'Tag','histogram','XLim',[1 storage_rec2d.numclasses]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  display reconstruction                                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function disp_reconst(volume)

tmpobj = findobj('Tag','reconst');
axes(tmpobj);
subplot(2,2,3);
%set(gca,'position',[0.33 0.05 0.327023 0.5]);
axis ij;
a=round(size(volume,3)./2)-5;
b=round(size(volume,3)./2)+5;
t_vol=volume(:,:,a:b);
ixy=sum(t_vol,3);
a=round(size(volume,1)./2)-5;
b=round(size(volume,1)./2)+5;
t_vol=volume(a:b,:,:);
iyz=(squeeze(sum(t_vol,1))');
a=round(size(volume,1)./2);%a=round(size(volume,1)./2)-5;
b=round(size(volume,1)./2);%b=round(size(volume,1)./2)+5;
t_vol=volume(:,a:b,:);
ixz=(squeeze(sum(t_vol,2)));
all=zeros(size(ixy,1)+20+size(iyz,1),size(ixy,2)+20+size(ixz,2));
all(1:size(ixy,1),1:size(ixy,2))=ixy(:,:);
all(size(ixy,1)+21:size(ixy,1)+20+size(iyz,1),1:size(iyz,2))=iyz(:,:);
all(1:size(ixz,1),size(ixy,2)+21:size(ixy,2)+20+size(ixz,2))=ixz(:,:);
imagesc(all);
set(gca,'Tag','reconst');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  draw isosurface                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function update_surface(euler_in)

global storage_rec2d;

tmp_obj=findobj('Tag','im_isosurf');
axes(tmp_obj);axis vis3d;axis off;
t = findobj('Tag','hgt');
m = eye(4);
set(t,'Matrix',m);
m = makehgtform('zrotate',euler_in(1).*pi./180,'yrotate',euler_in(2).*pi./180,'xrotate',euler_in(3).*pi./180);
set(t,'Matrix',m);
set(gcf,'Renderer','opengl');
view([0 90]); camlight; light;daspect('auto');
lighting phong;
