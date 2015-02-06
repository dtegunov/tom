function varargout = disp_reconstruction2d(varargin)
%GUI object for alignment demo mode
%
%
%
%Copyright (c) 2005
%TOM toolbox for Electron Tomography
%Max-Planck-Institute for Biochemistry
%Dept. Molecular Structural Biology
%82152 Martinsried, Germany
%http://www.biochem.mpg.de/tom
%
%Created: 12/01/06 AK
%
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
set(gcf,'Renderer','opengl');

switch func
    
    case 'init'
        storage_rec2d.isoisinited = 0;
        storage_rec2d.showstackdummy = 0;
        try
            a = load('/fs/sally/pool-baumeister/demo/demofilm4/frame.mat');
        catch
            a.frameno = 1;
        end
        storage_rec2d.movieframeno = a.frameno;
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
        storage_rec2d.isoisinited = 1;
        set(findobj('Tag','panel_isosurf'),'Visible','on');
        mask = tom_spheremask(ones(size(varargin{2})),size(varargin{2},1)-size(varargin{2},1)./10);
        storage_rec2d.vol = -varargin{2} .* mask;
        storage_rec2d.thresh = varargin{3};
        storage_rec2d.views = varargin{4};
        tmp_obj=findobj('Tag','im_isosurf');
        axes(tmp_obj);cla;axis vis3d;axis off;
        for i=1:size(storage_rec2d.views,2)
            eu=tom_sum_rotation([0 0 storage_rec2d.views(2,i); 270 90 storage_rec2d.views(1,i)],[0 0 0 ; 0 0 0]);
            tmpvol = tom_rotate(storage_rec2d.vol,[0 0 storage_rec2d.views(2,i)]);
            tmpvol = tom_norm(tmpvol,1);
            isonormals(tmpvol,patch(isosurface(tmpvol, storage_rec2d.thresh),'FaceColor', [230./255 213./255 17./255], 'EdgeColor', 'none','Tag',['patch_' num2str(i)],'AmbientStrength',0.5));
            set(findobj('Tag',['patch_' num2str(i)]),'Visible','off');            
        end
        
    case 'update_surface'
        
        euler_in = varargin{2};
        
        update_surface(euler_in);
        
    case 'destroy'
        frameno=storage_rec2d.movieframeno;
        save('/fs/sally/pool-baumeister/demo/demofilm4/frame.mat','frameno');
        clear storage_rec2d;
        tmpobj = findobj('Tag','disp_rec2d');
        if tmpobj ~= 0
            delete(tmpobj);
        end
        return;
        
    case 'movieopen'
        %Prepare the avi handle
%         storage_rec2d.mov = avifile(varargin{2});
%         storage_rec2d.mov.Fps=varargin{3};
%         if isequal(computer,'PCWIN')
%             storage_rec2d.mov.Compression='Cinepak';
%         else
%             storage_rec2d.mov.Compression='none';
%         end
        % 'Cinepak', 'Indeo3', 'Indeo5', 'MSVC',', 'RLE', 'None'
%        storage_rec2d.mov.Quality=varargin{4};
    case 'moviepicture'
        unix(['import -window disp_reconstruction2d -silent /fs/sally/pool-baumeister/demo/demofilm4/demo_' sprintf('%05.0f',storage_rec2d.movieframeno) '.tif']);
        storage_rec2d.movieframeno = storage_rec2d.movieframeno + 1;
        %set(gcf,'Renderer','zbuffer');
        %F = getframe(gcf);
        %storage_rec2d.mov = addframe(storage_rec2d.mov,F);
        %set(gcf,'Renderer','opengl');
    case 'moviesave'
%         storage_rec2d.mov = close(storage_rec2d.mov); 
        
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%mod = tom_emreadc(['model/model_' num2str(i)]);
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

update_hist(number);

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

global storage_rec2d;

if strcmp(pos_flag,'upper')
    
    set(findobj('Tag','upper_1'),'String',title{1});
    set(findobj('Tag','upper_2'),'String',title{2});
    set(findobj('Tag','upper_3'),'String',title{3});
    
    tmpobj = findobj('Tag','matching_particle');
    axes(tmpobj);
    imagesc(im1');colormap(gray);axis ij;axis off;
    set(tmpobj,'Tag','matching_particle');

    tmpobj = findobj('Tag','matching_projection');
    axes(tmpobj);
    imagesc(im2');colormap(gray);axis ij;axis off;
    set(tmpobj,'Tag','matching_projection');

    tmpobj = findobj('Tag','matching_xcorrelation');
    axes(tmpobj);
    imagesc(im3');colormap(gray);axis ij;axis off;
    hold on; plot(peak(1),peak(2),'ro'); hold off; drawnow;
    set(tmpobj,'Tag','matching_xcorrelation');
    
else
    set(findobj('Tag','matchresult_xcorrfunction'),'Visible','on');
    if storage_rec2d.showstackdummy == 0
        for i=4:-1:1
            set(findobj('Tag',['axes_dummy' num2str(i)]),'Visible','on');
            axes(findobj('Tag',['axes_dummy' num2str(i)]));
            if i > 1
                imagesc(imread('border.png'));axis off;
            else
                imagesc(imread('border2.png'));axis off;
            end
            storage_rec2d.showstackdummy = 1;
        end
    end
    
    set(findobj('Tag','lower_1'),'String',title{1});
    set(findobj('Tag','lower_2'),'String',title{2});
    set(findobj('Tag','lower_3'),'String',title{3});
    
    tmpobj = findobj('Tag','matchresult_particle');
    axes(tmpobj);
    imagesc(im1');colormap(gray);axis ij;axis off;
    set(tmpobj,'Tag','matchresult_particle');

    tmpobj = findobj('Tag','matchresult_projection');
    axes(tmpobj);
    imagesc(im2');colormap(gray);axis ij;axis off;
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
function disp_reconst(stack)

global storage_rec2d;

set(findobj('Tag','panel_reconstruction'),'Visible','on');

tmpobj = findobj('Tag','reconst');

for ii=1:1
    for i=1:size(stack,3)
        axes(tmpobj);
        imagesc(stack(:,:,i)');axis ij;axis off; colormap gray;
        unix(['import -window disp_reconstruction2d -silent /fs/sally/pool-baumeister/demo/demofilm4/demo_' sprintf('%05.0f',storage_rec2d.movieframeno) '.tif']);
        storage_rec2d.movieframeno = storage_rec2d.movieframeno + 1;
    end
end
% subplot(2,2,3);
% %set(gca,'position',[0.33 0.05 0.327023 0.5]);
% axis ij;
% a=round(size(volume,3)./2)-5;
% b=round(size(volume,3)./2)+5;
% t_vol=volume(:,:,a:b);
% ixy=sum(t_vol,3);
% a=round(size(volume,1)./2)-5;
% b=round(size(volume,1)./2)+5;
% t_vol=volume(a:b,:,:);
% iyz=(squeeze(sum(t_vol,1))');
% a=round(size(volume,1)./2);%a=round(size(volume,1)./2)-5;
% b=round(size(volume,1)./2);%b=round(size(volume,1)./2)+5;
% t_vol=volume(:,a:b,:);
% ixz=(squeeze(sum(t_vol,2)));
% all=zeros(size(ixy,1)+20+size(iyz,1),size(ixy,2)+20+size(ixz,2));
% all(1:size(ixy,1),1:size(ixy,2))=ixy(:,:);
% all(size(ixy,1)+21:size(ixy,1)+20+size(iyz,1),1:size(iyz,2))=iyz(:,:);
% all(1:size(ixz,1),size(ixy,2)+21:size(ixy,2)+20+size(ixz,2))=ixz(:,:);
% imagesc(all);

set(gca,'Tag','reconst');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  draw isosurface                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function update_surface(idx)

global storage_rec2d;

if storage_rec2d.isoisinited == 1
    tmp_obj=findobj('Tag','im_isosurf');
    axes(tmp_obj);axis vis3d;axis off;
    set(findobj('-regexp','Tag','patch_*'),'Visible','off');      
    set(findobj('Tag',['patch_' num2str(idx)]),'Visible','on');
    view([90 90]); 
    light('Position',[1 0 0],'Style','infinite');
    light('Position',[-1 0 0],'Style','infinite');
    lighting phong;
end