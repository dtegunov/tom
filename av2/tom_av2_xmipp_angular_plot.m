function tom_av2_xmipp_angular_plot(doc_filename,vol,flag,base_rad,sp_radius,offset,samp)
%TOM_AV2_XMIPP_ANGULAR_PLOT displays angular distribution
%
%   tom_av2_xmipp_angular_plot(doc_filename,vol_filename)
%
%PARAMETERS
%
%  INPUT
%   iput_doc             name of doc or struct
%   vol                 (opt.) volume 2 be displayed (filename or matrix)
%   flag                ('chimera'),'surfc','mesh' or 'ps'  
%   base_rad            (opt.) triangle size for the PS/bild file
%   sp_radius           (opt.) sphere radius for the PS/bild file (normally size of volume)  
%   offset              (opt.) shift coordinates center for (chim) bild (normally size of volume)  
%   samp                (opt.) sampling 4 surfc or mesh   
%
%
%
%  OUTPUT
%
%EXAMPLE
%
% tom_av2_xmipp_angular_plot('Iter_12_current_angles_1_cut.doc','test_vol.vol','chimera',4,240,240);
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by FB 01/24/06
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

if (nargin < 2)
    vol='';
end;

if (nargin < 3)
    flag='chimera';
end;

if (nargin < 4)
    base_rad=6;
end;

if (nargin < 5)
    sp_radius=100;
end;

if (nargin < 6)
    offset=100;
end;

if (nargin < 7)
    samp=0.01;
    if (strcmp(flag,'mesh') || strcmp(flag,'meshp'))
        samp=0.04;  
    end;
end;


if (isnumeric(vol))
    tom_spiderwrite('tvolxx4plot.spi',vol);
    vol='tvolxx4plot.spi';
end;


if (isempty(vol)==0)
    vtmp=tom_spiderread(vol);
    sp_radius=round(size(vtmp.Value,1));
    offset=round(size(vtmp.Value,1)./2);
end;

if (isstruct(doc_filename)==0)
    doc_in=tom_xmippdocread(doc_filename);
else
    doc_in=doc_filename;
end;



if (isfield(doc_in,'weight')==0)
    tmp_doc=['tangxx4plot.doc'];
    if (strcmp(flag,'chimera') ||  strcmp(flag,'ps'))
        st_tmp=tom_av2_xmipp_doc2angcldoc(doc_in,tmp_doc);
    else
        st_tmp=tom_av2_xmipp_doc2angcldoc(doc_in);
    end;
    doc_filename=tmp_doc;
    
end;


if (strcmp(flag,'chimera') || strcmp(flag,'ps'))
    
    params=[' -shift_center ' num2str(offset) ' -R ' num2str(sp_radius)  ' -r ' num2str(base_rad) ' -wcol 5'];
    if ( strcmp(flag,'ps'))
        params=[params ' -ps tangxx4plot.ps'];
    end;
    call_base='xmipp_angular_distribution_show -ang ';
    
    [a b c]=fileparts(doc_filename);
    
    call=[call_base doc_filename ' -bild ' a b '.bild'  params ];
    disp(call);
    unix(call);
    
    base_call='chimera ';
    
    if (strcmp(flag,'chimera'))
        if (isempty(vol))
            call=[base_call a b '.bild &'];
        else
            call=[base_call a b '.bild ' vol ' &'];
        end;
    end;
    
    if (strcmp(flag,'ps'))
        call=['okular tangxx4plot.ps'];
    end;
    
    disp(call);
    unix(call);
end;

if (strcmp(flag,'surfc') ||  strcmp(flag,'contourf') || strcmp(flag,'mesh') || strcmp(flag,'meshp'))
    
    
    %transfer variables
    refnum=[st_tmp(:).weight];
    ang_classes = [[st_tmp(:).rot] ; [st_tmp(:).tilt]]';
    ang_classes = ang_classes.*pi/180;
    
    %transform 2 cart coord sys
    xx = ang_classes(:,2).*cos(ang_classes(:,1));
    yy = ang_classes(:,2).*sin(ang_classes(:,1));
    ang_classes = [ang_classes,refnum'];
    
    
    min_x = min(xx);
    max_x = max(xx);
    min_y = min(yy);
    max_y = max(yy);
    [XI,YI] = meshgrid(min_x:samp:max_x,min_y:samp:max_y);
    %ZI = griddata(xx,yy,ang_classes(:,3),XI,YI);
    F = TriScatteredInterp(xx,yy,refnum');
    ZI=F(XI,YI); 

    
    sampling = unique(ang_classes(:,2))./pi.*180;
    rticks = 90./sampling(2);
    
    if (strcmp(flag,'mesh'))
        if (strcmp(flag,'mesh'))
            figure; mesh(ZI); axis image;
        else
            rr=(ZI./1000);
            figure; mesh(XI,YI,rr); axis image;
            hold on;
            rr=refnum./1000;
            plot3(xx,yy,rr,'ko');
            hold off;
        end;
     end;
    
    if (strcmp(flag,'surfc'))
         %figure;surfc(XI,YI,double(ZI)); set(gcf,'DoubleBuffer','on'); shading interp; colormap hot; 
         figure;surfc(double(ZI)); set(gcf,'DoubleBuffer','on'); shading interp; colormap hot; axis image;
    end;
    if (strcmp(flag,'contourf'))
        figure('color',[0 0 0],'Position',[2565 481 1086 1040]);cax = axes;axis off;hold on;
        set(gcf,'DoubleBuffer','on');
        contourf(XI,YI,ZI,200);shading flat;

        hold(cax,'on');
        maxrho = max_y;
        hhh=line([-maxrho -maxrho maxrho maxrho],[-maxrho maxrho maxrho -maxrho],'parent',cax);
        set(cax,'dataaspectratio',[1 1 1],'plotboxaspectratiomode','auto')
        v = [get(cax,'xlim') get(cax,'ylim')];
        ticks = sum(get(cax,'ytick')>=0);
        delete(hhh);
        rmin = 0; rmax = v(4); 
        
        % define a circleglobal reflookup;
        tc = [1 1 1];
        th = 0:pi/50:2*pi;
        ls = get(cax,'gridlinestyle');
        xunit = cos(th);
        yunit = sin(th);
        % now really force points on x/y axes to lie on them exactly
        inds = 1:(length(th)-1)/4:length(th);
        xunit(inds(2:2:4)) = zeros(2,1);
        yunit(inds(1:2:5)) = zeros(3,1);
        c82 = cos(82*pi/180);
        s82 = sin(82*pi/180);
        rinc = (rmax-rmin)/rticks;
        l = 2;
        
        for i=(rmin+rinc):rinc:rmax
%             hhh = line(xunit*i,yunit*i,'linestyle',ls,'color',tc,'linewidth',1,...
%                 'handlevisibility','off','parent',cax);
%             text((i+rinc/20)*c82,(i+rinc/20)*s82, ...
%                     ['  ' num2str(sampling(l))],'verticalalignment','bottom',...
%                     'handlevisibility','off','parent',cax,'color',tc)
%             l = l + floor(length(sampling./11));
        end
        %set(hhh,'linestyle','-') % Make outer circle solid
        % plot spokes
        th = (1:6)*2*pi/12;
        cst = cos(th); snt = sin(th);
        cs = [-cst; cst];
        sn = [-snt; snt];
        %line(rmax*cs,rmax*sn,'linestyle',ls,'color',tc,'linewidth',1,...
        %    'handlevisibility','off','parent',cax)
        
        rt = 1.1*rmax;
%         for i = 1:length(th)
%             text(rt*cst(i),rt*snt(i),int2str(i*30),...
%                 'horizontalalignment','center',...
%                 'handlevisibility','off','parent',cax,'color',tc);
%             if i == length(th)
%                 loc = int2str(0);
%             else
%                 loc = int2str(180+i*30);
%             end
%             text(-rt*cst(i),-rt*snt(i),loc,'horizontalalignment','center',...
%                 'handlevisibility','off','parent',cax,'color',tc)
%         end
        
        % set view to 2-D
        view(cax,2);
        % set axis limits
        axis(cax,rmax*[-1 1 -1.15 1.15]);
        
        view(cax,2);
     end;
    
end;

%clean up
% if (strcmp(flag,'chimera') ||  strcmp(flag,'ps'))   
%     unix('rm tangxx4plot.bild');
%     unix('rm tangxx4plot.doc');
%     if (strcmp(flag,'ps'))
%         unix('rm tangxx4plot.ps');
%     end;
%     try 
%         [a b]=unix('rm  tvolxx4plot.spi');
%     catch Me
%     end;
% end;


