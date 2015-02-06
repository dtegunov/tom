function tom_angular_plot(docfile,xmppath,anno_flag)

%tom_angular_plot shows the angular distribution of the projection class
%   averages. Precession and Nutation angles are shown. The in-plane
%   rotation angle is not displayed.
%
%   tom_angular_plot(doc_file,xmppath)
%
%PARAMETERS
%
%  INPUT
%   doc_file            Xmipp-Doc File of Xmipp Projection Matching
%                       procedure
%   xmppath             path to ReferenceLibrary
%   anno_flag           (0) 1 for full annotation 
%                          -1 no red dots
%
%  OUTPUT
%   interactive angular plot - select angle to show corresponding projection class
%
%EXAMPLE
%  tom_angular_plot('Iter_10_current_angles.doc','./ReferenceLibrary/ref0')
%
%REFERENCES
%
%SEE ALSO
%
%   created by AK
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


if (nargin<3)
   anno_flag=0; 
end;


st = tom_xmippdocread(docfile);

if (isfield(st,'ref'))
    maxref = max([st.ref]);
else
    maxref = length(st);
end;

reflookup = zeros(maxref,3);
for i=1:maxref
    if (isfield(st,'ref'))
        idx = find([st.ref]==i,1);
    else
        idx = i;
    end;
    
    if ~isempty(idx)
        ang1 = st(idx(1)).rot;
        ang2 = st(idx(1)).tilt;
    else
        ang1 = NaN;
        ang2 = NaN;
    end
    reflookup(i,:) = [i,ang1,ang2];
end


rot = [st.rot];
tilt = [st.tilt];

angs = [rot;tilt];

clear rot
clear tilt;

ang_classes = unique(angs', 'rows');

refnum = zeros(size(ang_classes,1),1);

for i=1:size(angs,2)
    angpair = angs(:,i);
    for j=1:size(ang_classes,1)
        if abs(ang_classes(j,1) - angpair(1)) < 0.0001 && abs(ang_classes(j,2) - angpair(2)) < 0.0001
            refnum(j) = refnum(j)+1;
        end
    end
end

ang_classes = ang_classes.*pi/180;

%[pointidx, pointcoords, distance] = tom_nearestpoint([0 0],ang_classes);
%refnum(pointidx)=300;

xx = ang_classes(:,2).*cos(ang_classes(:,1));
yy = ang_classes(:,2).*sin(ang_classes(:,1));


ang_classes = [ang_classes,refnum];


min_x = min(xx);
max_x = max(xx);
min_y = min(yy);
max_y = max(yy);

[XI,YI] = meshgrid(min_x:0.01:max_x,min_y:0.01:max_y);
ZI = griddata(xx,yy,ang_classes(:,3),XI,YI);

sampling = unique(ang_classes(:,2))./pi.*180;
rticks = 90./sampling(2);


figure('color',[0 0 0],'Position',[2565 481 1086 1040]);cax = axes;axis off;hold on;
set(gcf,'DoubleBuffer','on');
contourf(XI,YI,ZI,200);shading flat;
hold on;

rohandle = zeros(length(xx),1);


%tc = get(cax,'xcolor');
tc = [1 1 1];

ls = get(cax,'gridlinestyle');
% make a radial grid
    hold(cax,'on');
    maxrho = max_y;
    hhh=line([-maxrho -maxrho maxrho maxrho],[-maxrho maxrho maxrho -maxrho],'parent',cax);
    set(cax,'dataaspectratio',[1 1 1],'plotboxaspectratiomode','auto')
    v = [get(cax,'xlim') get(cax,'ylim')];
    ticks = sum(get(cax,'ytick')>=0);
    delete(hhh);
% check radial limits and ticks
    rmin = 0; rmax = v(4); %rticks = max(ticks-1,2);
%     if rticks > 5   % see if we can reduce the number
%         if rem(rticks,2) == 0
%             rticks = rticks/2;
%         elseif rem(rticks,3) == 0
%             rticks = rticks/3;
%         end
%     end

% define a circleglobal reflookup;
    th = 0:pi/50:2*pi;
    xunit = cos(th);
    yunit = sin(th);
% now really force points on x/y axes to lie on them exactly
    inds = 1:(length(th)-1)/4:length(th);
    xunit(inds(2:2:4)) = zeros(2,1);
    yunit(inds(1:2:5)) = zeros(3,1);
% plot background if necessary
    if ~ischar(get(cax,'color')),
       %patch('xdata',xunit*rmax,'ydata',yunit*rmax, ...
       %      'edgecolor',tc,'facecolor',get(cax,'color'),...
       %      'handlevisibility','off','parent',cax);
    end

% draw radial circles
    c82 = cos(82*pi/180);
    s82 = sin(82*pi/180);
    rinc = (rmax-rmin)/rticks;
    l = 2;
    for i=(rmin+rinc):rinc:rmax
        hhh = line(xunit*i,yunit*i,'linestyle',ls,'color',tc,'linewidth',1,...
                   'handlevisibility','off','parent',cax);
        
        if (anno_flag==1)       
        text((i+rinc/20)*c82,(i+rinc/20)*s82, ...
            ['  ' num2str(sampling(l))],'verticalalignment','bottom',...
            'handlevisibility','off','parent',cax,'color',tc)
        end;

        l = l + 1;
    end
    set(hhh,'linestyle','-') % Make outer circle solid

% plot spokes
    th = (1:6)*2*pi/12;
    cst = cos(th); snt = sin(th);
    cs = [-cst; cst];
    sn = [-snt; snt];
    line(rmax*cs,rmax*sn,'linestyle',ls,'color',tc,'linewidth',1,...
         'handlevisibility','off','parent',cax)

% annotate spokes in degrees
  

    rt = 1.1*rmax;
    for i = 1:length(th)
        text(rt*cst(i),rt*snt(i),int2str(i*30),...
             'horizontalalignment','center',...
             'handlevisibility','off','parent',cax,'color',tc);
        if i == length(th)
            loc = int2str(0);
        else
            loc = int2str(180+i*30);
        end
        text(-rt*cst(i),-rt*snt(i),loc,'horizontalalignment','center',...
             'handlevisibility','off','parent',cax,'color',tc)
    end

% set view to 2-D
    view(cax,2);
% set axis limits
    axis(cax,rmax*[-1 1 -1.15 1.15]);

view(cax,2);

if (anno_flag==-1)
    return;
end;

for i=1:length(xx)
    pointidx = tom_nearestpoint(ang_classes(i,1:2)./pi*180, reflookup(:,2:3));
    rohandle(i) = plot(xx(i),yy(i),'r*','Userdata',{pointidx,xmppath});
end
   
for i=1:length(xx)

    set(rohandle(i),'ButtonDownFcn',@button_down);
end


function button_down(src,evnt)


ud = get(gco,'UserData');

%im = tom_spiderread([ud{2} sprintf('%05.0f',ud{1}) '.proj']);
im = tom_spiderread([ud{2} sprintf('%05.0f',ud{1}) '.xmp']);
figure;tom_imagesc(im.Value)

% src - the object that is the source of the event
% evnt - empty for this property
%    sel_typ = get(gcbf,'SelectionType');
%    switch sel_typ 
%       case 'normal'
%          disp('User clicked left-mouse button')
%          set(src,'Selected','on')
%       case 'extend'
%          disp('User did a shift-click')
%          set(src,'Selected','on')
%       case 'alt'
%          disp('User did a control-click')
%          set(src,'Selected','on')
%          set(src,'SelectionHighlight','off')
%    end
