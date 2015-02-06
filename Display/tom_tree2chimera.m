function part_ord=tom_tree2chimera(tree,volume_list,thresholds,spacing_fact,colormap,tmp_dir)
%TOM_TREE2CHIMERA displays a tree from matlab linkage in Chimera
%  
%
%   tom_tree2chimera(tree,volume_list)
%
%
%  INPUT
%   tree           tree matrix from in linkage style
%   volume_list    (opt.) as sel  ...file order !!   
%   thresholds     ('mean+1std') threshold for iso rep
%   spacing_fact   (1.5) factor for offset between models
%   map4Score      ('quadr') maps score of Y-Axis 2 occupy the space as x
%                            to optain a quadratic field of view
%   colormap       ([1 0 0]) colormpa or maps for the models
%   tmp_dir        ('./xxx_tmp4dend_xxx') directrory where tmp files are stored  
%  
%
%  OUTPUT
%   part_ord     particle order (first line of 3d dendogram)
%
%EXAMPLE
%  
%  mkdir('4testXX');
%  sp_rad=[4 4 4 8 8 8];
%  for i=1:6
%    sp=tom_spheremask(ones(64,64,64),sp_rad(i))+(0.1.*rand(64,64,64)) ;
%    tom_emwrite(['4testXX/v_' num2str(i) '.em'],sp);
%  end;
%  [cc vol_list] = tom_av3_ccmat('4testXX/v','em',6);
%  tree = linkage(squareform(1-cc),'single');
%  tom_tree2chimera(tree,vol_list,0.28,2.1);
%  dendrogram(tree);
%
%REFERENCES
%
%SEE ALSO
%     
%  linkage,tom_av3_ccmat
%   
% NOTE:
% 
% tree should be in this format: 
% (linkage from: MATLAB Version: 8.0.0.783 (R2012b))
% 
% 
% tree =
% 
%     5.0000    6.0000    0.0944
%     4.0000    7.0000    0.0948
%     1.0000    2.0000    0.4652
%     3.0000    9.0000    0.4655
%     8.0000   10.0000    0.7553
% 
% 
%  created by FB 22/11/12
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
   volume_list='';
end;

if (nargin < 3 )
   thresholds='mean+1std';
end;

if (nargin < 4)
   spacing_fact=1.5;
end;

if (nargin < 5)
   map4Score='quadr';
end;

if (nargin < 6)
   colormap=[1 0 0];
end;

if (nargin < 7)
   tmp_dir='xxx_tmp4dend_xxx';
end;


warning off;
mkdir(tmp_dir);
warning on;

tmp_name=[tmp_dir '/vol_'];
if (isempty(volume_list))
    for i=1:size(tree,1)+1
        tom_emwrite([tmp_name num2str(i) '.em'],tom_spheremask(ones(32,32,32),14));
        volume_list{i}=[tmp_name num2str(i) '.em'];
        thresholds(i)=0.5;
    end;
end;
    
if (strcmp(thresholds,'mean+1std') )
    for i=1:length(volume_list)
        tmp=tom_emreadc(volume_list{i},'resample',[4 4 4]);
        tmp=tmp.Value;
        mea=mean(tmp(:));
        stda=std(tmp(:));
        thresholds_tmp(i)=mea+(1.*stda);
    end;
    thresholds=thresholds_tmp;
end;

num_of_mod=size(tree,1)+1;


tmp_vol=tom_emread(volume_list{1});
sz=size(tmp_vol.Value);
off4dln=round(sz(1)/2);

scaleX=round(sz(1).* spacing_fact);
if (strcmp(map4Score,'quadr'))
    scaleY=(scaleX .*num_of_mod)  ./ (max(tree(:,3))-min(tree(:,3)));
end;

%call matlab function rip of dendogram to get the plotting coordinates
[~,part_ord,X,Y]=dendrogram_rip(tree,num_of_mod);

%find minimum ~0 from Y for 1. row
idx_no_zeoro=find(Y>0);
mi_no_zero=min(Y(idx_no_zeoro));
%a little hacky renorm ./scaleY but makes things better
first_row_y=((mi_no_zero.*scaleY)-scaleX)./scaleY;

%ini LUT with the first row
LUT=ini_LUT(tree,X,Y,volume_list,part_ord);

label_c=0;
zz=0;
dlc=0;
co4dendlin=[];
coord4label=[];
mod_count=num_of_mod;

%iterate over all tree entries
for i=1:size(tree,1)
    %transfer one upsidedown u to xl and yl 
    xl=X(:,i);
    yl=Y(:,i);
    
    %extract min y points (root coords)
    [x,y]=get_root_points(xl,yl); 
    
    %get root coord form lookup
    root_coord=q_LUT(LUT,x,y,'coord');
    
    %2understand: figure; plot(xl,yl); hold on;plot(root_coord(:,1),root_coord(:,2),'ro'); hold off; 
    
    %calculate coords 4 labels (last line with give root coords)
    [coord4label,label_c]=calc_label_coord(coord4label,root_coord,first_row_y,scaleX,scaleY,off4dln,label_c);
    
    %scale coords 4 chimera
    root_coord=scale_coords(root_coord,scaleX,scaleY,first_row_y);
   
    %get names from lookup table
    names=q_LUT(LUT,x,y,'name');
    
    %left low
    zz=zz+1;
    coord_list_all(:,zz)=root_coord(1,:);
    volume_list_all{zz}=names{1};
    
    %right low
    zz=zz+1;
    coord_list_all(:,zz)=root_coord(2,:);
    volume_list_all{zz}=names{2};
    
    %update LUT wit middle up ...next left low
    mod_count=mod_count+1;
    luc=size(LUT,1)+1;
    LUT{luc,1}=luc;
    LUT{luc,2}=mean(x);
    LUT{luc,3}=max(yl);
    LUT{luc,4}=[tmp_name num2str(mod_count) '.em'];
    %generate new volume;
    v1=tom_emreadc(names{1});
    v2=tom_emreadc(names{2});
    new_vol=(v1.Value+v2.Value)./2;
    tom_emwritec([tmp_name num2str(mod_count) '.em'],new_vol);
   
    %generate coord list 4 .bild file
    [co4dendlin,dlc]=build_dln_coord(root_coord(1,:),root_coord(2,:),[mean(x).*scaleX max(yl).*scaleY 0],off4dln,dlc,co4dendlin);
    
    %debug plot 2 check coordinates in matlab ...faster!
    %to switch on set last param 1
    db_plot_int(coord_list_all,co4dendlin,dlc,zz,0);
    
end;

%final middle up
zz=zz+1;
volume_list_all{zz}=LUT{end,4};
coord_list_all(:,zz)=[LUT{end,2} LUT{end,3} 0].*[scaleX scaleY 0];

%correct 4 vol center ...chimera uses top left as center (check center button)
coord_list_all=coord_list_all-off4dln;

%calc mid of volumes and dendo-tree 
mid=calc_mid(coord_list_all,co4dendlin);

%shift all coords 2 mid
[coord_list_all,co4dendlin,coord4label]=shift_all_coords(mid,coord_list_all,co4dendlin,coord4label);

%write chimera .build files 2 disk
gen_build_file(tmp_dir,co4dendlin);
gen_build_file_label(tmp_dir,coord4label,part_ord,off4dln);

all_bild{1}=[tmp_dir '/lines4dend.bild'];
all_bild{2}=[tmp_dir '/labels.bild'];

tom_chimera_volumes('show.py',volume_list_all,'',thresholds,colormap,'',coord_list_all,'',all_bild);


%*****************************************************
%Helper functions
%****************************************************
function root_coord=scale_coords(root_coord,scaleX,scaleY,first_row_y)

if (root_coord(1,2)==0)
    root_coord(1,2)=first_row_y;
end;
if (root_coord(2,2)==0)
    root_coord(2,2)=first_row_y;
end;
root_coord(:,1)=root_coord(:,1).*scaleX;
root_coord(:,2)=root_coord(:,2).*scaleY;

function [coord4label,label_c]=calc_label_coord(coord4label,root_coord,first_row_y,scaleX,scaleY,off4dln,label_c)

if (root_coord(1,2)==0)
    root_coord(1,2)=first_row_y;
    label_c=label_c+1;
    coord4label(:,label_c)=[(root_coord(1,1).*scaleX-(off4dln/3)) ((root_coord(1,2).*scaleY)-(2.5*off4dln)) 0];
end;
if (root_coord(2,2)==0)
    root_coord(2,2)=first_row_y;
    label_c=label_c+1;
    coord4label(:,label_c)=[(root_coord(2,1).*scaleX-(off4dln/3)) ((root_coord(2,2).*scaleY)-(2.5*off4dln)) 0];
end;

function LUT=ini_LUT(tree,X,Y,volume_list,part_ord)
zz=1;
for i=1:size(tree,1)
    xl=X(:,i);
    yl=Y(:,i);
    idx=find(yl==0);
    if (zz>part_ord)
        break;
    end;
    for ii=1:length(idx)
        LUT{zz,1}=zz; %running number 
        LUT{zz,2}=zz; %x coordinate just count up 
        LUT{zz,3}=yl(idx(ii)); % y-coordinate (0)
        LUT{zz,4}=volume_list{part_ord(zz)}; %particle name
        zz=zz+1;
    end;
end;

function [x,y]=get_root_points(xl,yl)

[val,idx]=sort(yl);

x=xl(idx(1:2));
y=yl(idx(1:2));

[val,idx]=sort(x);
x=x(idx);
y=y(idx);

function [coord_list_all,co4dendlin,coord4label]=shift_all_coords(mid,coord_list_all,co4dendlin,coord4label)

for i=1:size(coord_list_all,2)
    coord_list_all(:,i)=coord_list_all(:,i)-mid;
end;
for i=1:size(co4dendlin,2)
    co4dendlin(:,i)=co4dendlin(:,i)-mid;
end;
for i=1:size(coord4label,2)
    coord4label(:,i)=coord4label(:,i)-mid;
end;


function res=q_LUT(LUT,x,y,flag)

%extracts entries form LUT according 2 given x and y root coords

mat= [LUT{:,2} ;LUT{:,3}]';

for i=1:length(x)
    xM=mat(:,1)==x(i);
    yM=mat(:,2)==y(i);
    idx(i)=find(xM.*yM);
end;

for i=1:length(idx)
    if (strcmp(flag,'coord') )
        res(i,:)=[LUT{idx(i),2} LUT{idx(i),3} 0];
    end;
    if (strcmp(flag,'name') )
        res{i}=LUT{idx(i),4};
    end;
end;




function  mid=calc_mid(coord_list_all,co4dendlin)

all_coords=cat(2,coord_list_all,co4dendlin);

minX=min(all_coords(1,:));
maxX=max(all_coords(1,:));

minY=min(all_coords(2,:));
maxY=max(all_coords(2,:));

mid=[(minX+maxX)/2  (minY+maxY)/2 0]';

function gen_build_file_label(tmp_dir,coord4label,order,font_sz)

name=[tmp_dir '/labels.bild' ];

fid=fopen(name,'wt');
fprintf(fid,'.color 1 0 0\n');
zz=0;
tmp=sort(coord4label(1,:)); %as partorder is used coord just count up!
for i=1:size(coord4label,2)
    zz=zz+1;
    fprintf(fid,'.font sans %d\n',font_sz); 
    fprintf(fid,'.cmov %f %f %f\n',tmp(zz),coord4label(2,zz), coord4label(3,zz));
    fprintf(fid,'%d\n',order(i));
end;
fclose(fid);



function gen_build_file(tmp_dir,co4dendlin)
name=[tmp_dir '/lines4dend.bild' ];

fid=fopen(name,'wt');
fprintf(fid,'.color 1 0 0\n');
zz=0;
for i=1:(size(co4dendlin,2)/3)
    zz=zz+1;
    fprintf(fid,'.move %f %f %f\n',co4dendlin(1,zz),co4dendlin(2,zz), co4dendlin(3,zz));
    zz=zz+1;
    fprintf(fid,'.draw %f %f %f\n',co4dendlin(1,zz),co4dendlin(2,zz), co4dendlin(3,zz));
    zz=zz+1;
    fprintf(fid,'.draw %f %f %f\n',co4dendlin(1,zz),co4dendlin(2,zz), co4dendlin(3,zz));
end;
fclose(fid);


function [co4dendlin,dlc]=build_dln_coord(p1,p2,p3,off4dln,dlc,co4dendlin)


%swap points 2 go from left 2 right
if (p1(1)>p2(1))
    tmp_p2=p2;
    p2=p1;
    p1=tmp_p2;
end;

%**************************
%left branch
%**************************
%left bottom
dlc=dlc+1;
if ((p1(2)+off4dln) <= p3(2) )
    co4dendlin(:,dlc)=[p1(1) p1(2)+off4dln 0];
else
    co4dendlin(:,dlc)=[p1(1) p3(2) 0];
end;
%middle 
dlc=dlc+1;
co4dendlin(:,dlc)=[p1(1) p3(2) 0];

%right up
dlc=dlc+1;
if (p3(1)-off4dln >= p1(1))
    co4dendlin(:,dlc)=[p3(1)-off4dln p3(2) 0];
else
    co4dendlin(:,dlc)=[p1(1) p3(2) 0];
end;

%**************************
%right branch
%**************************
%right bottom
dlc=dlc+1;
if ((p2(2)+off4dln) <= p3(2) )
    co4dendlin(:,dlc)=[p2(1) p2(2)+off4dln 0];
else
    co4dendlin(:,dlc)=[p2(1) p3(2) 0]; 
end;
%middle
dlc=dlc+1;
co4dendlin(:,dlc)=[p2(1) p3(2) 0];
%left up
dlc=dlc+1;
if ((p3(1)+off4dln) <= p2(1) )
    co4dendlin(:,dlc)=[p3(1)+off4dln p3(2) 0];
else
    co4dendlin(:,dlc)=[p2(1) p3(2) 0];
end;

%debug plot
db_plot=0;
if (db_plot)
    figure; plot3(p1(1),p1(2),p1(3),'r+')
    hold on; plot3(p2(1),p2(2),p2(3),'b+')
    hold on; plot3(p3(1),p3(2),p3(3),'ko')
    start_db=dlc-5;
    idx=start_db:start_db+2;
    hold on;plot3(co4dendlin(1,idx),co4dendlin(2,idx),co4dendlin(3,idx),'r-');
    idx=start_db+3:start_db+5;
    hold on; plot3(co4dendlin(1,idx),co4dendlin(2,idx),co4dendlin(3,idx),'b-');
end;

function db_plot_int(coord_list_all,co4dendlin,dlc,zz,db_plot)

if (zz==2)
    figure;
end;

if (db_plot==1)
    hold on; 
    %line(xl,yl);
    plot(coord_list_all(1,(zz-1):zz),coord_list_all(2,(zz-1):zz),'ro');
    start_db=dlc-5;
    idx4pl=start_db:start_db+2;
    plot(co4dendlin(1,idx4pl),co4dendlin(2,idx4pl),'g-');
    idx4pl=start_db+3:start_db+5;
    plot(co4dendlin(1,idx4pl),co4dendlin(2,idx4pl),'g-');
    hold off;
end; 


%**************************************************
%helper function from matlab dendogram command
%**********************************************
function [T,perm,A,B] = dendrogram_rip(Z,varargin)
%DENDROGRAM Generate dendrogram plot.
%   DENDROGRAM(Z) generates a dendrogram plot of the hierarchical binary
%   cluster tree represented by Z.  Z is an (M-1)-by-3 matrix, generated by
%   the LINKAGE function, where M is the number of objects in the original
%   dataset.
%
%   A dendrogram consists of many U-shaped lines connecting objects in a
%   hierarchical tree.  The height of each U represents the distance
%   between the two objects being connected.  If there were 30 or fewer
%   data points in the original dataset, each leaf in the dendrogram
%   corresponds to one data point.  If there were more than 30 data points,
%   the complete tree can look crowded, and DENDROGRAM collapses lower
%   branches as necessary, so that some leaves in the plot correspond to
%   more than one data point.
%
%   DENDROGRAM(Z,P) generates a dendrogram with no more than P leaf nodes,
%   by collapsing lower branches of the tree.  To display the complete
%   tree, set P = 0. The Default value of P is 30.
%
%   H = DENDROGRAM(...) returns a vector of line handles.
%
%   [H,T] = DENDROGRAM(...) generates a dendrogram and returns T, a vector
%   of size M that contains the leaf node number for each object in the
%   original dataset.  T is useful when P is less than the total number of
%   objects, so some leaf nodes in the display correspond to multiple
%   objects.  For example, to find out which objects are contained in leaf
%   node k of the dendrogram, use find(T==k). When there are fewer than P
%   objects in the original data, all objects are displayed in the
%   dendrogram.  In this case, T is the identity map, i.e., T = (1:M)',
%   where each node contains only a single object.
%
%   [H,T,OUTPERM] = DENDROGRAM(...) generates a dendrogram and returns a
%   vector of the node labels of the leaves shown in the dendrogram,
%   ordered from left to right on a horizontal dendrogram and bottom to
%   top for a vertical dendrogram. OUTPERM is a permutation of the vector
%   1:P, where P is the number of nodes shown.
%   
%   [ ... ] = DENDROGRAM(..., 'PARAM1',val1, 'PARAM2',val2, ...) specifies
%   optional parameter name/value pairs:
%  
%      'Reorder'         - A numeric vector PERM representing the requested
%                          ordering of the nodes in the complete tree,
%                          ordered from left to right on a horizontal
%                          dendrogram and bottom to top for a vertical
%                          dendrogram. PERM must be a permutation of the
%                          vector 1:M. 
%
%      'CheckCrossing'   - A logical value defining whether to check if
%                          PERM will cause crossing branches in the plot.
%                          It's only useful when the 'Reorder' option is
%                          provided. Default is true. If true, a warning
%                          will be issued if PERM causes crossing branches
%                          in the plot. When the dendrogram doesn't draw
%                          the complete tree, a warning won't be given if
%                          PERM causes crossing branches in the complete
%                          tree but not in the dendrogram shown in the plot.
%
%      'ColorThreshold'  - A threshold T. Dendrogram assigns a unique color
%                          to each group of nodes within the dendrogram
%                          whose linkage is less than T. Values are:
%                            * A string 'default'. T will be set to 70% of
%                              the maximum linkage i.e. 0.7 * max(Z(:,3)).
%                            * A numerical scalar in the range of
%                              0 < T < max(Z(:,3)). If T is not given, or
%                              T is less than or equal to zero, or T is
%                              greater than the maximum linkage then the
%                              dendrogram will be drawn using only one
%                              color.
%
%      'Orientation'     - A string that orients the dendrogram within the
%                          figure window. Values are:
%                            * 'top'      -- top to bottom (default)
%                            * 'bottom'   -- bottom to top
%                            * 'left'     -- left to right
%                            * 'right'    -- right to left
%
%      'Labels'          - A character array or cell array strings S with
%                          one label for each observation. Any leaves in
%                          the tree containing a single observation are
%                          labeled with that observation's label.
%
%   Example:
%      X = rand(100,2);
%      Y = pdist(X,'cityblock');
%      Z = linkage(Y,'average');
%      [H, T] = dendrogram(Z);
%   
%      rng('default')
%      X=rand(10,3);
%      Z=linkage(X,'ave');
%      % Draw the dendrogram using the default setting
%      subplot(2,1,1); dendrogram(Z);
%      % Reorder the nodes shown in the plot. 
%      subplot(2,1,2); dendrogram(Z,'reorder',[ 9 5 8 10 2 4 1 6 7 3 ]);
%
%   See also LINKAGE, PDIST, CLUSTER, CLUSTERDATA, COPHENET, INCONSISTENT,
%   SILHOUETTE.

%   Copyright 1993-2011 The MathWorks, Inc.
%   $Revision: 1.1.8.7 $

numLeaves = size(Z,1)+1; %the number of observations

if nargin < 2
    p = 30;
end

if nargin == 2
    p = varargin{1};
end

color = false;
orientation = 't'; %default top
obslabels = [];
threshold = 0.7 * max(Z(:,3));
leafOrder = [];

if nargin > 2
    if isnumeric(varargin{1})
        p = varargin{1};
        offset = 1;
          
    else
        p = 30;
        offset = 0;
    end
    
    pnames = {'orientation' 'colorthreshold' 'labels'  'reorder','checkcrossing'};
    dflts =  {orientation   'default'        obslabels leafOrder, true};
    [orientation,threshold,obslabels,leafOrder,check,setFlag] =  ...
        internal.stats.parseArgs(pnames, dflts, varargin{1+offset:end});
    
    if ~isscalar (check) || ~islogical (check)
         error(message('stats:dendrogram:BadCheck'));
    end
    
    if ~isempty(orientation) && ischar(orientation)
        orientation = lower(orientation(1));
    else
        orientation = 0;    % bad value
    end
    
    if  orientation == 0 || ~ismember(orientation,{'t','b','r','l'}) 
      warning(message('stats:dendrogram:BadOrientation'));
    end
   
    
    color = setFlag.colorthreshold; % threshold argument was given
    if color
        if ischar(threshold) && strncmpi(threshold,'default',length(threshold))
            threshold = 0.7 * max(Z(:,3));
        elseif ~isnumeric(threshold)
            error(message('stats:dendrogram:BadThreshold'));
        end
    end
    
    if ~isempty(obslabels)
        if ischar(obslabels)
            obslabels = cellstr(obslabels);
        elseif ~iscellstr(obslabels)
            error(message('stats:dendrogram:BadLabels'));
        end
        if ~isvector(obslabels) || numel(obslabels) ~= numLeaves
            error(message('stats:dendrogram:InputSizeMismatch'));
        end
        obslabels = obslabels(:);
    end
    
    if ~isempty(leafOrder)
        if (~isvector(leafOrder) || numel(leafOrder)~=numLeaves )
            error(message('stats:dendrogram:BadLeafOrder'));
        else
            leafOrder = leafOrder(:)'; %make it to be a row vector
            if ~isequal(sort(leafOrder),1:numLeaves)
                error(message('stats:dendrogram:BadLeafOrder'));
            end
        end
    end
end

if ~isscalar(p) || p < 0 || p == 1
            error(message('stats:dendrogram:BadP'));
end    

% For each node currently labeled numLeaves+k, replace its index by
% min(i,j) where i and j are the nodes under node numLeaves+k.
Z = transz(Z);
T = (1:numLeaves)'; 

% If there are more than p nodes, the dendrogram looks crowded.
% The following code will make the last p link nodes into leaf nodes,
% and only these p nodes will be visible.
if (numLeaves > p) && (p ~= 0)
    
    Y = Z((numLeaves-p+1):end,:);         % get the last nodes
    
    R = unique(Y(:,1:2));
    Rlp = R(R<=p);
    Rgp = R(R>p);
    W(Rlp) = Rlp;                 % use current node number if <=p
    W(Rgp) = setdiff(1:p, Rlp);   % otherwise get unused numbers <=p
    W = W(:);
    T(R) = W(R);
    
    % Assign each leaf in the original tree to one of the new node numbers
    for i = 1:p
        c = R(i);
        T = clusternum(Z,T,W(c),c,numLeaves-p+1,0); % assign to its leaves.
    end
    
    % Create new, smaller tree Z with new node numbering
    Y(:,1) = W(Y(:,1));
    Y(:,2) = W(Y(:,2));
    % At this point, it's possible that Z(i,1) < Z(i,i) for some rows.
    %The newly formed cluster will always be represented by the number in Z(i,1);
    Z = Y; 
    
    %Assign each leaf in the leafOrder to one of the new node numbers.
    if ~isempty(leafOrder)
        leafOrder = T(leafOrder)';
        d = diff(leafOrder);
        d = [1 d];
        leafOrder = leafOrder (d~=0);
        if numel(leafOrder) ~= p
            error(message('stats:dendrogram:InvalidLeafOrder'));
        end       
    end
    
    numLeaves = p; % reset the number of nodes to be p (row number = p-1).
end

A = zeros(4,numLeaves-1);
B = A;
X = 1:numLeaves; %the initial points for observation 1:n
Y = zeros(numLeaves,1);

if isempty(leafOrder)
    r = Y;
    % arrange Z into W so that there will be no crossing in the dendrogram.
    W = zeros(size(Z));
    W(1,:) = Z(1,:);
    
    nsw = zeros(numLeaves,1); rsw = nsw;
    nsw(Z(1,1:2)) = 1; rsw(1) = 1;
    k = 2; s = 2;
    
    while (k < numLeaves)
        i = s;
        while rsw(i) || ~any(nsw(Z(i,1:2)))
            if rsw(i) && i == s
                s = s+1;
            end
            i = i+1;
        end
        
        W(k,:) = Z(i,:);
        nsw(Z(i,1:2)) = 1;
        rsw(i) = 1;
        if s == i
            s = s+1;
        end
        k = k+1;
    end
    
    % initialize X based on W
    g = 1;
    for k = 1:numLeaves-1
        i = W(k,1); %the left node in W(k,:)
        if ~r(i),
            X(i) = g;
            g = g+1;
            r(i) = 1;
        end
        i = W(k,2); %the right node in W(k,:)
        if ~r(i),
            X(i) = g;
            g = g+1;
            r(i) = 1;
        end
    end
    perm(X) = 1:numLeaves;
    
else % ~isempty(leafOrder) if a leaf order is specified use
    %get X based on the specified order
    X(leafOrder) = 1:numLeaves;
    %need to check whether leafOrder will have crossing branch
    if (check)
       checkCrossing(Z(:,[1 2]), leafOrder);
    end
    perm = leafOrder;
end


label = num2str(perm');
if ~isempty(obslabels)
    label = cellstr(label);
    % label(:) = {''};   % to make sure non-singletons get an empty label
    singletons = find(histc(T,1:numLeaves)==1);
    for j=1:length(singletons)
        sj = singletons(j);
        label(perm==sj) = obslabels(T==sj);
    end
end
% set up the color

theGroups = 1;
groups = 0;
cmap = [0 0 1];

if color
    groups = sum(Z(:,3)< threshold);
    if groups > 1 && groups < (numLeaves-1)
        theGroups = zeros(numLeaves-1,1);
        numColors = 0;
        for count = groups:-1:1
            if (theGroups(count) == 0)
                P = zeros(numLeaves-1,1);
                P(count) = 1;
                P = colorcluster(Z,P,Z(count,1),count);
                P = colorcluster(Z,P,Z(count,2),count);
                numColors = numColors + 1;
                theGroups(logical(P)) = numColors;
            end
        end
        cmap = hsv(numColors);
        cmap(end+1,:) = [0 0 0];
    else
        groups = 1;
    end
    
end


col = zeros(numLeaves-1,3);

for n = 1:(numLeaves-1)
    i = Z(n,1); j = Z(n,2); w = Z(n,3);
    A(:,n) = [X(i) X(i) X(j) X(j)]';
    B(:,n) = [Y(i) w w Y(j)]';
    X(i) = (X(i)+X(j))/2; Y(i)  = w;
    if n <= groups
        col(n,:) = cmap(theGroups(n),:);
    else
        col(n,:) = cmap(end,:);
    end
end


% ---------------------------------------
function T = clusternum(X, T, c, k, m, d)
% assign leaves under cluster c to c.

d = d+1;
n = m; flag = 0;
while n > 1
    n = n-1;
    if X(n,1) == k % node k is not a leave, it has subtrees
        T = clusternum(X, T, c, k, n,d); % trace back left subtree
        T = clusternum(X, T, c, X(n,2), n,d);
        flag = 1; break;
    end
end

if flag == 0 && d ~= 1 % row m is leaf node.
    T(X(m,1)) = c;
    T(X(m,2)) = c;
end
% ---------------------------------------
function T = colorcluster(X, T, k, m)
% find local clustering

n = m;
while n > 1
    n = n-1;
    if X(n,1) == k % node k is not a leave, it has subtrees
        T = colorcluster(X, T, k, n); % trace back left subtree
        T = colorcluster(X, T, X(n,2), n);
        break;
    end
end
T(m) = 1;
% ---------------------------------------
function Z = transz(Z)
%TRANSZ Translate output of LINKAGE into another format.
%   This is a helper function used by DENDROGRAM and COPHENET.

%   In LINKAGE, when a new cluster is formed from cluster i & j, it is
%   easier for the latter computation to name the newly formed cluster
%   min(i,j). However, this definition makes it hard to understand
%   the linkage information. We choose to give the newly formed
%   cluster a cluster index M+k, where M is the number of original
%   observation, and k means that this new cluster is the kth cluster
%   to be formed. This helper function converts the M+k indexing into
%   min(i,j) indexing.

numLeaves = size(Z,1)+1;

for i = 1:(numLeaves-1)
    if Z(i,1) > numLeaves
        Z(i,1) = traceback(Z,Z(i,1));
    end
    if Z(i,2) > numLeaves
        Z(i,2) = traceback(Z,Z(i,2));
    end
    if Z(i,1) > Z(i,2)
        Z(i,1:2) = Z(i,[2 1]);
    end
end


function a = traceback(Z,b)

numLeaves = size(Z,1)+1;

if Z(b-numLeaves,1) > numLeaves
    a = traceback(Z,Z(b-numLeaves,1));
else
    a = Z(b-numLeaves,1);
end
if Z(b-numLeaves,2) > numLeaves
    c = traceback(Z,Z(b-numLeaves,2));
else
    c = Z(b-numLeaves,2);
end

a = min(a,c);


function checkCrossing(Z, order)
%check whether the give Tree will have crossing branches with the given
%permutation vector

numBranches = size(Z,1);
numLeaves = numBranches + 1;
%numLabels = numBranches + numLeaves;
% reorder the tree
perm = order(:);
%XPos is the position indices for leaves 1:numLeaves. 
XPos(perm) = 1:numLeaves;
Z0 = Z; %keep the original tree
%renumber the leave nodes in Z such that number N represents the Nth nodes
%in the plot
Z = XPos(Z);

% Check if the reordered tree structure leads to a tree with no crossing branches
minPos = 1:numLeaves;
maxPos = 1:numLeaves;
sz = ones(numLeaves,1);

for i = 1:numBranches
    currentMinPos = min(minPos(Z(i,:)));
    currentMaxPos = max(maxPos(Z(i,:)));
    currentSize = sum(sz(Z(i,:)));
    if currentMaxPos - currentMinPos + 1 ~= currentSize 
        warning(message('stats:dendrogram:CrossingBranches'));
        break;
    end
    
    j =XPos(Z0(i,1));% j is the cluster number for the newly formed cluster.
    %Note that we can't use j = XPos(min(Z0(i,:))), Because when not all of
    % the points are shown, the value in the first column of Z may be bigger
    % than the value in the second column.
    minPos(j)= currentMinPos;
    maxPos(j) = currentMaxPos;
    sz(j) = currentSize;
    
end





