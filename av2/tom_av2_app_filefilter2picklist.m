function  filt_picklist=tom_av2_app_filefilter2picklist(align2d,filefilter,f_pickList_out,app_img,app_poly,verbose)
%tom_av2_app_filefilter2picklist filters a picklist according 2 a given filefilter
%                                ...removes particles inside polygones and on bad images
%
%   f_picklist=tom_av2_app_filefilter2picklist(pickList,filefilter,f_pickList_out,app_img,app_poly)
%
%PARAMETERS
%
%  INPUT
%   align2d             filename or var  picklist
%   filefilter          filename or var  filefilter
%   f_pickList_out      (opt.) filename of output picklist
%   app_img             (1) apply good bad of filefilter 
%   app_poly            (1) apply poly of filefilter
%   verbose             (1) verbose flag    
%  
%  OUTPUT
%   filt_picklist       filtered picklist
%
%
%EXAMPLE
%   
%   tom_av2_app_filefilter2picklist('pick_raw.mat','filefilter_poly.mat');
%
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by SN/FB 01/24/06
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

if (nargin < 3)
    f_pickList_out='';
end;

if (nargin < 4)
    app_img=1;
end;

if (nargin < 5)
    app_poly=1;
end;

if (nargin < 6)
    verbose=1;
end;


if (isstruct(align2d)==0)
    load(align2d)
end;

if (isstruct(filefilter)==0)
    load(filefilter)
    filefilter=particlepicker;
end;

for i=1:100
    try
        h=tom_reademheader(align2d(1,i).filename);
        sz_img=h.Header.Size(1:2)';
        break;
    catch Me
    end;
end;

%build index
for i=1:size(align2d,2)
    [a b c]=fileparts(align2d(1,i).filename);
    all_filnames_al{i}=[b c];
end;

sum_ff=0;
sum_pl=0;

picklist_new=[];
for i=1:length(filefilter.filelist)
    if (app_img==0 || filefilter.filefilter{i}==1)
        ind=find(ismember(all_filnames_al,filefilter.filelist{i}));
        if (isempty(ind)==0)
            if (app_poly==1 && i <= length(particlepicker.maskcell) )
                pick_outside=outside_poly_mask(align2d(1,ind),particlepicker.maskcell{i},sz_img);
            else
                pick_outside=align2d(1,ind);
            end;
            picklist_new=cat(2,picklist_new,pick_outside);
            if (size(pick_outside,2)~= size(align2d(1,ind),2))
                if (verbose~=0)
                    disp([filefilter.filelist{i} ':(ff_idx=' num2str(i) ') ' num2str(size(align2d(1,ind),2)) ' to ' num2str(size(pick_outside,2)) ' by polygon']);
                end;
            end;
            sum_pl=sum_pl+((size(align2d(1,ind),2)) - size(pick_outside,2));
        end;
       
     else
            tmpp_ind=find(ismember(all_filnames_al,filefilter.filelist{i}));
            if (verbose~=0 && isempty(tmpp_ind)==0)
                sum_ff=sum_ff+length(tmpp_ind);
                disp([filefilter.filelist{i} ':(ff_idx=' num2str(i)  ') '   num2str(length(tmpp_ind)) ' to 0 by filefilter']);
            end;
     end;
end;

if (verbose~=0)
    disp(' ');
    disp('Summary:');
    disp([num2str(sum_ff) ' parts removed by filefilter']);
    disp([num2str(sum_pl) ' parts removed by polygon']);
end;
    
if (isempty(f_pickList_out)==0)
    align2d=picklist_new;
    save(f_pickList_out,'align2d');
    clear('align2d');
end;

filt_picklist=picklist_new;

function pick_ouside=outside_poly_mask(pick_img,polis,sz_img)

if (isempty(polis))
    pick_ouside=pick_img;
else
    inside=zeros(size(pick_img,2),1);
    for i=1:size(pick_img,2)
        pos_tmp=[floor(pick_img(1,i).position.x) floor(pick_img(1,i).position.y)];
        for ii=1:length(polis)
            gg=polis{ii};
            in_tmp=inpolygon(pos_tmp(1),pos_tmp(2),gg(:,1),gg(:,2));
            if (in_tmp)
                inside(i)=1;
            end;
        end;
    end;
    pick_ouside=pick_img(1,find(inside==0));
end;

