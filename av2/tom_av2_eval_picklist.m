function [pl_false_pos pl_false_neg]=tom_av2_eval_picklist(pick_hand,pick_auto)
%TOM_AV2_EVAL_PICKLIST compares 2 picklists
%
%   [pl_false_pos pl_false_neg]=tom_av2_eval_picklist(pick_hand,pick)
%
%PARAMETERS
%
%  INPUT
%   pick_hand           hand picklist 
%   pick_auto           auto picklist
%
%
%  OUTPUT
%    pl_false_pos       pl false negative
%    pl_false_neg       pl fals positive
% 
%EXAMPLE
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


%parse inputs


lookup=transform_pl(pick_hand,pick_auto);



pl_false_pos=align2d;
pl_false_neg=align2d;


for i=1:length(lookup.filename)

    %all_pos=
    
    pos_auto=lookup.pos_and_id_auto{i};
    pos_hand=lookup.pos_and_id_hand{i};
    for ii=1:size(pos_auto,2)
        [pointidx, pointcoords, distance] = tom_nearestpoint(pos_auto(i,:),pos_hand(:,1:2));
        if (distance >  round(pick_hand(1,1).radius./2) )
            pl_false_pos(i)=pick_auto(i);
            zz_fp=zz_fp+1;
        end;
        if (distance <  round(pick_hand(1,1).radius./2) )
            
        end;
            
    end;
  %  [pointidx, pointcoords, distance] = tom_nearestpoint(pos_auto,);
    
    
end;




function lookup=transform_pl(align_low,align_high)


%new_idx=zeros(size(align2d,2),1);


%get unique list of filenames
for i=1:size(align_low,2)
    all_names_low{i}=align_low(1,i).filename;
end;
name_list_low=unique(all_names_low);

for i=1:length(name_list_low)
    [a b c]=fileparts(name_list_low{i});
    name_list_low{i}=[b c];
end;

clear('all_names_low');

%get unique list of filenames
for i=1:size(align_high,2)
    all_names_high{i}=align_high(1,i).filename;
end;
name_list_high=unique(all_names_high);
for i=1:length(name_list_high)
    [a b c]=fileparts(name_list_high{i});
    name_list_high{i}=[b c];
end;

clear('all_names_high');


look_count=1;
%build up lookuptable
five_p=round(length(name_list_low)./20);
zz_p=1;
tic;


for i=1:length(name_list_low)
    
    
    idx=find(ismember(name_list_high,name_list_low{i}));
    if (isempty(idx))
        continue;
    end;
    
    
    pos_low=get_pos_per_img(align_low,name_list_low{i});
    pos_high=get_pos_per_img(align_high,name_list_low{i});
    
    
    
   
    
    %save in lookup table
    lookup.filename{look_count}=name_list_low{i};
    lookup.pos_and_id_hand{look_count}=pos_low;
    lookup.pos_and_id_auto{look_count}=pos_low;
    look_count=look_count+1;
    
    if (mod(i,five_p)==0)
        toc;
        disp([num2str(zz_p.*5) '% done'] );
        zz_p=zz_p+1;
        tic;
    end;
    
end;

function pos=get_pos_per_img(align,img_name)


pos=zeros(40000,3);
zz=1;
for i=1:size(align,2)
    if (strfind(align(1,i).filename,img_name))
        pos(zz,:)=[align(1,i).position.x align(1,i).position.y i];
        zz=zz+1;
    end;
end;

if (zz==1)
    pos=-1;
    return;
end;

tmp=pos(1:zz-1,:);
clear('pos');

pos=tmp;



