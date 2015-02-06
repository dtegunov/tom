function [nr,pos_out,val_out,angle_out,shift_out,ratio_vect]=tom_av2_index_search(index_stack,lookup,part,align_flag,mask,mask_cc,demo_flag)
%TOM_AV2_INDEX_SEARCH performs multiref alilgnment with a index search
%
%   [idx,c,sumd]=tom_pca_and_k_means(stack,num_of_eigs,num_of_classes,binning,equal_flag,demo)
%PARAMETERS
%
%  INPUT
%   index_stack          index stack computed of the ref stack by tom_av2_index_calc
%   lookup               lookup table for the org pos in ref stack by tom_av2_index_calc 
%   part                 particle which should be aligned to ref stack
%   align_flag           (no alignment) flag for in-plane alignment ...use 'pre alignment' for in plane alignment 
%   mask_cc              mask for cross correlation function ...to get rid of side peaks
%   demo_flag            (no_demo) use 'demo' or 'demo_final' 
%   
%  OUTPUT
%   nr                  number of best matching reference
%   pos_out             position of cc-peak      
%   val_out             value of cc-peak
%   angle_out           output angle only computet if align_flag is set to 'pre alignment'
%                       angle_out is set 0 if align_flag is set to 'no alignment'   
% 
%
%EXAMPLE
%   [ref_num peak val]=tom_av2_index_search(idx_stack.Value,lookup,tom_bandpass(tmp_im,3,70),'no_alignment',mask_cc);
%   
%
%REFERENCES
%
%SEE ALSO
%   tom_av2_index_calc, tom_av2_index_bintree_not_2_index.m,
%   tom_av2_index_plot_indexstack.m
%
%   created by fb (eckster)
%   updated by ...
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




if nargin < 4
    align_flag='no alignment';
end;

if nargin < 5
    mask=ones(size(part));
end;

if nargin < 6
    mask_cc=ones(size(part));
end;

if nargin < 7
    demo_flag='no_demo';
end;

level=[];
sz_ind_stack=size(index_stack);

%demo_flag='demo';

if (strcmp(demo_flag,'demo'))
    figure;
end;



angle_out=0;
rott=[0 0];
shift_out=[0 0];
middle_part=[round(size(part,1)./2) round(size(part,1)./2)]+1;

for i=1:100000
    
    %names for left and right
    left_name=strrep([num2str(level) num2str(0) ],' ','');
    right_name=strrep([num2str(level) num2str(1) ],' ','');
    
   
    %check for end of tree 
    sind_left=tom_av2_index_bintree_not_2_index(left_name);
    sind_right=tom_av2_index_bintree_not_2_index(right_name);
    
    if (sind_left > sz_ind_stack(3) || sind_right > sz_ind_stack(3))
        break;
    end;
    
    if (std2(index_stack(:,:,sind_left))==0  || std2(index_stack(:,:,sind_right))==0)
        break;
    end;
    
    
    idx_left=index_stack(:,:,sind_left);
    idx_right=index_stack(:,:,sind_right);
    
    
    %left branch
    if (strcmp(align_flag,'pre Alignment'))
        [t_angle_out t_shift_out cc part_alg_left]=tom_av2_align(idx_left,part,mask,tom_cart2polar(ones(size(part))),mask_cc,'',5,0);
        rott(1)=t_angle_out(1);
        shiftt(1,:)=t_shift_out;
        %[a b]=tom_peak(cc.*mask_cc);
        b=cc;
        a=shiftt(1:2) + middle_part;
    else
        cc=tom_corr(part,idx_left,'norm');
        [a b]=tom_peak(cc.*mask_cc);
        shiftt(1,:)=a(1,:)-middle_part;
    end;
    vals(1)=b; pos(1,:)=a(1,:);
    
    %right branch
    if (strcmp(align_flag,'pre Alignment'))
        [t_angle_out t_shift_out cc part_alg_right]=tom_av2_align(idx_right,part,mask,tom_cart2polar(ones(size(part))),mask_cc,'',5,0);
        rott(2)=t_angle_out(1);
        shiftt(2,:)=t_shift_out;
        %[a b]=tom_peak(cc.*mask_cc);
        b=cc;
        a=shiftt(1:2) + middle_part;
    else
        cc=tom_corr(part,idx_right,'norm');
        [a b]=tom_peak(cc.*mask_cc);
        shiftt(2,:)=a(1,:)-middle_part;
    end;
    vals(2)=b;  pos(2,:)=a(1,:);
    
    %get max
    [a b]=max(vals);
    
    level=[level num2str(b-1)];
   
    disp(['level:' num2str(i)  '  ' left_name ' ' right_name  '  cc-ratio%: '  num2str(abs(round((vals(1)./vals(2)-1).*10000))/100)  '  cc1:' num2str(vals(1))  ' cc2:' num2str(vals(2)) ]);
    ratio_vect(i)=abs(round((vals(1)./vals(2)-1).*10000))/100;
    if (strcmp(demo_flag,'demo'))
        disp(['level:' num2str(i)  '  ' left_name ' ' right_name  '  cc-ratio%: '  num2str(abs(round((vals(1)./vals(2)-1).*10000))/100 ) ]);
        subplot(2,3,1); tom_imagesc(idx_left); title('left');
        subplot(2,3,2); tom_imagesc(idx_right); title('right');
        subplot(2,3,3); plot(vals); title('ccs');
        subplot(2,3,5); tom_imagesc(part_alg_left); title('particle aligned left');
        subplot(2,3,6); tom_imagesc(part_alg_right); title('particle aligned right');
        disp('');
    end;
    
    
    
end;

%find position in original stack
for i=1:size(lookup,1);
     tmp(i)=strcmp(lookup{i,1},level);  
end;

%transfer return values
nr=lookup{find(tmp),2};
pos_out=pos(b,:);
val_out=vals(b);
angle_out=rott(b);
shift_out=shiftt(b,:);

if (strcmp(demo_flag,'demo_final'))
    figure;
    subplot(2,3,1  ); tom_imagesc(idx_left); title('left');
    subplot(2,3,2); tom_imagesc(idx_right); title('right');
    subplot(2,3,3); plot(vals); title('ccs');
    subplot(2,3,5); tom_imagesc(part); title('particle');
end;
