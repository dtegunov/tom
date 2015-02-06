function inp_ref_stack=tom_av2_index_inp_sample_refstack(ref_stack,increment,mask,norm_flag)
%TOM_AV2_INDEX_INP_SAMPLE_REFSTACK performs in-plane rotation of given
%reference stack
%
%   inp_ref_stack=tom_av2_index_inp_sample_refstack(ref_stack,increment,mask,norm_flag)
%PARAMETERS
%
%  INPUT
%   ref_stack            stack of referecnes
%   increment            increment for rotation
%   mask                 mask for rotaion
%   mask_cc              mask for cross correlation function ...to get rid of side peaks
%   norm_flag            flag for norming the refereces 
%   
%  OUTPUT
%   inp_ref_stack       in plane rotation sampled stack
%       
%   
%
%EXAMPLE
%   inp_ref_stack=tom_av2_index_inp_sample_refstack(rf.Value,5);
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


sz_rs=size(ref_stack);
if (length(size(ref_stack))==2)
    sz_rs(3)=1;
end;

if nargin < 2
    increment=5;
end;

if nargin < 3
    mask=tom_spheremask(ones(sz_rs(1),sz_rs(2)),round(sz_rs(1)./2)-2,2);
end

if nargin < 4
    norm_flag='phase';
end



rotations=[0:increment:(360-increment)];


%allocate some memory
inp_ref_stack=zeros(sz_rs(1),sz_rs(2),length(rotations).*sz_rs(3));

%loop over all references
zz=1;
for i=1:size(ref_stack,3)
    %loop over all in-plane rotions
    for ii=1:length(rotations)
        tmp=tom_norm((ref_stack(:,:,i)+1).*2,norm_flag).*mask;
        inp_ref_stack(:,:,zz)=tom_rotate(tmp,rotations(ii));
        zz=zz+1;
    end;
end;


