function vol_sym=tom_symref_part(vol,st,vol_out,verbose,debug_drop)
%TOM_SYMREF_PART does n-fold symmetrization of a 3D reference
%
%    vol_sym=tom_symref_part(vol,st,vol_out)
%
%   
%   differnt symmetries for different masks can be applied
%
%PARAMETERS
%
%  INPUT
%   vol                 3D volume to be symmterized
%   st                  matlab struct containig masks and syms
%   vol_out             output volume name
%   verbose             (1) output control
%   debub_drop          (0) volume is written before sym with .org 
%  
%  OUTPUT
%   
%
%EXAMPLE
%   
% mod=tom_spheremask(ones(128,128,128),5,0,[80 80 71]) + tom_spheremask(ones(128,128,128),5,0,[93 93 59]);
%  
% tmp=zeros(128,128,128);
% tmp(:,:,65:end)=1;
% st.mask{1}=tmp;
%
% tmp=zeros(128,128,128);
% tmp(:,:,1:64)=1;
% st.mask{2}=tmp;
%
% st.sym=[8 21];
% 
% vol_sym=tom_symref_part(mod,st);
%
%   
%NOTE:
%
% st.mask cell containing masks which are applied before sym
% st.sym  array containin symmetries
%
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by FF 07/28/03
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


if (nargin<4)
    verbose=1;
end;

if (nargin<5)
    debug_drop=0;
end;



if (ischar(st))
    load(st);
end;

if (ischar(vol))
    if (tom_isemfile(vol))
        vol_in=vol;
        vol=tom_emread(vol);
        em_flag=1;
        sp_flag=0;
    end;
    if (tom_isspiderfile(vol))
        if (exist(strrep(vol,'reconstruction','reconstructionLess'),'file') )
           vol2=tom_spiderread(strrep(vol,'reconstruction','reconstructionLess'));             
           disp('reconstructionLess read');
        end;
        vol_in=vol;
        vol=tom_spiderread(vol);
       
        sp_flag=1;
        em_flag=0;
    end;
    vol=vol.Value;
else
    vol_in='in-memory';
end;


if (verbose==1)
    disp(' ');
end;
    
tmp=zeros(size(st.mask{1}));
for i=1:length(st.mask)
    if (verbose==1)
        disp([vol_in ' * st.mask{' num2str(i) '} ==> sym: ' num2str(st.sym(i))  ' + ' ]);
    end;
    if (st.sym(i)>1)   
        if (i>2 && exist('vol2','var') )
            tmp=tmp+tom_symref(vol2.Value.*st.mask{i},st.sym(i));
            disp('rec less used');
        else
            tmp=tmp+tom_symref(vol.*st.mask{i},st.sym(i));
        end;
    else
        tmp=tmp+(vol.*st.mask{i});
    end;
end;

if (isfield(st,'mask_all'))
    if (verbose==1)
        disp('background permutation with st.mask_all');
    end;
    tmp=tom_permute_bg(tmp,st.mask_all);
end;
    

if (nargin<3)
    vol_sym=tmp;
    return;
end;

if (sp_flag==1)
    tom_spiderwrite(vol_out,tmp);
    if (debug_drop)
         tom_spiderwrite([vol_out '.org'],tmp);
    end;
end;

if (em_flag==1)
    tom_emwrite(vol_out,tmp);
    if (debug_drop)
         tom_emwrite([vol_out '.org'],tmp);
    end;
end;

if (verbose==1)
    disp(['==> ' vol_out ' written !' ]);
    disp(' ');
end;



