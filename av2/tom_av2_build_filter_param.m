function filter_param=tom_av2_build_filter_param(filter_param,sz,default_flag,application)
%TOM_AV2_BUILD_FILTER_PARAM creates ...
%
%   filter_param=tom_av2_build_filter_param(filter_param,sz,default_flag,application)
%
%PARAMETERS
%
%  INPUT
%   filter_param        ...
%   sz                  ...
%   default_flag        ...
%   application         ...
%  
%  OUTPUT
%   filter_param		...
%
%EXAMPLE
%   ... = tom_av2_build_filter_param(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by .. mm/dd/yy
%   updated by ..
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

if (strcmp(default_flag,'default'))
    default=2;
else
    default=0;
end;


if (strcmp(application,'tom_av2_multi_ref_alignment'))
    
    
    if (isstruct(filter_param)==0 )
        filter_param.mask.classify1.Apply=default; filter_param.mask.classify1.Value(1)=sz(1); filter_param.mask.classify1.Value(2)=sz(2);
        filter_param.mask.classify2.Apply=default; filter_param.mask.classify2.Value(1)=sz(1); filter_param.mask.classify2.Value(2)=sz(2);
        filter_param.mask.align.Apply=default;  filter_param.mask.align.Value(1)=sz(1); filter_param.mask.align.Value(2)=sz(2);
        filter_param.mask.ccf_rot.Apply=default; szz=size(tom_cart2polar(ones(sz(1)))); filter_param.mask.ccf_rot.Value(1)=szz(1); filter_param.mask.ccf_rot.Value(2)=szz(2);
        filter_param.mask.ccf_trans.Apply=default; filter_param.mask.ccf_trans.Value(1)=sz(1); filter_param.mask.ccf_trans.Value(2)=sz(2);
        filter_param.filter.align.Apply=default;
        filter_param.filter.classify.Apply=default;
    end;
    
    if (isfield(filter_param,'mask')==0)
        filter_param.mask.classify1.Apply=default; filter_param.mask.classify1.Value(1)=sz(1); filter_param.mask.classify1.Value(2)=sz(2);
        filter_param.mask.classify2.Apply=default; filter_param.mask.classify2.Value(1)=sz(1); filter_param.mask.classify2.Value(2)=sz(2);
        filter_param.mask.align.Apply=default;  filter_param.mask.align.Value(1)=sz(1); filter_param.mask.align.Value(2)=sz(2);
        filter_param.mask.ccf_rot.Apply=default; sz=size(tom_cart2polar(ones(sz(1)))); filter_param.mask.ccf_rot.Value(1)=sz(1); filter_param.mask.ccf_rot.Value(2)=sz(2);
        filter_param.mask.ccf_trans.Apply=default; filter_param.mask.ccf_trans.Value(1)=sz(1); filter_param.mask.ccf_trans.Value(2)=sz(2);
    end;
   
    if (isfield(filter_param,'filter')==0)
        filter_param.filter.align.Apply=default;
        filter_param.filter.classify.Apply=default;
    end;
    
    
    if (isfield(filter_param.mask,'classify1')==0)
         filter_param.mask.classify1.Apply=default; filter_param.mask.classify1.Value(1)=sz(1); filter_param.mask.classify1.Value(2)=sz(2);
    end;
    
    if (isfield(filter_param.mask,'classify2')==0)
         filter_param.mask.classify2.Apply=default; filter_param.mask.classify2.Value(1)=sz(1); filter_param.mask.classify2.Value(2)=sz(2);
    end;
    
    if (isfield(filter_param.mask,'align')==0)
         filter_param.mask.align.Apply=default;  filter_param.mask.align.Value(1)=sz(1); filter_param.mask.align.Value(2)=sz(2);
    end;
    
    if (isfield(filter_param.mask,'ccf_rot')==0)
         filter_param.mask.ccf_rot.Apply=default; szz=size(tom_cart2polar(ones(sz(1)))); filter_param.mask.ccf_rot.Value(1)=szz(1); filter_param.mask.ccf_rot.Value(2)=szz(2);
    end;
    
    if (isfield(filter_param.mask,'ccf_trans')==0)
                filter_param.mask.ccf_trans.Apply=default; filter_param.mask.ccf_trans.Value(1)=sz(1); filter_param.mask.ccf_trans.Value(2)=sz(2);
    end;
    
    if (isfield(filter_param.filter,'align')==0)
        filter_param.filter.align.Apply=default;
    end;
    
    if (isfield(filter_param.filter,'classify')==0)
        filter_param.filter.classify.Apply=default;
    end;
    
end;
