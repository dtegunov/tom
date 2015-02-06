function filter_param=tom_av2_build_filter_param(filter_param,sz,default_flag,application)

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
