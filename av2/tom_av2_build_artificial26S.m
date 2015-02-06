function vol=tom_av2_build_artificial26S(flag);
%TOM_AV2_BUILD_ARTIFICIAL26S creates ...
%
%   vol=tom_av2_build_artificial26S(flag);
%
%PARAMETERS
%
%  INPUT
%   flag                ...
%  
%  OUTPUT
%   vol                 ...
%
%EXAMPLE
%   ... = tom_av2_build_artificial26S(...);
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

if (nargin==0)
    flag='2caps';
end;

%load core 
core=tom_emread('/fs/bmsan/pub/4Florian/20S_core_4.1766A.em');

%rescale core
fact=4.176./3.6;
new_size=round(size(core.Value).*fact);
core=tom_rescale3d(core.Value,new_size);
core=tom_rotate(core,[270 90 90],'linear');
core=tom_filter(core,12,'quadr','real');
[thresh,mass,vol_back]=tom_calc_isosurface(-tom_norm(core,1),720,4.2,0.1);
%core=((core < 0).*core);
%core=tom_filter(core,3,'quadr','real');
core=((-tom_norm(core,1))>thresh);

core=-tom_norm(core,1);

%load 26S
%whole_part=tom_emread('/fs/sally01/lv01/pool/pool-baumeister/data/all_7mue/log/rec2_sym/data/step3/model/model_67');
whole_part=tom_emread('/fs/sally/pool-baumeister/data/all_7mue/log/rec2_sym/data/step3/model/model_67');
whole_part=tom_filter(whole_part.Value,4,'quadr','real').*1;

[thresh,mass,vol_back]=tom_calc_isosurface(-tom_norm(whole_part,1),2700,3.6,0.01);

whole_part=((-tom_norm(whole_part,1))>thresh);

%cut out caps
cap1=tom_cut_out(whole_part,[100 60 50],[50 60 60],'no-fill');
cap2=tom_cut_out(whole_part,[12 50 60],[50 60 60],'no-fill');
%cap1=tom_filter(cap1,4,'quadr','real').*1;
%cap2=tom_filter(cap2,4,'quadr','real').*1;

cap1=-tom_norm(cap1,1);
cap2=-tom_norm(cap2,1);



%build new 26S
vol=zeros(160,160,160);

vol=tom_paste(vol,core,[26 26 26]);

if (strcmp(flag,'2caps')==1 | strcmp(flag,'1cap')==1)
    vol=tom_paste(vol,cap2,[13 47 57]);
end;

if (strcmp(flag,'2caps')==1)
    vol=tom_paste(vol,cap1,[99 57 47]);
end;

%vol=tom_norm(vol,3000);

[thresh,mass,vol_back]=tom_calc_isosurface(-tom_norm(vol,1),720,4.2,0.1);
vol=-((-tom_norm(vol,1))>thresh);

