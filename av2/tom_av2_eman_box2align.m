function align2d=tom_av2_eman_box2align(f_box,base_p_imges,ext_in)
%tom_av2_eman_box2align transfers a eman box 2 align2d struct 
%
%   tom_av2_eman_box2align(f_box,base_p_imges,ext_in)
%
%  TOM_AV2_EMAN_BOX2ALIGN transfers a eman box 2 align2d struct
% 
%   
%
%PARAMETERS
%
%  INPUT
%   f_box            box file names
%   base_p_imges     base path of the images
%   ext_in           (*.em) extension of the images
%
%  EXAMPLE
%     
%  
% (1) 
% micrographs are numbered like: 1.em,2.em 3.em ...
%
% align2d=tom_av2_eman_box2align('hstpp_eman_*.box','/fs/pool/pool-tpp/TPP_tom_acquisition/TPP_240709_clean_2/792009/high/');
%   
% (2) 
% micrographs are numbered like: test_1.em,test_2.em test_3.em ...
%
% align2d=tom_av2_eman_box2align('hstpp_eman_*.box','/fs/pool/pool-tpp/TPP_tom_acquisition/TPP_240709_clean_2/792009/high/test_');
%
% (3) 
% micrographs are numbered like: test_1.mrc,test_2.mrc test_3.mrc ...
%
% align2d=tom_av2_eman_box2align('hstpp_eman_*.box','/fs/pool/pool-tpp/TPP_tom_acquisition/TPP_240709_clean_2/792009/high/test_','.mrc');  
%
%
%
%  save('pick.mat','align2d');  ...to use in tom_av2_particlepicker   
%
%  NOTE: the boxer files and the micrographs have to have the same
%  numbering sceme
%      
%
%REFERENCES
%
%SEE ALSO
%   tom_av2_particlepicker
%
%   created  by fb
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
    ext_in='.em';
end;

st=dir(f_box);

path_d=pwd;

[base f_name ext]=fileparts(f_box);

align(1,1).dataset = '';
align(1,1).filename = '';
align(1,1).position.x = 0;
align(1,1).position.y = 0;
align(1,1).class = 'default';
align(1,1).radius = 0;
align(1,1).color = [0 1 0];
align(1,1).shift.x = 0;
align(1,1).shift.y = 0;
align(1,1).angle = 0;
align(1,1).isaligned = 0;
align(1,1).ccc = 0;
align(1,1).quality = 0;
align(1,1).normed = 'none';
align(1,1).ref_class = 0;


align2d(1,1:300000)=align;
zz=1;

for i=1:length(st)
    
    if (isempty(base))
        tmp=textread(st(i).name);
    else
        tmp=textread([base '/' st(i).name]);
    end;
    
    for ii=1:size(tmp,1)
        [a b c]=fileparts(st(i).name);
        idx=strfind(b,'_');
        num= b(idx(length(idx))+1:end);
        
        align2d(1,zz).filename=[base_p_imges '' num ext_in];
        if (tom_isemfile(align2d(1,zz).filename)==0)
            disp(['Warning cannot read '  align2d(1,zz).filename]);
        end;
        align2d(1,zz).position.x=tmp(ii,1)+round(tmp(ii,3)./2);
        align2d(1,zz).position.y=tmp(ii,2)+round(tmp(ii,3)./2);
        align2d(1,zz).radius=round(tmp(ii,3)./2);
        zz=zz+1;
    end;
    if (mod(zz,500)==0)
        disp([align2d(1,zz-1).filename ' count ' num2str(zz-1)  ' processed !']);
    end;
end;

tmpp=align2d(1,1:zz-1);

clear('align2d');

align2d=tmpp;



