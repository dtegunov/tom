function res=tom_av2_measure_modell(align2d,num_of_shells,iter,step)
%TOM_AV2_MEASURE_MODELL creates ...
%
%   res=tom_av2_measure_modell(align2d,num_of_shells,iter,step)
%
%PARAMETERS
%
%  INPUT
%   align2d             ...
%   num_of_shells       ...
%   iter                ...
%   step                ...
%  
%  OUTPUT
%   res                 ...
%
%EXAMPLE
%   ... = tom_av2_measure_modell(...);
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

if (nargin==2)
    iter(1)=1;
    iter(2)=size(align2d,1);
end;


disp('test');

sz_align=size(align2d,2);
zz=1;

for i=iter(1):step:iter(2)
%     align2d_1=align2d(i,1:round(sz_align./2));
%     align2d_1(1,1).stack_size=[160 160 round(sz_align./2)];
%     align2d_2(1,1)=align2d(i,1);
%     align2d_2(1,1).angleclass=align2d(i,round(sz_align./2)+2).angleclass;
%     align2d_2(1,2:sz_align./2)=align2d(i,round(sz_align./2)+2:sz_align);
%     align2d_2(1,1).stack_size=[160 160 round(sz_align./2)];
    
    mod1=tom_av2_backproj_test(align2d,i,0,[1 round(sz_align./2)]);
    mod2=tom_av2_backproj_test(align2d,i,0,[round(sz_align./2)+2 sz_align]);
    res(:,:,zz)=tom_compare(mod1,mod2,num_of_shells);
    save('res_s','res');
    zz=zz+1;
    zz
end;

% zz=1;
% 
%  for i=iter(1):step:iter(2)
%      mod1=tom_emread(['model/model_' num2str(i)] );
%      mod2=tom_emread(['model/model_' num2str(i+1)] );
%      res(:,:,zz)=tom_compare(mod1.Value,mod2.Value,20);
%      figure;
%      zz=zz+1;
%  end;
%  
%  
 
 
 
 
 
 
 
 
 
 
 
 