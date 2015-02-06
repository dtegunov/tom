function [ang_back shift_back sym_struct aligned_sym_struct]=tom_get_sym_axis(volume,sym,iter,phi,psi,theta,num_ref,ini_ang,ini_shift,mask,output)
%TOM_GET_SYM_AXIS calculates symmetry axis
%
%   [ang shift aligned_struct sym_struct aligned_sym_struct]=tom_get_sym_axis(volume,sym,iter,phi,psi,theta,num_ref,ini_ang,ini_shift,mask,output)
%
%  tom_get_sym_axis calculates symmetry axis by iterativ symmetrisation and fitting 
%  
%  
%  
%
%PARAMETERS
%
%  INPUT
%   volume                input volume
%   sym                   symmetry
%   iter                  number of iterations
%   phi                   start incre stop
%   psi                   start incre stop
%   theta                 start incre stop
%   num_ref               number of refinement steps
%   ini_ang               (initial angle) 
%   ini_shift             (intital shift)  
%   mask                   mask for fitting 
%   output                (opt.) outputfolder for intermed results
%
% OUTPUT 
% 
%   ang_back            angle of sym axis 
%   shift_back          shift of sym axis 
%   sym_struct          sym struct     
%   aligned_sym_struct  aligned sym struct
%
%EXAMPLE
%    
% [ang_back shift_back sym_struct aligned_sym_struct]=tom_get_sym_axis(vol,7,3,[-8 4 8],[-8 4 8],[-8 4 8],5);
% 
% %apply axis and shift
%  my_aligned_struct=tom_rotate(tom_shift(sym_struct,shift_back),[ang_back(2) ang_back(1) ang_back(3)]);
% 
%
%
% %MODEL EXAMPLE:
%
% %build 8-fold model
% cyl=tom_cylindermask(ones(64,64,64),8);
% cyl(:,:,1:15)=0;
% cyl(:,:,65-15:end)=0;
% sp=tom_spheremask(ones(64,64,64),5,0,[42 42 33]);
% symm=tom_symref(sp+cyl,8)>0;
% 
% % orient along y-axis
% symm_rot=tom_rotate(tom_rotate(symm,[270 90 90]),[90 0 0]);
%
% %apply shift and rotation
%
% symm_rot_trans=tom_shift(tom_rotate(symm,[0 0 17]),[3 2 0]);
% [ang_back shift_back sym_struct aligned_sym_struct]=tom_get_sym_axis(symm_rot_trans,8,3,[-32 16 32],[-32 16 32],[-32 16 32],5);
% %apply axis and shift
% 
% my_aligned_struct=tom_shift(tom_rotate(sym_struct,[ang_back(1) ang_back(2) ang_back(3)]),shift_back); 
%
%
%NOTE:
% 
% symmetry axis has to be y-axis!! (zxz notation)
% 
%
%REFERENCES
%
%SEE ALSO
%   
%
%   created by fb 
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


sz=size(volume);

if (nargin < 7)
    num_ref=5;
end;

if (nargin < 8)
   ini_ang=[0 0 0];
end;

if (nargin < 9)
    ini_shift=[0 0 0];
end;

if (nargin <10 || strcmp(mask,'default'))
    mask=tom_spheremask(ones(sz),round(sz(1)./2)-2,1);
end;

if (nargin <11)
    output='out_sym';
end;
    
try
    warning off;
    mkdir(output); 
    mkdir([output '/debug/']);
    mkdir([output '/output/']);
    warning on;
catch ME
    
end;

%transfor 2 y-axis
volume_org=volume;
volume=tom_rotate(volume,[-90 180 90]);

al_part=volume;


tom_emwrite([output '/aling_' num2str(0) '_' num2str(0) '.em'],al_part);


phi_work=phi;
psi_work=psi;
theta_work=theta;

all_angles=ini_ang;
all_shiftvec=ini_shift;

for ii=1:num_ref
    [ang_sum shift_sum rott]=tom_sum_rotation(all_angles,all_shiftvec);
    org=tom_shift(tom_rotate(volume,-ang_sum),-shift_sum);
    disp(['Refinement ' num2str(ii) ':  phi: ' num2str(phi_work(1)) ':' num2str(phi_work(2)) ':' num2str(phi_work(3)) ...
                                    ':  psi: ' num2str(psi_work(1)) ':' num2str(psi_work(2)) ':' num2str(psi_work(3)) ...             
                                    ':  theta: ' num2str(theta_work(1)) ':' num2str(theta_work(2)) ':' num2str(theta_work(3)) ] );
    al_part=org;
    for i=1:iter
        
        %apply symmetry
        new_sym=tom_symref(al_part,sym,'y-axis');
        tom_emwrite([output '/sym_' num2str(ii) '_' num2str(i) '.em'],new_sym);
        
        %fit structure back
        [angles,shiftvec,ccc,al_part,ccf_max,angle_max]=tom_av3_align(new_sym,org,phi_work,psi_work,theta_work,mask);
        
        tom_emwrite([output '/aling_' num2str(ii) '_' num2str(i) '.em'],al_part);
        
        disp([ '    '  num2str(i) ':  Shifts: ' num2str(shiftvec) '  Angles: ' num2str(angles)  '  ccf: ' num2str(ccc)  ]);
    end;
    disp(' ');
    all_angles(ii,:)=angles;
    all_shiftvec(ii,:)=shiftvec;
    phi_work=phi_work./2;
    psi_work=psi_work./2;
    theta_work=theta_work./2;
end;


[ang_sum shift_sum rott]=tom_sum_rotation(all_angles,all_shiftvec);
shift_sum=shift_sum';

%ang_sum=all_angles;
%shift_sum=all_shiftvec;

%transfer return values!
%[ang_back shift_back rott]=tom_sum_rotation([ang_sum(2) ang_sum(1) ang_sum(3); -180 90 -90],[shift_sum;0 0 0]);

[ang_back shift_back rott]=tom_sum_rotation([-90 180 90; ang_sum(2) ang_sum(1) ang_sum(3); -180 90 -90],[0 0 0 ;shift_sum;0 0 0]);
% ang_back=ang_sum;
% shift_back=shift_sum;

sym_struct=tom_rotate(new_sym,[-180 90 -90]);

%aligned_sym_struct_tmp=tom_rotate(tom_shift(new_sym,-shift_back),[ang_back(1) ang_back(2) ang_back(3)]);
aligned_sym_struct_tmp=tom_shift(tom_rotate(sym_struct,[ang_back(1) ang_back(2) ang_back(3)]),shift_back);
aligned_sym_struct=aligned_sym_struct_tmp;

%mesure shift finally quality control
cc=tom_corr(volume_org,aligned_sym_struct_tmp,'norm');
[pos val]=tom_peak(cc,'spline');

sh=pos-(floor(size(volume)./2)+1);
disp(' ');
disp(['back fit shift: ' num2str(sh)]);

disp('done!');


%lagacy

%apply transformation
%al_part=tom_shift(org,-shiftvec);
%al_part=tom_rotate(al_part,[-angles(2) -angles(1) -angles(3)]);
%al_part= double(tom_rotate(tom_shift(org,-shiftvec),[-angles(1)
%-angles(2) -angles(3)]));

%tom_emwrite(['tmp_models/upper_s' num2str(app_sym) '_iter_' num2str(i)' '.em' ],new_sym_cut); %tom_emwrite(['tmp_models/upper_s' num2str(3) '_iter_' num2str(i)' '.em' ],new_sym_cut_3t);
%st_upper(app_sym).ccf(i)=ccc;
%st_upper(app_sym).angles(i,:)=angles;%st_upper(app_sym).shifts(i,:)=shiftvec;

%   sym_cc=tom_corr(tom_rotate(al_part,round(360./sym)),al_part,'norm');
%         [aa bb]=tom_peak(sym_cc);



 %correct for b-factor
%         [bfactor decay_restore corrected_org]=tom_fit_bfactor(new_sym,4.42,[5 20],[11 Inf],2500000,1,0);
%         tom_emwrite(['out/sym_corr_' num2str(ii) '_' num2str(i)  '.em'],corrected_org);



%shift_back=shift_back-sh';
%aligned_sym_struct=tom_shift(aligned_sym_struct,-sh);
%aligned_sym_struct=tom_shift(tom_rotate(new_sym,[ang_sum(2) ang_sum(1) ang_sum(3)]),shift_back); 
%org=tom_shift(tom_rotate(volume,-ang_sum),-shift_sum);

        