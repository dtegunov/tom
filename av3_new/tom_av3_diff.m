function [diff vols_aligned all]=tom_av3_diff(vols,range,flag_save,flag)
% TOM_AV3_DIFF calcs diff between n volumes
%  
%     [diff cc sh]=tom_av3_diff(vol1,vol2)
%  
%  PARAMETERS
%  
%    INPUT
%     vol1                 cell file with volnames or wildcard string,
%                          supports EM and Spider file format
%     flag_save            save difference and aligned volumes in EM format
%     flag                 0=rotation+shift on, 1=rotation off, 2=shift off
%     
%
%    OUTPUT
%     diff                 4D-difference
%     all                  all alignment info data
%    
%  
%  EXAMPLE
%     range=[-3 3 3]; % search range, from -3 degree to +3 degree,
%                       % stepsize 3 degree.
%     vols{1}='1.vol'   % 1st volume name 
%     vols{2}='v2_10_22_33.vol'   % 2nd volume name 
%     [diff vols_aligned all]=tom_av3_diff(vols,range);
%
%     % another example:
%     vol=tom_emread('1.vol');
%     v2=tom_shift(tom_mex_rotate(vol.Value,[-10 20 30]),[5 6 7]);
%     mask = tom_spheremask(ones(size(v2)),size(v2,1)./2-3,2,[size(v2)./2+1]);
%     tom_emwrite(['v2_shift_5_6_7_rot_m10_20_30.vol'],v2)
%     m{1}='1.vol'
%     m{2}='v2_shift_5_6_7_rot_m10_20_30.vol'
%     range=[-30 10 30]; % search range, from -30 degree to +30 degree
%     [diff vols_aligned all]=tom_av3_diff(m,range);
%
%  REFERENCES
%  
%  SEE ALSO
%     ...
%  
%     created by SN/FB 03/30/10
%  
%     Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%     Journal of Structural Biology, 149 (2005), 227-234.
%  
%     Copyright (c) 2004-2007
%     TOM toolbox for Electron Tomography
%     Max-Planck-Institute of Biochemistry
%     Dept. Molecular Structural Biology
%     82152 Martinsried, Germany
%     http://www.biochem.mpg.de/tom

if nargin<4
    flag=0;
end;
if nargin<3
    flag_save=0;
end;

start=range(1);
step=range(2);
stop=range(3);

sh_out=0;
% SYNOPSIS: [out] = tom_mex_corr3d([datatype], size, templates, particles, angles, mask, wedge_ref, wedge_part, type, [...])
% mex-file compiled from file "/fs/sandy01/lv03/pool/bmsan/haller/DA/tomc/src/mex/tom_mex_corr3d.cpp" at Feb 11 2008, 16:59:35

idx=0;
for phi=start:step:stop
    for psi=start:step:stop
        for theta=start:step:stop
            idx=idx+1;
            angles(idx,:)=[phi, psi, theta];
        end;
    end;
end;
disp('---------------------------------------------------------------------------');
disp('---------------------------Calculate alignment-----------------------------');
disp(' ');

disp(['Total number of angles screened: ' num2str(idx)]);
disp(' ');


idx_all=0;
for i=1:length(vols)
        [type]=tom_determine_file_type2(vols{i});
        switch type
            case 'em'
                in=tom_emreadc(vols{i});
                in1=in.Value;
            case 'spider'
                in=tom_spiderread(vols{i});
                in1=in.Value;
        end;                
    for ii=1:length(vols)
        [type]=tom_determine_file_type2(vols{ii});
        switch type
            case 'em'
                in=tom_emreadc(vols{ii});
                in2=in.Value;
            case 'spider'
                in=tom_spiderread(vols{ii});
                in2=in.Value;
        end;
        disp(['Align: ' vols{i} ]);
        disp(['and  : ' vols{ii} ]);
        disp(' ');

        mask = tom_spheremask(ones(size(in1)),size(in1,1)./2-3,2,[size(in1)./2+1]);
        [out] = tom_mex_corr3d('single', [size(in1)],in1, in2, angles',mask,[],[],'peak');
        idx_all=idx_all+1;
        all{idx_all}.peak=double(out.peak);
        all{idx_all}.ncc_val=out.ncc_val;
        all{idx_all}.angle=double(out.angle);
        all{idx_all}.angle_idx=out.angle_idx;
        all{idx_all}.template=vols{ii};
        all{idx_all}.particle=vols{i};
        disp(['Peak: ' num2str(all{idx_all}.peak) ' NCC: ' num2str(all{idx_all}.ncc_val) ' Angle: ' num2str(all{idx_all}.angle)]);
    end;
end;

disp('---------------------------------------------------------------------------');
disp('-------------------------Calculate differences-----------------------------');
disp(' ');

idx_all=0;
for i=1:length(vols)
        [type]=tom_determine_file_type2(vols{i});
        switch type
            case 'em'
                in=tom_emreadc(vols{i});
                in1=in.Value;
            case 'spider'
                in=tom_spiderread(vols{i});
                in1=in.Value;
        end;                
    for ii=1:length(vols)
        [type]=tom_determine_file_type2(vols{ii});
        switch type
            case 'em'
                in=tom_emreadc(vols{ii});
                in2=in.Value;
            case 'spider'
                in=tom_spiderread(vols{ii});
                in2=in.Value;
        end;
        disp(' ');
        disp(['Calculate difference: ' vols{i} ]);
        disp(['and                 : ' vols{ii} ]);        
        idx_all=idx_all+1;
        angles=all{idx_all}.angle;
        peak=all{idx_all}.peak;
        %shifts=size(in1)-peak;
        shifts=peak;
        disp(['Volume : ' vols{i} ' Angles: ' num2str(angles) ' Shifts: '  num2str(shifts)]);
        disp(' ');
        in1_rot=tom_mex_rotate(in1,angles);
        in1_shift=tom_shift(double(in1_rot),shifts);
%        in1_rot=tom_rotate(in1,[-angles(2) -angles(1) -angles(3)],'linear');
        diff(:,:,:,idx_all)=(tom_norm(in1_shift,'mean0+1std',mask)-tom_norm(in2,'mean0+1std',mask)).*mask;        
        vols_aligned(:,:,:,idx_all)=in1_shift.*mask;
        if flag_save==1
            tom_emwrite(['diff_vol_' num2str(i) '_' num2str(ii) '.em'],diff(:,:,:,idx_all));
            tom_emwrite(['aligned_vol_' num2str(i) '_' num2str(ii) '.em'],vols_aligned(:,:,:,idx_all));
            disp('Volumes saved.')
            disp(' ')
        end;
    end;
end;
