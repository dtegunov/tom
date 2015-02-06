function [angles,shifts,ccf_peak,rot_part,ccf_max,angle_max]=tom_av3_align(ref,part,phi,psi,theta,mask,filter,mask_ccf,wedge_ref,wedge_part)
%TOM_AV3_ALIGN performs three-dimensional alignment of two volumes ...
% taking the missing wedge of both volumes (particles) into account.
% Cross-correlation is normed between -1:1 and calculated inside a
% mask-volume.
%
%   [angles,shifts,ccf_peak,rot_part,ccf_max,angle_max]=tom_av3_align(ref,part,phi,psi,theta,mask,filter,mask_ccf,wedge_ref,wedge_part)
%
%PARAMETERS
%
%  INPUT
%   ref                 reference particle
%   part                particle to be aligned
%   phi                 vector of Euler-angles to be scanned
%   psi                 vector of Euler-angles to be scanned
%   theta               vector of Euler-angles to be scanned
%   mask                mask volume inside the ccf is calculated
%   filter              ...
%   mask_ccf            mask volume applied to the ccf
%   wedge_ref           missing wedge volume of reference in real space
%   wedge_part          missing wedge volume of particle in real space
%  
%  OUTPUT
%   angles              corresponding angles to maximum ccc
%   shifts              corresponding shifts to maximum ccc
%   ccf_peak            maximum ccc value
%   rot_part            aligned part
%   ccf_max             full maximum ccf function
%   angle_max           full angles (index-encoded) corresponding to
%                        maximum ccf function
%
%EXAMPLE
%   v=tom_av2_build_artificial26S;
%   v2=-tom_bin(v,2);
%   yyy=zeros(size(v2));
%   wedge_ref=tom_wedge(yyy,30);
%   v2_rot_wedge=tom_apply_weight_function(tom_shift(tom_rotate(v2,[0 10 20]),[0 0 0]),wedge_ref);
%   mask=tom_spheremask(ones(size(v2)));
%   [angles,shifts,ccf_peak,rot_part,ccf_max,angle_max]=tom_av3_align(v2,v2_
%   rot_wedge,[0 10 20],[0 10 20],[0 10 20],mask,filter,mask,wedge_ref,wedge_ref);
%
%   Example 2:
%   where=which('20S_full_res.em');
%   s20=tom_emread(where);
%   b=-tom_bin(s20.Value,1);
%   b_rot_shift=tom_shift(tom_rotate(b,[10 -10 20]),[1 2 3 ]);
%   [angles,shifts,ccf_peak,rot_part,ccf_max,angle_max]=tom_av3_align(b,b_rot_shift,[-10 10 20],[-10 10 20],[0 10 20]);
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by SN/FB 06/27/06
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


%parse inputs

if (sum(size(part)==size(ref))~=3)
    error('volumes have different size!');
end;

flag='wedges';

switch nargin
    case 5
        flag='no_wedges';
        mask=ones(size(part));
        mask_ccf=ones(size(part));
        filter.Apply=0;
    case 6
        flag='no_wedges';
        mask_ccf=ones(size(part));
        filter.Apply=0;
    case 7
        mask_ccf=ones(size(part));
        filter.Apply=0;
        flag='no_wedges';
    case 8
        wedge_part=wedge_ref;
    case 10
        flag='wedges';
    otherwise
        error('wrong number of arguments!');
end;

if (isempty(mask))
    mask=tom_spheremask(ones(size(part)));
end;
if (isempty(mask_ccf))
    mask_ccf=ones(size(part));
end;
 

if strcmp(flag,'wedges')
    yyy=zeros(size(part));
    yyy(1,1,1) =1;
    psf_part=real(tom_ifourier(ifftshift(fftshift(tom_fourier(yyy)).*wedge_part)));
end;

ccf_max=-ones(size(part)).*2;
angle_max=-ones(size(part)).*2;


demo_mode=0;
if (demo_mode==1)
h1=figure;set(h1,'DoubleBuffer','on');set(h1,'Position',[106   543   627   507]); set(h1,'Name','ref_rot');
h2=figure; set(h2,'DoubleBuffer','on');set(h2,'Position',[104    34   629   426]); set(h2,'Name','part');
h3=figure; set(h3,'DoubleBuffer','on');set(h3,'Position',[739   544   622   506]); set(h3,'Name',['ccf '  ]);
%h4=figure; 
end;

part=tom_apply_filter(part,filter);

zz=1;
%loop over all angles
for i_phi=phi(1):phi(2):phi(3)
    for i_psi=psi(1):psi(2):psi(3)
        for i_theta=theta(1):theta(2):theta(3)

            [rotmatrix]=tom_angles2rotmatrix([i_phi i_psi i_theta]);
            ref_rot=tom_rotate(ref,[rotmatrix]);
%            tom_emwrite(['out2/ref_' num2str(zz) '.em'],ref_rot);
            zz=zz+1;
            
            if strcmp(flag,'no_wedges')
                ccf=tom_corr(ref_rot,part,'norm',mask);
                 if (demo_mode==1)
                    figure(h1);tom_dspcub(ref_rot); 
                    figure(h2);tom_dspcub(part); 
                    figure(h3);tom_dspcub(ccf); 
                    %h4=figure; tom_dspcub(rot_part); set(h4,'Position',[736    37   632   423]); set(h4,'Name',['rot ref: angles ' num2str(angles) ' shifts  ' num2str(shifts)]);
                    drawnow;
                end;
            
            else
                wedge_ref_rot=tom_rotate(wedge_ref,[i_phi i_psi i_theta]);
                psf_ref=real(tom_ifourier(ifftshift(fftshift(tom_fourier(yyy)).*wedge_ref_rot)));
                ref_rot=tom_apply_filter(ref_rot,filter);
                ref_rot=ref_rot.*mask;
                part=part.*mask;
                ccf=tom_corr(ref_rot,part,'norm',mask,psf_ref,psf_part);
                

                if (demo_mode==1)
                    figure(h1);tom_dspcub(ref_rot); 
                    figure(h2);tom_dspcub(part); 
                    figure(h3);tom_dspcub(ccf); 
                    %h4=figure; tom_dspcub(rot_part); set(h4,'Position',[736    37   632   423]); set(h4,'Name',['rot ref: angles ' num2str(angles) ' shifts  ' num2str(shifts)]);
                    drawnow;
                end;


            end;

            ccf_max=(ccf>ccf_max).*ccf+(ccf<=ccf_max).*ccf_max;

            angle_index=tom_av3_angle2index([i_phi i_psi i_theta], [phi(1) psi(1) theta(1)], [phi(2) psi(2) theta(2)],[phi(3) psi(3) theta(3)]);
            angle_max=(ccf>=ccf_max).*angle_index+(ccf<ccf_max).*angle_max;

        end;
    end;
end;

[pos ccf_peak]=tom_peak(ccf_max.*mask_ccf);

[angles] = tom_av3_index2angle(angle_max(pos(1),pos(2),pos(3)), [phi(1) psi(1) theta(1)], [phi(2) psi(2) theta(2)],[phi(3) psi(3) theta(3)]);


%get subpixel
rot_ref=tom_rotate(ref,angles);
ccf_s=tom_corr(rot_ref,part,'norm',mask);
[pos ccf_peak]=tom_peak(ccf_s.*mask_ccf,'spline');


shifts=pos-(floor(size(part)./2)+1);
  %  error([num2str(shifts(1)) ' ' num2str(shifts(2)) ' ' num2str(shifts(3))])

rot_part= double(tom_rotate(tom_shift(part,-shifts),[-angles(1) -angles(2) -angles(3)]));





