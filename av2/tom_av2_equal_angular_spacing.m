function angles=tom_av2_equal_angular_spacing(theta_range,phi_range,increment,method)

% TOM_AV2_EQUAL_ANGULAR_SPACING computes an quasi-equal angular sampling
%
%   angles=tom_av2_equal_angular_spacing(theta_range,phi_range,increment,method)
%
% PARAMETERS
%   INPUT
%
%   theta_range        min and max of theta (0 to 90)
%   phi_range          min and max of phi (0 to 359.9)
%   increment          angular increment
%   method             'spider' or 'xmipp', 'xmipp_all'
%
%   OUTPUT
%
%   angles              two dimensional vector with angles
%
%   TOM_AV2_EQUAL_ANGULAR_SPACING computes an quasi-equal angular sampling suitable
%   for projection matching purposes.
%
%EXAMPLE
%   angles=tom_av2_equal_angular_spacing([0 90],[0 359.9],15,'spider');
%   polar(angles(:,1).*pi./180,angles(:,2).*pi./180,'o');
%
%SEE ALSO
%   adapted from a Spider routine written by P. Penczek by SN
%
%   created by SN 10/05/07
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

%PSI = 0.0;

if isequal(method,'spider')

    if (theta_range(1)==0.0 && theta_range(2)==00) || theta_range(1)>theta_range(2)
        theta_range(1)=0.0;
        theta_range(2)=90.0;
    end;
    if (phi_range(1)==0.0 && phi_range(2)==00) || phi_range(1)>phi_range(2)
        phi_range(1)=0.0;
        phi_range(2)=359.9;
    end;

    if theta_range(1)<90 && theta_range(2)==90.0 && phi_range(1)==0.0 && phi_range(2)>0
        SKIP=1;
    else
        SKIP=0;
    end;
    LITER=0;
    NLIST=4;

    for THETA=theta_range(1):increment:theta_range(2)
        if THETA==0.0 || THETA==180.0
            DETPHI=360.0;
            LT=1;
        else
            DETPHI=increment./sind(THETA);
            LT=max(floor((phi_range(2)-phi_range(1))./DETPHI)-1,1 ); % ???
            DETPHI=(phi_range(2)-phi_range(1))./LT;
        end;
        for I=1:LT
            PHI=phi_range(1)+(I-1).*DETPHI;
            if SKIP && THETA==90.0 && PHI > 180.0
               break;
            end;
            LITER=LITER+1;
            angles(LITER,1)=PHI;
            angles(LITER,2)=THETA;
        end;
    end;

elseif isequal(method,'xmipp')
    
    sampling=increment;
    irot=0;

    tilt_nstep=round(180./sampling)+1;

    tilt_sam=(180./tilt_nstep);

    for tilt_step=0:tilt_nstep
        tilt=(tilt_step./(tilt_nstep-1)).*180;

        if (tilt>0)
            rot_nstep=ceil(360.*sind(tilt)./sampling);

        else
            rot_nstep=1;
        end;
        rot_sam=360./rot_nstep;

        for rot=0:rot_sam:360-rot_sam


            irot=irot+1;
            angles(irot,1)=rot;
            angles(irot,2)=tilt;
        end;
    end;
    
   elseif isequal(method,'xmipp_all') 
    
    sampling=increment;
    irot=0;

    tilt_nstep=round(180./sampling)+1;

    tilt_sam=(180./tilt_nstep);

    for tilt_step=0:tilt_nstep
        tilt=(tilt_step./(tilt_nstep-1)).*180;

        if (tilt>0)
            rot_nstep=ceil(360.*sind(tilt)./sampling);

        else
            rot_nstep=1;
        end;
        rot_sam=360./rot_nstep;

        for rot=0:rot_sam:360-rot_sam


            irot=irot+1;
            angles(irot,1)=rot;
            angles(irot,2)=tilt;
        end;
    end;
elseif isequal(method,'xmipp_doc')
    if (theta_range(1)==0.0 && theta_range(2)==00) || theta_range(1)>theta_range(2)
        theta_range(1)=0.0;
        theta_range(2)=90.0;
    end;
    if (phi_range(1)==0.0 && phi_range(2)==00) || phi_range(1)>phi_range(2)
        phi_range(1)=0.0;
        phi_range(2)=359.9;
    end;
    
    if theta_range(1)<90 && theta_range(2)==90.0 && phi_range(1)==0.0 && phi_range(2)>0
        SKIP=1;
    else
        SKIP=0;
    end;
    LITER=0;
    NLIST=4;
    
    for THETA=theta_range(1):increment:theta_range(2)
        if THETA==0.0 || THETA==180.0
            DETPHI=360.0;
            LT=1;
        else
            DETPHI=increment./sind(THETA);
            LT=max(floor((phi_range(2)-phi_range(1))./DETPHI)-1,1 ); % ???
            DETPHI=(phi_range(2)-phi_range(1))./LT;
        end;
        for I=1:LT
            PHI=phi_range(1)+(I-1).*DETPHI;
            if SKIP && THETA==90.0 && PHI > 180.0
                break;
            end;
            LITER=LITER+1;
            angles(LITER,1)=PHI-180;
            angles(LITER,2)=THETA;
        end;
    end;
    
    
    
    
else
    error('Method not supported');
end


