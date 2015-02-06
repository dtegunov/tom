function tom_av3_build_phantom_particle(dest_folder,num_of_particles,wedge)
%TOM_AV3_BUILD_PHANTOM_PARTICLE creates ...
%
%   tom_av3_build_phantom_particle(dest_folder,num_of_particles,wedge)
%
%PARAMETERS
%
%  INPUT
%   dest_folder         ...
%   num_of_particles    ...
%   wedge               ...
%  
%  OUTPUT
%
%EXAMPLE
%   tom_av3_build_phantom_particle(...);
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



part=tom_av2_build_artificial26S('2caps');
%save part part
%load part
part=tom_bin(part,2);

w_func=tom_wedge(ones(size(part)),wedge);
w_func=w_func.*tom_spheremask(ones(size(w_func)));

yyy=zeros(size(part));
yyy(1,1,1) =1;
psf=real(tom_ifourier(ifftshift(fftshift(tom_fourier(yyy)).*w_func)));

for i=11:10+num_of_particles
    name=[dest_folder '/part_' num2str(i) '.em'];
    shifts_and_angles(i,1:3)=round((rand(3,1).*360)./60).*60;

    shifts_and_angles(i,4:6)=round(rand(3,1).*3);
%    shifts_and_angles(i,4:6)=[0 0 0];
    
    tmp_part=tom_rotate(part,shifts_and_angles(i,1:3));
    rotmatrix=tom_angles2rotmatrix([shifts_and_angles(i,1:3)]);
    tmp_part=tom_shift(tmp_part,shifts_and_angles(i,4:6));

%    tmp_part=tom_apply_weight_function(tmp_part,w_func);

    tmp_part=tom_ifourier(tom_fourier(tmp_part).*tom_fourier(psf));
    
    tmp_part=tom_emheader(tmp_part);
    
    tmp_part.Header.Comment=[num2str(shifts_and_angles(i,1)) ' ' num2str(shifts_and_angles(i,2)) ' ' num2str(shifts_and_angles(i,3)) '  ' num2str(shifts_and_angles(i,4)) ' ' num2str(shifts_and_angles(i,5)) ' ' num2str(shifts_and_angles(i,6))  ]; 
    tom_emwrite(name,tmp_part.Value+(rand(size(tmp_part.Value))./3 )) ;
    tom_emwrite([dest_folder '/rotmatrix_' num2str(i) '.em'],rotmatrix);
    Align(1,i).Angle.Phi=shifts_and_angles(i,1);
    Align(1,i).Angle.Psi=shifts_and_angles(i,2);
    Align(1,i).Angle.Theta=shifts_and_angles(i,3);
    Align(1,i).Shift.X=shifts_and_angles(i,4);
    Align(1,i).Shift.Y=shifts_and_angles(i,5);
    Align(1,i).Shift.Z=shifts_and_angles(i,6);
    
    
    
    disp(name);
    
end;

tom_emwrite('angles',shifts_and_angles);
save('p_align3','Align');
