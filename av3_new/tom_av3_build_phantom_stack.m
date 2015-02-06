function Align=tom_av3_build_phantom_stack(num_of_particles)
%TOM_AV3_BUILD_PHANTOM_STACK creates ...
%
%   Align=tom_av3_build_phantom_stack(num_of_particles)
%
%PARAMETERS
%
%  INPUT
%   num_of_particles    ...
%  
%  OUTPUT
%   Align               ...
%
%EXAMPLE
%   ... = tom_av3_build_phantom_stack(...);
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


%build artificial templates
wedge=30;


part1=tom_bin(tom_av2_build_artificial26S('nocap'),3);
part2=tom_bin(tom_av2_build_artificial26S('1cap'),3);
%part3=tom_bin(tom_av2_build_artificial26S('2caps'),3);
part3=tom_spheremask(ones(20,20,20),5)==0;

zz=1;

for i=1:num_of_particles
    
    ang=(rand(1,3).*360);
    sh=(rand(1,3).*(size(part1,1)./15));
%     ang=[30 30 30];
%     sh=[2 1 1];
    
    new_part=tom_rotate(part1,ang);
    new_part=tom_apply_weight_function(new_part,tom_wedge(new_part,30));
    %new_part=tom_shift(new_part,sh);
    new_part=tom_rotate(new_part,[-ang(2) -ang(1) -ang(3)]);
    new_part=new_part+(rand(size(new_part))./5);
    new_part=tom_norm(new_part+100,'phase');
    new_part=new_part.*tom_spheremask(ones(size(new_part)),((size(new_part,1))./2-1));
    tom_emwrite(['particle_' num2str(zz) '.em'],new_part);
    Align(1,zz).Angle.Phi=ang(1);
    Align(1,zz).Angle.Psi=ang(2);
    Align(1,zz).Angle.Theta=ang(3);
    Align(1,zz).Shift.X=sh(1);
    Align(1,zz).Shift.Y=sh(2);
    Align(1,zz).Shift.Z=sh(3);
    
    Align(1,zz).Angle.Phi=0;
    Align(1,zz).Angle.Psi=0;
    Align(1,zz).Angle.Theta=0;
    Align(1,zz).Shift.X=0;
    Align(1,zz).Shift.Y=0;
    Align(1,zz).Shift.Z=0;
    
    zz=zz+1;
    
    
    ang=(rand(1,3).*360);
    sh=(rand(1,3).*(size(part2,1)./15));
    %ang=[30 30 30];
    %sh=[2 1 1];
    
    new_part=tom_rotate(part2,ang);
    new_part=tom_apply_weight_function(new_part,tom_wedge(new_part,30));
    %new_part=tom_shift(new_part,sh);
    new_part=tom_rotate(new_part,[-ang(2) -ang(1) -ang(3)]);
    new_part=new_part+(rand(size(new_part))./5);
    new_part=tom_norm(new_part+100,'phase');
    new_part=new_part.*tom_spheremask(ones(size(new_part)),((size(new_part,1))./2-1));
    tom_emwrite(['particle_' num2str(zz) '.em'],new_part);
    Align(1,zz).Angle.Phi=ang(1);
    Align(1,zz).Angle.Psi=ang(2);
    Align(1,zz).Angle.Theta=ang(3);
    Align(1,zz).Shift.X=sh(1);
    Align(1,zz).Shift.Y=sh(2);
    Align(1,zz).Shift.Z=sh(3);
    
    Align(1,zz).Angle.Phi=0;
    Align(1,zz).Angle.Psi=0;
    Align(1,zz).Angle.Theta=0;
    Align(1,zz).Shift.X=0;
    Align(1,zz).Shift.Y=0;
    Align(1,zz).Shift.Z=0;
    
    zz=zz+1;
    
    ang=(rand(1,3).*360);
    sh=(rand(1,3).*(size(part3,1)./15));
%     ang=[30 30 30];
%     sh=[2 1 1];
    
    new_part=tom_rotate(part3,ang);
    new_part=tom_apply_weight_function(new_part,tom_wedge(new_part,30));
    %new_part=tom_shift(new_part,sh);
    new_part=tom_rotate(new_part,[-ang(2) -ang(1) -ang(3)]);
    new_part=new_part+(rand(size(new_part))./5);
    new_part=tom_norm(new_part+100,'phase');
    new_part=new_part.*tom_spheremask(ones(size(new_part)),((size(new_part,1))./2-1));
    tom_emwrite(['particle_' num2str(zz) '.em'],new_part);
    Align(1,zz).Angle.Phi=ang(1);
    Align(1,zz).Angle.Psi=ang(2);
    Align(1,zz).Angle.Theta=ang(3);
    Align(1,zz).Shift.X=sh(1);
    Align(1,zz).Shift.Y=sh(2);
    Align(1,zz).Shift.Z=sh(3);
    
    Align(1,zz).Angle.Phi=0;
    Align(1,zz).Angle.Psi=0;
    Align(1,zz).Angle.Theta=0;
    Align(1,zz).Shift.X=0;
    Align(1,zz).Shift.Y=0;
    Align(1,zz).Shift.Z=0;
    
    zz=zz+1;
    
    
%     ang=(rand(1,3).*360);
%     sh=(rand(1,3).*(size(part3,1)./15));
%     ang=[30 30 30];
%     sh=[2 1 1];
%     
%     new_part=tom_rotate(part3,ang);
%     new_part=tom_apply_weight_function(new_part,tom_wedge(new_part,30));
%     new_part=tom_shift(new_part,sh);
%     new_part=new_part+(rand(size(new_part))./15);
%     tom_emwrite(['particle_' num2str(zz) '.em'],new_part);
%     Align(1,zz).Angle.Phi=ang(1);
%     Align(1,zz).Angle.Psi=ang(2);
%     Align(1,zz).Angle.Theta=ang(3);
%     Align(1,zz).Shift.X=sh(1);
%     Align(1,zz).Shift.Y=sh(2);
%     Align(1,zz).Shift.Z=sh(3);
%     zz=zz+1;
    
    
    
    
end;
    