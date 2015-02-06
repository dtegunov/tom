function [Align_new,new_average]=tom_av3_applyalignment(Align,iteration_nr,particle_list,new_particle_path,new_particle_extension)
%TOM_AV3_APPLYALIGNMENT creates an average. Apply alignment to ...
%e.g. unbinned particles
%
%   [Align_new,new_average]=tom_av3_applyalignment(Align,iteration_nr,particle_list,new_particle_path,new_particle_extension)
%
%preselect particles by their index: serie1_list=0;ii=1;for i=1:size(serie1,2);if serie1(i)==1;serie1_list(ii)=i;ii=ii+1;end;end;
%
%check for Tempfilename and Filename ... TODO
%
%PARAMETERS
%
%  INPUT
%   Align                       ...
%   iteration_nr                ...
%   particle_list               ...
%   new_particle_path           ...
%   new_particle_extension      ...
%  
%  OUTPUT
%   Align_new                   ...
%   new_average                 ...
%
%EXAMPLE
% [anew,average]=tom_av3_applyalignment(Align,4,[1:140],'/fs/montreal/pub/2
% 6S/alignment3D/data_Oana_1/particles/bin1/part26S_','.em');
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by SN mm/dd/yy
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

part_alig=tom_reademheader(Align(1,1).Tempfilename);
part_new=tom_reademheader([new_particle_path num2str(particle_list(1)) new_particle_extension]);
binning_factor=part_new.Header.Size(1)./part_alig.Header.Size(1);


ii=1;
for i=1:size(particle_list,2)
    
    idx=particle_list(i);
    
    unterstrich=findstr(Align(1,idx).Tempfilename,'_');

    try
        punkt=findstr(Align(1,idx).Tempfilename,'.');
    catch
        punkt=size(Align(1,idx).Tempfilename);
    end;
    if isempty(punkt); punkt=size(Align(1,idx).Tempfilename)+1; end;
    unterstrich_index=unterstrich(size(unterstrich,2));
    punkt_index=punkt(size(punkt,2));
    
    fn = Align(1,idx).Tempfilename;
    
    filenr=fn(unterstrich_index+1:punkt_index-1);
    
    
    Align_new(1,ii)=Align(iteration_nr,idx);
    Align_new(1,ii).Filename=[new_particle_path num2str(filenr) new_particle_extension];

    for iteration=1:iteration_nr
        angles(iteration,1)=Align(iteration,idx).Angle.Phi;
        angles(iteration,2)=Align(iteration,idx).Angle.Psi;
        angles(iteration,3)=Align(iteration,idx).Angle.Theta;
        shifts(iteration,1)=Align(iteration,idx).Shift.X.*binning_factor;
        shifts(iteration,2)=Align(iteration,idx).Shift.Y.*binning_factor;
        shifts(iteration,3)=Align(iteration,idx).Shift.Z.*binning_factor;
    end;    
    
    [euler_out shift_out rott]=tom_sum_rotation(angles,shifts);
    
    Align_new(1,ii).Angle.Phi=euler_out(1);
    Align_new(1,ii).Angle.Psi=euler_out(2);
    Align_new(1,ii).Angle.Theta=euler_out(3);
    Align_new(1,ii).Angle.Rotmatrix=rott;
    
    Align_new(1,ii).Shift.X=shift_out(1);
    Align_new(1,ii).Shift.Y=shift_out(2);
    Align_new(1,ii).Shift.Z=shift_out(3);
    
    
    
    ii=ii+1;
    
end;

new_average=tom_av3_average(Align_new,'sum');

%new_average=tom_av3_average(Align_new);