function [Align]=tom_av3_align_particles(ref,wedge_ref,Align,phi,psi,theta,iteration,mask,filter,mask_ccf,parallel)
%TOM_AV3_ALIGN_PARTICLES creates ...
%
%   [Align]=tom_av3_align_particles(ref,wedge_ref,Align,phi,psi,theta,iteration,mask,filter,mask_ccf,parallel)
%
%PARAMETERS
%
%  INPUT
%   ref                 ...
%   wedge_ref           ...
%   Align               ...
%   phi                 ...
%   psi                 ...
%   theta               ...
%   iteration           ...
%   mask                ...
%   filter              ...
%   mask_ccf            ...
%   parallel            ...
%  
%  OUTPUT
%   Align               ...
%
%EXAMPLE
%   ... = tom_av3_align_particles(...);
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



num_of_tasks=parallel.number_of_tasks;
max_lost_packages=parallel.packageloss;
jobmanager=parallel.jobmanager;
file_path=fileparts(Align(1,1).Filename);

if (1||num_of_tasks==1)
   
    [Align]=tom_av3_align_particles_worker(ref,wedge_ref,Align,phi,psi,theta,iteration,mask,filter,mask_ccf,[1 size(Align,2)])

else
    
    [a,b]=system(['chmod -R ugo+rwx ' file_path '/tmp/']);
     jm = findResource('scheduler','type','jobmanager','lookupurl',jobmanager);
  
    result_p=tom_get_paraell_results(jm,'hosts');
    disp(['Alignment on ' num2str(size(result_p,1)) ' hosts']); 
    for iii=1:size(result_p,1)
        disp(result_p{iii});
    end;    
    clear('result_p');
    j = createJob(jm,'Name','demojob');

    %add path of functions and data
    pdd={'/fs/pool/pool-bmsan-apps/tom_dev/IOfun' '/fs/pool/pool-bmsan-apps/tom_dev/Filtrans/'  '/fs/pool/pool-bmsan-apps/tom_dev/Analysis/' '/fs/pool/pool-bmsan-apps/tom_dev/Misc/' '/fs/pool/pool-bmsan-apps/tom_dev/Sptrans/' ...
        '/fs/pool/pool-bmsan-apps/tom_dev/Geom/' '/fs/pool/pool-bmsan-apps/tom_dev/av2/' '/fs/pool/pool-bmsan-apps/tom_dev/av3_new/' '/fs/sandy01/lv03/pool/bmsan/apps/tom_dev/Display/' };
    set(j,'FileDependencies',pdd)
    pdd2={[file_path]};
    set(j,'PathDependencies',pdd2);
    rehash toolbox;

    %set(j,'FileDependencies',{'/fs/bmsan/apps/tom_dev/av2/tom_av2_multi_ref_alignment.m'})

    packages=tom_calc_packages(num_of_tasks,size(Align,2));
    for i=1:num_of_tasks
        createTask(j,@tom_av3_align_particles_worker,1,{ref,wedge_ref,Align,phi,psi,theta,iteration,mask,filter,mask_ccf,packages(i,1:2)});
    end;
    submit(j);
    %tom_disp_paraell_progress(j,packages(:,3));
    waitForState(j);
    out = getAllOutputArguments(j);
    [result_p errorsum]=tom_get_paraell_results(j);
    if (errorsum > 0);
        for i=1:size(result_p,1)
            disp(['Hostname: ' result_p{i,1} '  Error: ' result_p{i,2}]);
        end;
    end;
    destroy(j);
    [Align lost error_m]=built_Align(out,packages,max_lost_packages);
    

    if (error_m==1)
        return;
    end;
    
    
end;









function [Align]=tom_av3_align_particles_worker(ref,wedge_ref,Align,phi,psi,theta,iteration,mask,filter,mask_ccf,package)

%error([num2str(package(1)) ' ' num2str(package(2))])
Align(iteration+1,:)=Align(iteration,:);

demo_mode=0;
if (demo_mode==1)
    h1=figure;
    h2=figure; 
    h3=figure; 
    h4=figure; 
end;

for i=package(1):package(2)

    disp(num2str(i))
    part=tom_emread(Align(iteration,i).Filename);
    part=part.Value;

    %wedge_part=tom_av3_wedge(ones(size(part)), Align(iteration,i).Tomogram.AngleMin, Align(iteration,i).Tomogram.AngleMax);
    
    wedge_part=tom_wedge(ones(size(part)),30);
    
    [angles,shifts,ccf_peak,rot_part,ccf_max,angle_max]=tom_av3_align(ref,part,phi,psi,theta,mask,filter,mask_ccf,wedge_ref,wedge_part);
    
   
    
    if (demo_mode==1)
        figure(h1);tom_dspcub(ref); set(h1,'Position',[106   543   627   507]); set(h1,'Name','ref'); 
        figure(h2);tom_dspcub(part); set(h2,'Position',[104    34   629   426]); set(h2,'Name','part');
        figure(h3);tom_dspcub(ccf_max); set(h3,'Position',[739   544   622   506]); set(h3,'Name',['ccf ' num2str(ccf_peak) ]);
        figure(h4);tom_dspcub(rot_part); set(h4,'Position',[736    37   632   423]); set(h4,'Name',['rot ref: angles ' num2str(angles) ' shifts  ' num2str(shifts)]);
        drawnow;
    end;
    
    
    Align(iteration+1,i).Angle.Phi=angles(1);
    Align(iteration+1,i).Angle.Psi=angles(2);
    Align(iteration+1,i).Angle.Theta=angles(3);
    rotmatrix=tom_angles2rotmatrix([angles(1) angles(2) angles(3)]);
    Align(iteration+1,i).Angle.Rotmatrix=rotmatrix;
    Align(iteration+1,i).Shift.X=shifts(1);
    Align(iteration+1,i).Shift.Y=shifts(2);
    Align(iteration+1,i).Shift.Z=shifts(3);
% error([num2str(Align(iteration+1,i).Shift.X) ' ' num2str(Align(iteration+1,i).Shift.Y) ' ' num2str(Align(iteration+1,i).Shift.Z)])      
    Align(iteration+1,i).CCC=ccf_peak;
end;


function [Align lost error_m]=built_Align(out,packages,max_lost_packages)

error_m=0;

lost=0;

if (isempty(out)==1)
    error_m=1;
    disp(['ERROR: All Packages lost !']);
    disp(['...restarting current iteration !']);
    error_m=1;
    Align='';
    return;
end;

tmp_st=out{1};
Align=tmp_st;
for i=1:size(out)
    if (isempty(out{i})==1)
        for ii=1:packages(i,3);
            lost=lost+1;
            disp('warning: Package lost');
        end;
        continue;
    end;
    tmp_st=out{i};
    Align(:,packages(i,1):packages(i,2))=tmp_st(:,packages(i,1):packages(i,2));
    
end;

if (lost > max_lost_packages)
    error_m=1;
    disp(['ERROR: Maximum number of lost packages reached !']);
    disp(['...restarting current iteration !']);
   
end;
