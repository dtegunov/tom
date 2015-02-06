function model=tom_av2_reconstruction(align2d)
%TOM_AV2_RECONSTRUCTION creates ...
%
%   model=tom_av2_reconstruction(align2d)
%
%**************************************************************************
% Welcome to main program
%**************************************************************************
%
%PARAMETERS
%
%  INPUT
%   align2d             ...
%  
%  OUTPUT
%   model               ...
%
%EXAMPLE
%   ... = tom_av2_reconstruction(...);
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


%transfer often used data;
error_m=0;
error_count=0;
count_st.hist=size(align2d,1); count_st.step=1; count_st.iteration=1;
num_of_steps=align2d(count_st.hist,1).rec.file.Run.Num_of_steps;
output_dir=align2d(count_st.hist,1).rec.file.Outputdir;
model=tom_emread(align2d(count_st.hist,1).rec.file.Startmodel_Path); model=model.Value;
demo=align2d(count_st.hist,1).rec.file.DemoMode;
use_index_search=align2d(count_st.hist,1).rec.classify.indexs.on;

%build file structure
tom_av2_build_file_struct(output_dir,'rec',num_of_steps);

for i=1:num_of_steps

    iterations=str2num(align2d(count_st.hist,1).rec.file.Run.Iterations_per_step{i});
    
    count_st.step=i; 
    for ii=1:iterations
        count_st.iteration=ii;    
        disp(['2d reconstruction start: step nr' num2str(i) '  iteration nr: ', num2str(ii) , '  start time: ' datestr(clock)]);


        if (error_m==0)
            [align2d error_m]=tom_av2_create_projections(model,align2d,count_st,demo);
        end;

        if (error_m==0)
            if (use_index_search==0)    
                [align2d error_m]=tom_av2_angular_classification(align2d,count_st,demo); % correlate projections and particles
            else
                [align2d error_m]=tom_av2_index_angular_classification(align2d,count_st,demo); % correlate projections and particles ... using index search
            end
        end;

        if (error_m==0)
            [model error_m]=tom_av2_backproj(align2d,count_st,demo);
        end;

        if (error_m==0)
            [model error_m]=tom_av2_post_processing(align2d,count_st,model);
        end;

        if (error_m > 0 )
            switch error_count
                case 2
                    tom_purgequeue(align2d(1,1).classify.paraell.jobmanager);
                    unix('chmod -R ugo+rwx /fs/sally/pool-baumeister/');
                case 5
                    % unix('sudo /etc/init.d/matlab');
                case 100
            end;
            error_m=0;
            error_count=error_count+1;
        end;

        %save(['align/align2d_save'],'align2d');
        count_st.hist=count_st.hist+1;
        align2d(count_st.hist,1)=align2d(count_st.hist-1);
    end;
        
end;











