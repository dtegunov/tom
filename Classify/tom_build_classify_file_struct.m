function tom_build_classify_file_struct(basedir,run,num_of_classes,dim_flag)
%TOM_BUILD_CLASSIFY_FILE_STRUCT creates ...
%
%   tom_build_classify_file_struct(basedir,run,num_of_classes,dim_flag)
%
%PARAMETERS
%
%  INPUT
%   basedir             ...
%   run                 ...
%   num_of_classes      ...
%   dim_flag            ...
%  
%  OUTPUT
%
%EXAMPLE
%   tom_build_classify_file_struct(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by ... (author date)
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

error(nargchk(0, 4, nargin, 'struct'))

if (isunix==1)
    if (exist([basedir '/run_' num2str(run)],'dir')==0)  
        mkdir([basedir '/run_' num2str(run)]);
    end;
else
    if (exist([basedir '\run_' num2str(run)],'dir')==0)  
        mkdir([basedir '\run_' num2str(run)]);
    end;

end;

for i=1:num_of_classes
    
    if (isunix==1)
         if (exist([basedir '/run_' num2str(run) '/class_' num2str(i) ],'dir')==0)
            mkdir([basedir '/run_' num2str(run) '/class_' num2str(i) ]);
        end;  
    else
        if (exist([basedir '\run_' num2str(run) '\class_' num2str(i) ],'dir')==0)
            mkdir([basedir '\run_' num2str(run) '\class_' num2str(i) ]);
        end;
    end;

end;

if (strcmp(dim_flag,'3d')==1)
    if (exist([basedir '\run_' num2str(run) '\avg' ])==0)
        mkdir([basedir '\run_' num2str(run) '\avg' ])
    end;
end;
