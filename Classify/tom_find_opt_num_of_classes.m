function [dists_a,num_of_classes]=tom_find_opt_num_of_classes(dataset,interv,cluster_function,cluster_params,disp_flag)
%TOM_FIND_OPT_NUM_OF_CLASSES creates ...
%
%   [dists_a,num_of_classes]=tom_find_opt_num_of_classes(dataset,interv,cluster_function,cluster_params,disp_flag)
%
%PARAMETERS
%
%  INPUT
%   dataset             ...
%   interv              ...
%   cluster_function    ...
%   cluster_params      ...
%   disp_flag           ...
%  
%  OUTPUT
%   dists_a       		...
%   num_of_classes		...
%
%EXAMPLE
%   ... = tom_find_opt_num_of_classes(...);
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

error(nargchk(0, 5, nargin, 'struct'))

if nargin < 5
    disp_flag = 0;
end

if disp_flag == 1
    h = waitbar(0,'Clustering...');
end

zz=1;

for i=interv
    
    %cluster it
    if (strcmp(cluster_function,'k-means'))
         [classes centriod sum_dist dists] = kmeans(dataset,i);
    end;
    
    %calc mean dist
    silh=silhouette(dataset,classes);
    
    m_sh=mean(silh);
    
    dists_a(zz,1)=m_sh;
    dists_a(zz,2)=i;
    
    zz=zz+1;
    
    if disp_flag == 1
        waitbar(i./size(interv,2),h,[num2str(i), ' of ', num2str(size(interv,2)), ' combinations done.']);
    end

end;


%opt num of classes

[num,ind]=max(dists_a(:,1));

num_of_classes=dists_a(ind,2);

if disp_flag == 1
    close(h);
end