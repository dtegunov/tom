function data_clean=tom_eliminate_outlayers(data,n_std,iterations,verbose_flag)
%TOM_ELIMINATE_OUTLAYERS eliminates outlayers iteratively
%
%   data_clean=tom_eliminate_outlayers(data,n_std,iterations,verbose)
%
%PARAMETERS
%
%  INPUT
%   data                dataset (1D 2D 3D ... and back)
%   n_std               n times standardeviation (2)
%   iterations          iterations (20)    
%   verbose_flag        verbose_flag (0)                        
%
%  OUTPUT
%   data_clean          filtered data            
%
%
%EXAMPLE
% 
% im=rand(100,100); figure; tom_imagesc(im); 
% im2=im; im2(round(rand(100,1).*1000))=1000000; figure; tom_imagesc(im2);
% im_out=tom_eliminate_outlayers(im); figure; tom_imagesc(im_out);
% im_out2=tom_eliminate_outlayers(im2); figure; tom_imagesc(im_out2); 
%
%
%REFERENCES
%
%SEE ALSO
%   -
%
%   
%   fb 03/12/007
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



if (nargin < 2)
    n_std=2;
end;    

if (nargin < 3)
    iterations=20;
end;

if (nargin < 4)
    verbose_flag='hurz';
end;


index=find( (data >= min(min(min(data))))  );

for i=1:iterations
    [mean_out,max_out,min_out,std_out,var_out]=tom_dev(data(index),'noinfo');
    index=find(abs(data-mean_out) < n_std.*std_out); 
    index2=find(abs(data-mean_out) >= n_std.*std_out); 
    if (strcmp(verbose_flag,'noinfo')==0)
        disp(['Iteration: '  num2str(i)  ' mean: ' num2str(mean_out) '  std: ' num2str(std_out) '  number of datapoints: ' num2str(length(index)) ]);
    end;
end;


data_clean=data;
data_clean(index2)=mean_out;


% figure; plot(data); hold on; plot(data_clean,'r-'); hold off;
% disp('end');

