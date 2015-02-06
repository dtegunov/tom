function [covM use_parts cc_cov cc_eig eigen_vals]=tom_cov_incremental(outfilestruct,num_of_iterations,num_of_eigs,start_num,multiplier,binning,mask,normstr,output)
%TOM_COV_INCREMENTAL calculates incremental cov and eigs ...to examine
%significant sample
%
%   [covM part_nr cc cc_eig eigen_vals]=tom_cov_incremental(outfilestruct,num_of_iterations,num_of_eigs,start_num,multiplier,binning,mask,output)
%PARAMETERS
%
%  INPUT
%   outfilestruct           reshaped stack (2d) with tom_reshape_stack
%   num_of_iterations       number of eigenimages (vectors)
%   num_of_eigs             pca or nlca or cpca
%   start_num               the align structure
%   multiplier              history number of alignment structure
%   binning                 0,1,2 ...
%   mask                    mask structure
%   normstr                 'mean0+1std' '3std' ...  or no_norm to switch off
%   output                   0 no output 1 text 2 gui 
%
%  
%  OUTPUT
%   covM      		projected data
%   use_parts       particle number per iteration
%   cc_cov   		correlation coefficien of cov(n) with cov(n-1)
%   cc_eig          correlation coefficien of  eigv(n) with eigv(n-1)   
%   eigen_vals      eigenvalues 
%
%   check statistics toolbox pca for more infos
%
%EXAMPLE
%   [covM part_nr cc_cov cc_eig
%   eigen_vals]=tom_cov_incremental(outfilestruct,20,10,10,1.5,1,ones(80,80),'no_norm',2);
%      
%  
%
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by fb(eckster)  2008
%   
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
    num_of_iterations=20;
end;

if (nargin < 3)
    num_of_eigs=10;
end;

if (nargin < 4)
    start_num=num_of_eigs;
end;

if (nargin < 5)
    multiplier=1.5;
end;

if (nargin < 6)
    binning=1;
end;

if (nargin < 7)
    mask=ones(outfilestruct.size(1),outfilestruct.size(2));
end;

if (nargin < 8)
    normstr='no_norm';
end;

if (nargin < 9)
    output=1;
end;

if (output > 1)
    figure;
end;

num_of_parts=outfilestruct.idx.num_of_files;

part_nr(1)=0;
cc_cov(1)=0;
cc_eig(1,:)=zeros(num_of_eigs,1);
eigen_vals(1,:)=zeros(num_of_eigs,1);
use_parts(1)=start_num;
zz=1;

for i=1:num_of_iterations
    
    part_nr=round(rand(use_parts(i),1).*num_of_parts);
    real_nr=outfilestruct.idx.written_index(part_nr);
    
    stack=tom_av2_stack_read(outfilestruct,real_nr);
    
    img_out=tom_reshape_stack_memory(stack,binning,mask,normstr);
    
    covM=cov(img_out');
    [scores,coefs,eigenvalues]=tom_calc_pca(img_out,num_of_eigs,'pca','','','','','','',0);
    eigen_vals(i,:)=(eigenvalues*ones(num_of_eigs,1));

    if (i>1)
        [cc_cov(i) cc_eig(i,:)]=calc_eigenval_correlation(scores,scores_old,covM,covM_old);
    end;
    covM_old=covM;
    scores_old=scores;
    
    if (output > 0)
        disp(['iteration: ' num2str(i)  ' cc cov: ' num2str(cc_cov(i)) ' used particles: ' num2str(use_parts(i)) ]);
        disp(['  cc eigs: ' num2str(cc_eig(i,:))  ]);
        disp(['  eigvals: ' num2str(eigen_vals(i,:))  ]);
        disp(' ');
    end;
    
    if (output > 1)
        zz=display_eigs(scores,covM,zz,i);
    end;
    
    use_parts(i+1)=round(use_parts(i).*multiplier);
    
 end;

 use_parts=use_parts(1:i);


 function [cc_cov cc_eig]=calc_eigenval_correlation(scores,scores_old,covM,covM_old)

sz=round(sqrt(size(scores(1,:),2)));

for i=1:size(scores,1)
    tmp=reshape(scores(i,:),sz,sz);
    tmp_old=reshape(scores_old(i,:),sz,sz); 
    cc_eig(i)=tom_ccc(tmp,tmp_old,'norm');
    if cc_eig(i) < 0
         cc_eig(i)=tom_ccc(tmp,-tmp_old,'norm');
    end;
end;

cc_cov=tom_ccc(covM,covM_old,'norm');


function zz=display_eigs(scores,covM,zz,i)

sz=round(sqrt(size(scores(1,:),2)));

subplot(5,5,zz); imagesc(covM); axis image;  title([num2str(i)]); colormap hot;  drawnow; zz=zz+1;
subplot(5,5,zz); imagesc(reshape(scores(1,:),sz,sz)); axis image;  title([num2str(i)]); colormap hot;  drawnow; zz=zz+1;
subplot(5,5,zz); imagesc(reshape(scores(2,:),sz,sz)); axis image;  title([num2str(i)]); colormap gray;  drawnow; zz=zz+1;
subplot(5,5,zz); imagesc(reshape(scores(3,:),sz,sz)); axis image;  title([num2str(i)]); colormap gray;  drawnow; zz=zz+1;
subplot(5,5,zz); imagesc(reshape(scores(4,:),sz,sz)); axis image;  title([num2str(i)]); colormap gray;  drawnow; zz=zz+1;

if (mod(i,5)==0)
    zz=1;
end;













