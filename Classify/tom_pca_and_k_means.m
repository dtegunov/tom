function [idx,centriod,sum_dist,dists]=tom_pca_and_k_means(stack,num_of_classes,eigs_start_stop,binning,equal_flag,verbose)
%TOM_PCA_AND_K_MEANS calculates pca and kmeans 
%
%   [idx,c,sumd]=tom_pca_and_k_means(stack,num_of_classes,eigs_start_stop,binning,equal_flag,demo)
%PARAMETERS
%
%  INPUT
%   stack                particle stack
%   num_of_classes       (5) pca or nlca or cpca
%   eigs_start_stop      (1,20)  vector which determines the first and last eigenvalue 
%   binning              (0) the align structure
%   equal_flag           (0)  0,1 distributes equal to centroids of k-means
%                         (works only for 2 classes)  
%   verbose             (0) output on/off 1 2 Arpack output levels
%  OUTPUT
%   idx                 vector of classes
%   c                   centroids
%   sumd                point-to-centroid distances in the 1-by-K vector
%   
%check statistics toolbox pca for more infos
%
%EXAMPLE
%   [idx,c,sumd]=tom_pca_and_k_means(stack,3,[1 20],2)
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   tom_calc_pca tom_pcagui
%
%   created by fb (eckster)
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

if nargin < 2
    num_of_classes=5;
end;

if nargin < 3
   eigs_start_stop=[1 20];
end;

if nargin < 4
   binning=0;
end;

if nargin < 5
    equal_flag=0;
end;

if nargin < 6
    verbose=0;
end;


normflag=1;
sz=size(stack);
if ischar(stack)
    dimension=3;
else
    dimension=2;
end;

if (dimension==2)
    
    reshapedStack=zeros(sz(3),(sz(1)./2^binning).*(sz(2)./2^binning) );
    for i=1:sz(3)
        im = stack(:,:,i);
        im = tom_bin(im,binning); 
        im = reshape(im,1,[]);
        if (normflag == 1)
            im = tom_norm((im+100).*2,'phase');
        end
        reshapedStack(i,:)=im;
    end;

else
    %3d
    disp('not implemented!');
end;


%calc pca
if (eigs_start_stop(2) > 0)
    [scores,coefs,eigenValues]=tom_calc_pca(double(reshapedStack),eigs_start_stop(2),'pca','','','','','','',verbose);
else
    scores=reshapedStack;
end;


[idx centriod sum_dist dists] = kmeans(scores',num_of_classes);


if (equal_flag==1) && (num_of_classes==2)
    
    
    if length(find(idx==1)) > length(find(idx==2))
        diff=round((length(find(idx==1))-length(find(idx==2)))./2);
        for i=1:diff
            [a idx_min]=max(dists(:,1));
            dists(idx_min,1)=min(dists(:,1))-1;
            idx(idx_min)=2;
        end;
        return;
    end;

    if length(find(idx==1)) < length(find(idx==2))
        diff=round((length(find(idx==2))-length(find(idx==1)))./2);
        for i=1:diff
            [a idx_min]=min(dists(:,2));
            dists(idx_min,2)=min(dists(:,2))-1;
            idx(idx_min)=1;
        end;
        return;
    end;
        
    
    
%     for i=1:num_of_classes
%         
%     end;
    
    
    
    
end;


% cl1=zeros(size(stack(:,:,1)));
% cl2=zeros(size(stack(:,:,1)));

%debug

% for i=1:size(idx,1)
%     if idx(i)==1
%         cl1=cl1+stack(:,:,i);
%     else
%         cl2=cl2+stack(:,:,i);
%     end;
% end;
