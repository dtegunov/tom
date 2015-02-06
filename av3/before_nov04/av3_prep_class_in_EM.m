function [eigstack, evstack] = av3_prep_class_in_EM(ccc,neig,particlefilename,motl,wedgelist,iclass,nx,idspcubmode)
%   [eigstack, evstack] = av3_prep_class_in_EM(ccc,neig,particlefilename,motl,wedgelist,iclass,nx,idspcubmode)
%
%   ccc                 correlation matrix
%   neig                no of eigenvectors
%   particlefilename    expected as <particlefilename>_index.em
%   motl                motive list
%   wedgelist
%   iclass              class of particles - asssigned by some previous
%                       classification. Only these particles are included
%                       into the eigenvectors (eigstack)
%   nx                  number of slices in x for dspcub
%   idspcubmode         dspcubmode - 0,1,2 for xy, xz, and yz slices
%
%   Procedure aims to create stacks that are suitable for subsequent
%   classification using the EM program ('CLASS' command). 
error(nargchk(7,8,nargin))
if nargin < 8
    idspcubmode = 0;
end;
[eigvec eigval] = eigs(ccc,neig);
evstack = zeros(neig+2,size(motl,2)+1);
evstack(1,2:size(motl,2)+1) = 1:size(motl,2);
for ieig=1:neig
    evstack(ieig+1,1) = eigval(ieig,ieig);
end;
totvar=0;
for ii=1:size(ccc,1)
    totvar = totvar+ccc(ii,ii);%trace
end;
evstack(neig+2,1)=totvar;
evstack(2:neig+1,2:size(motl,2)+1)=eigvec';
for ieig=1:neig
    average = av3_eigvec2vol(motl, eigvec(:,ieig), particlefilename,wedgelist,iclass);%rec eigenvector
    im = av3_vol2gallery(average,nx,idspcubmode,ones(size(average)),'nothing');
    if ieig==1
        eigstack = im;
        eigstack(:,:,neig)=0;
    end;
    eigstack(:,:,ieig) = im;
    evstack(ieig+1,2:size(motl,2)+1) = evstack(ieig+1,2:size(motl,2)+1)*sqrt(eigval(ieig,ieig));
    %proper normalization of eigenvectors - eigvec*eigvec' = eigval
end;
