function [eigstack, evstack] = av3_prep_class_in_EM(ccc,neig,particlefilename,...
    motl,wedgelist,iclass,nx,idspcubmode,flag)
% AV3_PREP_CLASS_IN_EM perpares stacks for classification in EM
%
%   [eigstack, evstack] = av3_prep_class_in_EM(ccc,neig,particlefilename,...
%       motl,wedgelist,iclass,nx,idspcubmode,flag)
%   Procedure aims to create stacks that are suitable for subsequent
%   classification using the EM program ('CLASS' command).
%
% PARAMETERS
%  INPUT
%   ccc                 correlation matrix (procedure AV3_CCCMAT)
%   neig                no of eigenvectors
%   particlefilename    expected as <particlefilename>_index.em
%   motl                motive list
%   wedgelist
%   iclass              class of particles - asssigned by some previous
%                       classification. Only these particles are included
%                       into the eigenvectors (eigstack)
%   nx                  number of slices in x for dspcub
%   idspcubmode         dspcubmode - 0,1,2 for xy, xz, and yz slices
%   flag                flag for normalization of eigenvectors: if set to
%                       'EM' eigenvectors (initially normalized to unity)
%                       are multiplied with sqrt(eigval). (default)
%                       Alternatively, the eigenvectors can be multiplied
%                       with eigval (any other flag than 'EM').
%
%  OUTPUT
%   eigstack            stack containing the eigenvectors
%   evstack             stack containg the eigenvalues and -factors
%                       Format:
%                       1st column: particle number
%                       2nd column: eigenfactors 1st EV
%                       neig+1    : eigenfactors of EV neig
%                       neig+2    : reserved for classification
%                       1st row   : eigenvalues
%
% SEE ALSO
%   AV3_CCCMAT
%
% last change 03/31/05 FF - updated docu
error(nargchk(7,9,nargin))
if nargin < 9
    flag = 'EM';
end;
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
    if strmatch(flag,'EM')
        evstack(ieig+1,2:size(motl,2)+1) = evstack(ieig+1,2:size(motl,2)+1)*sqrt(abs(eigval(ieig,ieig)));
        %normalization of eigenvectors according to EM 
        % eigvec*eigvec' =eigval
    else
        evstack(ieig+1,2:size(motl,2)+1) = evstack(ieig+1,2:size(motl,2)+1)*(abs(eigval(ieig,ieig)));
    end;
    
end;
