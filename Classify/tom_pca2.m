function tom_pca2(stack_in,pca_st)
%TOM_PCA2 creates ...
%
%   tom_pca2(stack_in,pca_st)
%
%PARAMETERS
%
%  INPUT
%   stack_in            ...
%   pca_st              ...
%  
%  OUTPUT
%
%EXAMPLE
%   tom_pca2(...);
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

error(nargchk(0, 2, nargin, 'struct'))

tom_reshape_stack(stack_in,['outt'],1,0);

r_stack=tom_emread('outt');

r_st_old=r_stack;

[r_stack means]=tom_rm_mean(r_stack.Value);
%r_stack=r_stack.Value;

covM=cov(r_stack);

[V,E]=eig(covM);

%new_data=V(:,1023:1024)'*r_stack';

new_data=V(:,:)'*r_stack';


%data_back=inv(V(:,:)')*new_data;
data_back=V*new_data;

t=tom_rm_mean(data_back',means);
t=reshape(t',32,32,90);
disp('end');