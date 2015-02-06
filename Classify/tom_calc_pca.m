function [scores,coefs,eigenvalues]=tom_calc_pca(data,num_of_eigs,pca_type,Align,iter_num,binning,wedge_size,mask_st,filter,verbose)
%TOM_CALC_PCA calculates principle component analysis
%
%   [scores,coefs,eigenvalues]=tom_calc_pca(data,num_of_eigs,pca_type,Align,iter_num,binning,wedge_size,mask_st,filter)
%PARAMETERS
%
%  INPUT
%   data                reshaped stack (2d) with tom_reshape_stack
%   num_of_eigs         number of eigenimages (vectors)
%   pca_type            pca or nlca or cpca
%   Align               the align structure
%   iter_num            history number of alignment structure
%   binning             0,1,2 ...
%   wedge_size          size of wedge
%   mask_st             mask structure
%   filter              filter structure
%   verbose             (1) 0 for no output     
%  
%  OUTPUT
%   scores      		projected data
%   coefs       		eigenvectors
%   eigenvalues 		eigenvalues
%
%   check statistics toolbox pca for more infos
%
%EXAMPLE
%   [scores,coefs,eigenvalues]=tom_calc_pca(data,20,'pca','','','','','','',0)
%    calcs pca without output
%  
%  [scores,coefs,eigenvalues]=tom_calc_pca(data,20)
%    calcs pca with Ritz output (Arpack db-level 1)
%
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

error(nargchk(0, 10, nargin, 'struct'))

if (nargin < 2)
   num_of_eigs=20;
end;

if (nargin < 3)
    pca_type='pca';
end;

if (isempty(pca_type))
    pca_type='pca';
end;

if (nargin < 10)
    verbose=1;
end;

db_flag=0;

if (strcmp(pca_type,'pca')==1)
    covM=cov(data);
    options.disp=verbose;
    [coefs,eigenvalues]=eigs(covM,num_of_eigs,'lm',options);
    scores=coefs'*data';
    if (db_flag==1)
        covM_part=cov(data');
        tom_emwrite('cov_part.em',covM_part);
    end;
        
end;

if (strcmp(pca_type,'nlca')==1)
    scores=NLCA(st_new',ones(size(st_new,2),1),[3 900]);
    coefs=0;
end;




if (strcmp(pca_type,'cpca')==1)

    if (nargin==6)
        wedge_size=30;
    end;

    if (nargin==8)
        filter=1;
    end;


    %calculate correlation Matrix
    covM=tom_calc_correlation_matrix(Align,iter_num,binning,'3d',wedge_size,mask_st,filter);
    %covM=tom_calc_correlation_matrix(Align,0,'3d');
    [tmp_coefs,eigenvalues]=eigs(double(covM),num_of_eigs);


    %transform eigs in pix space
    h=tom_reademheader(Align(iter_num,1).Filename);
    coefs=zeros(num_of_eigs,h.Header.Size(1).*h.Header.Size(1).*h.Header.Size(1));

    figure;

    for i=1:num_of_eigs
        for ii=1:size(Align(iter_num,:),2)
            vol_tmp=tom_emread(Align(iter_num,ii).Filename);
            phi=Align(iter_num,ii).Angle.Phi;
            psi=Align(iter_num,ii).Angle.Psi;
            the=Align(iter_num,ii).Angle.Theta;
            %tshift=[Align(iter_num,ii).Shift.X./(2^binning) Align(iter_num,ii).Shift.Y./(2^binning) Align(iter_num,ii).Shift.Z./(2^binning)];
            tshift=[Align(iter_num,ii).Shift.X Align(iter_num,ii).Shift.Y Align(iter_num,ii).Shift.Z];
            vol_tmp= double(tom_rotate(tom_shift(vol_tmp.Value,-tshift),[-psi -phi -the]));
            tom_dspcub(vol_tmp);
            vol_tmp=reshape(vol_tmp,1,h.Header.Size(1).*h.Header.Size(2).*h.Header.Size(3));
            coefs(i,:)=coefs(i,:)+tmp_coefs(ii,i).*vol_tmp;
        end;
    end;

    %project real data in reduced space
    for i=1:num_of_eigs

        for ii=1:size(Align(iter_num,:),2)
            phi=Align(iter_num,ii).Angle.Phi;
            psi=Align(iter_num,ii).Angle.Psi;
            the=Align(iter_num,ii).Angle.Theta;
            %tshift=[Align(iter_num,ii).Shift.X./(2^binning) Align(iter_num,ii).Shift.Y./(2^binning) Align(iter_num,ii).Shift.Z./(2^binning)];
            tshift=[Align(iter_num,ii).Shift.X Align(iter_num,ii).Shift.Y Align(iter_num,ii).Shift.Z];
            vol_tmp=tom_emread(Align(iter_num,ii).Filename);
            vol_tmp= double(tom_rotate(tom_shift(vol_tmp.Value,-tshift),[-psi -phi -the]));
            scores(i,ii)=coefs(i,:)*reshape(vol_tmp,1,h.Header.Size(1).*h.Header.Size(2).*h.Header.Size(3) )';
        end;
    end;

    coefs=coefs';

    disp('done');



end;

