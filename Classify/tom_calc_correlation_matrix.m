function corr_m=tom_calc_correlation_matrix(Align,iter_num,binning,dim,wedge_size,mask_st,filter)
%TOM_CALC_CORRELATION_MATRIX creates ...
%
%   corr_m=tom_calc_correlation_matrix(Align,iter_num,binning,dim,wedge_size,mask_st,filter)
%
%PARAMETERS
%
%  INPUT
%   Align               ...
%   iter_num            ...
%   binning             ...
%   dim                 ...
%   wedge_size          ...
%   mask_st             ...
%   filter              ...
%  
%  OUTPUT
%   corr_m      		...
%
%EXAMPLE
%   .. = tom_calc_correlation_matrix(...);
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

error(nargchk(0, 7, nargin, 'struct'))



if (strcmp(dim,'3d')==1)
    
    num_of_files=size(Align,2);
    mask = tom_create_mask(mask_st.mask);
%    mask=tom_bin(mask,binning);
    
    %allocate memory
    corr_m=zeros(num_of_files,num_of_files,'single');

    %sum up all rotations
    for i=1:num_of_files
        for iteration=1:iter_num
            angles(iteration,1)=Align(iteration,i).Angle.Phi;
            angles(iteration,2)=Align(iteration,i).Angle.Psi;
            angles(iteration,3)=Align(iteration,i).Angle.Theta;
            shifts(iteration,1)=Align(iteration,i).Shift.X./(2^binning);
            shifts(iteration,2)=Align(iteration,i).Shift.Y./(2^binning);
            shifts(iteration,3)=Align(iteration,i).Shift.Z./(2^binning);
        end;
        [euler_out(i,:) shift_out(i,:) rot_M(i,:,:)]=tom_sum_rotation(-angles,-shifts);
    end;

    
    %calculate Avg
    h=tom_reademheader(Align(iter_num,1).Filename);
    sz=(h.Header.Size')./(2^binning);
    
    
    mean_part=tom_av3_average(Align(1,:));
    
    
    yyy=zeros(h.Header.Size');
    yyy(1,1,1) =1;
    
    
    %calculate correlation  
    for i=1:num_of_files
       
        
        phi=Align(iter_num,i).Angle.Phi;
        psi=Align(iter_num,i).Angle.Psi;
        the=Align(iter_num,i).Angle.Theta;
%        tshift=[Align(iter_num,i).Shift.X./(2^binning) Align(iter_num,i).Shift.Y./(2^binning) Align(iter_num,i).Shift.Z./(2^binning)];
        tshift=[Align(iter_num,i).Shift.X Align(iter_num,i).Shift.Y Align(iter_num,i).Shift.Z];       
        part_A=tom_emread((Align(iter_num,i).Filename));
        
        part_A=tom_bin(part_A.Value,binning).*(mask);
        part_A = double(tom_rotate(tom_shift(part_A,-tshift),[-psi -phi -the]));
        wedge_A_rot=tom_rotate(tom_wedge(ones(size(part_A)),wedge_size),[-psi -phi -the]);
        psf_A=real(tom_ifourier(ifftshift(fftshift(tom_fourier(yyy)).*wedge_A_rot)));
        
        for ii=(i):num_of_files
            
            
            phi=Align(iter_num,ii).Angle.Phi;
            psi=Align(iter_num,ii).Angle.Psi;
            the=Align(iter_num,ii).Angle.Theta;
%            tshift=[Align(iter_num,ii).Shift.X./(2^binning) Align(iter_num,ii).Shift.Y./(2^binning) Align(iter_num,ii).Shift.Z./(2^binning)];
            tshift=[Align(iter_num,ii).Shift.X Align(iter_num,ii).Shift.Y Align(iter_num,ii).Shift.Z];
            part_B=tom_emread((Align(iter_num,ii).Filename));
            part_B=tom_bin(part_B.Value,binning).*mask;
            part_B = double(tom_rotate(tom_shift(part_B,-tshift),[-psi -phi -the]));
            
            wedge_B_rot=tom_rotate(tom_wedge(ones(size(part_B)),wedge_size),[-psi -phi -the]);
            psf_B=real(tom_ifourier(ifftshift(fftshift(tom_fourier(yyy)).*wedge_B_rot)));
            
            
            ccf=tom_corr(part_A,part_B,'norm',mask,psf_A,psf_B);
            
            %ccf=tom_corr(part_A,part_B,'norm',mask);
            [pos corr_m(i,ii)]=tom_peak(ccf);
            disp([num2str(i) ' ' num2str(ii) ' ccf:' num2str(corr_m(i,ii))]);
        end;
    end;
    
    
    %%fill up symmetric values
    for i=2:num_of_files
         for ii=1:(i-1)
            corr_m(i,ii)=corr_m(ii,i);
         end;
    end;
    
   
    
    
    else

    %to do !    
    

end;