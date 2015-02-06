function cc=tom_check_symmetry(vol,ang_samp,kernel_size)
%TOM_CHECK_SYMMETRY checks symmetry
%
%   cc=tom_check_symmetry(vol,ang_samp,kernel_size)
%
%PARAMETERS
%
%  INPUT
%   vol                 volume in memory
%   ang_samp            angular sampling [start increment stop]
%   kernel_size         (opt.) size of the kernel volume 
%
%  OUTPUT
%
%   cc                  normed cross correlation function
%                         colums == angles
%                         rows == sampling in z 
%
%EXAMPLE
%   cc=tom_check_symmetry(volume,[1 5 355],5);
%   figure; tom_imagesc(cc)
%
%   cc=tom_check_symmetry(volume,[1 2 358]);
%   figure; plot(cc); 
%
%
%
%
%
%REFERENCES
%
%SEE ALSO
%     TOM_SYMREF,TOM_CORR
%
%   created by FB (eckster) 08/09/08
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


%disp('');

 if  nargin < 3
     kernel_size=size(vol,3)-1;
end;

demo=0;

if (demo==1)
    figure;
end;

zz=1;
for ii=1:size(vol,3)-kernel_size
    v_part=vol(:,:,ii:ii+kernel_size);
    
    for i=ang_samp(1):ang_samp(2):ang_samp(3);
        ccf=tom_corr(v_part,tom_rotate(v_part,i),'norm');
        [a cc(ii,zz)]=tom_peak(ccf);
        %cc(ii,zz)=tom_ccc(v_part,tom_rotate(v_part,i),'norm');
       % sublplot
        
        zz=zz+1;
    end;
    zz=1;
    disp(['Slice number: ' num2str(ii) ' finished.']);
end;


