function tom_av2_xmipp_invert_sel(f_sel,norm)
%TOM_AV2_XMIPP_INVERT_SEL inverts particle of sel file
%
%   tom_av2_xmipp_invert_sel(f_sel,repl,repl_to,bin,new_sel)
%
%  TOM_AV2_XMIPP_INVERT_SEL inverts parts in a sel (inplace)
%  
%
%PARAMETERS
%
%  INPUT
%   f_sel       filename of the input sel
%   norm        1: mean0+1std 
%               2: Ramp + mean0+1std        
%
%EXAMPLE
%     
%  tom_av2_xmipp_invert_sel('13_corrf_high_128_clean.sel');
%   
%  tom_av2_xmipp_invert_sel('13_corrf_high_128_clean.sel',2);
%   
%
%REFERENCES
%
%SEE ALSO
%   tom_aling2d
%
%   created by fb ...ole !!
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



d=importdata(f_sel);

if (nargin <2)
    norm=0;
end;
%tic;
parfor i=1:length(d.textdata)
    if(exist(d.textdata{i},'file'))
        try
        f_name=d.textdata{i};
       
        im=tom_spiderread(f_name);    
        im=im.Value;
        if (norm~=3)
            im=-im;
        end;
        if (norm==1)
            im=tom_norm(im,'mean0+1std');
        end;
        
        if (norm>1)
            im=tom_xmipp_normalize(im,'Ramp');
            im=tom_norm(im,'mean0+1std');
        end;
        
        tom_spiderwrite(f_name,im);
        catch ME
            disp([f_name ' not readable skipping !']);
        end;
       
        
    else
        disp([d.textdata(i) 'not found ...skipping']);
    end;

    if (mod(i,1000)==0)
     %   toc;
        disp(num2str(i));
      %  tic;
    end;
    
end;



