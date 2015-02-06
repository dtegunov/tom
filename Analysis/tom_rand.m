function v=tom_rand(max_num,v_length,determ_flag)
%TOM_RAND calculates unique rand numbers 
%
%   img_out=tom_reshape_stack_memory(stack,binning,mask,normstr);
%PARAMETERS
%
%  INPUT
%   max_num             biggest rand num
%   v_length            number of rands
%   determ_flag         flag for determ behaviour (error,dublicate,stop)

%  
%  OUTPUT
%   v                   vector of rand numbers
%
%EXAMPLE
%   img_out=tom_reshape_stack_memory(stack,1);
%   
%  
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by ... fb(eckster) anno 2008
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


for i=1:20
    v_tmp=round(rand(v_length.*4,1)*(max_num-1))+1;
    v_tmp=unique(v_tmp);
    v_tmp=v_tmp(randperm(length(v_tmp)));
    if (length(v_tmp) > v_length)
        v=v_tmp(randperm(v_length));
        break;
    end;
end;