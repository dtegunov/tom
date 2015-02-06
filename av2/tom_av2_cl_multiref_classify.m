function [classes cc_out]=tom_av2_cl_multiref_classify(cent_st,stack)
%  TOM_AV2_CL_MULTIREF_CLASSIFY classifies using multiref alignment
%  
%     classes=tom_av2_cl_multiref_classify(cent_st,stack)
%  
%  PARAMETERS
%  
%    INPUT
%     cent_st             structure optained form tom_av2_cl_multiref_train
%     stack               stack of images which should be classified  
%    
%    OUTPUT
%    
%     classes             vector containing the classes
%     cc_out              cross correlation values    
%
%  EXAMPLE
%  
%   
%   [classes cc_out]=tom_av2_cl_multiref_classify(mult_cl,stack);
%   
%  
%
%  
%  REFERENCES
%  
%  SEE ALSO
%     tom_av2_cl_multiref_train
%  
%     created by FB 
%  
%     Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%     Journal of Structural Biology, 149 (2005), 227-234.
%  
%     Copyright (c) 2004-2007
%     TOM toolbox for Electron Tomography
%     Max-Planck-Institute of Biochemistry
%     Dept. Molecular Structural Biology
%     82152 Martinsried, Germany
%     http://www.biochem.mpg.de/tom

if (isfield(cent_st,'mult_alg')==0)
    cent_st.mult_alg=1;
end;
num_cl_good=cent_st.opt_num_cl_good;
num_cl_bad=cent_st.opt_num_cl_bad;

classes=zeros(size(stack,3),1);
v=cat(1,ones(num_cl_good,1),zeros(num_cl_bad,1));
ccc=zeros(size(v,1),1);


if (cent_st.mult_alg==1)
     mask_trans=tom_spheremask(ones(size(cent_st.cents,1),size(cent_st.cents,2)),round((cent_st.alg_max_sh./100).*size(cent_st.cents,1)));
     [stack_out all_tmpl cc_out cc_sum m_pos]=tom_os3_alignStack2(stack,cent_st.cents,'default',[0 0],'mean0+1std','',[1 3],mask_trans);
     classes=(m_pos<=num_cl_good) .* (cc_out > cent_st.cc_opt);  
else
    cc_out=zeros(size(stack,3),1);
    for i=1:size(stack,3)
        for ii=1:length(v)
             ccc(ii)=tom_ccc(stack(:,:,i),cent_st.cents(:,:,ii),'norm');
        end;
        [m_val m_pos]=max(ccc);
        classes(i)=m_pos<=num_cl_good;
        ccc=zeros(size(v,1),1);
        cc_out(i)=m_val;
    end;
    
end;
