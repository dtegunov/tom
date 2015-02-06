function [ccc_sum_all]=tom_scale3d_model(model,scale_range,stack)




zz=1;
for i=scale_range

    model_sc=tom_rescale3d(model,round((size(model).*i)));
    mod_tmp=ones(size(stack,1),size(stack,1),size(stack,1));
    mod_tmp=mod_tmp.*tom_dev(model_sc,'noinfo');
    
    if (size(model_sc,1) < size(mod_tmp,1))
        corr=round((size(mod_tmp)-size(model_sc))./2);
        mod_tmp=tom_paste2(mod_tmp,model_sc,corr);
    else
        corr=round((size(model_sc)-size(mod_tmp))./2);  
        mod_tmp=tom_cut_out(model_sc,corr,size(mod_tmp),'nofill');
    end;

    proj=sum(mod_tmp,3);

    ccc_sum=0;
    for ii=1:size(stack,3)
        ccf=tom_corr(tom_smooth(proj,size(proj,1)./10),tom_smooth(stack(:,:,ii),size(proj,1)./10),'norm');
        [pos ccc]=tom_peak(ccf);
        ccc_sum=ccc_sum+ccc;
    end;
    ccc_sum_all(zz)=ccc_sum;
    zz=zz+1;
end;