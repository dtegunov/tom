function stack = tom_os3_alignStack(stack,reference,mask)
    
    filter_param.Apply = 0;
    %transMask = tom_os3_pasteCenter(zeros(size(reference),'single'),tom_os3_sphereMask(zeros(10)));
    
    if(~exist('mask'))
        mask = tom_os3_sphereMask(reference);
    end;
    
%    reference = tom_norm(((reference)),'mean0+1std');
    szt=round(size(reference,1).*0.18);
    transMask=tom_spheremask(ones(size(reference)),szt);
    
    for i= 1:size(stack,3)
        normed_part = tom_norm(stack(:,:,i),'mean0+1std');
        [a b c d]=tom_av2_align(reference,normed_part,ones(size(reference)),[],transMask,filter_param,1,0);
          
        stack(:,:,i)= tom_norm(d,'mean0+1std',mask).*mask;
    end; 
    
    
    
    