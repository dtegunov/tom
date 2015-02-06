function variance=tom_calc_variance_stack(stack,flag)

if nargin < 2
   flag='use_all';     
end;

if (size(stack,3)==0)
    variance=zeros(size(stack,1),size(stack,2));
    return;
end;
    
    
mean_im=sum(stack,3)./size(stack,3);

variance=zeros(size(stack,1),size(stack,2));
for i=1:size(stack,3)
    if (std2(stack(:,:,i))~=0 || strcmp(flag,'use_not_const')==0 )
        diff=(stack(:,:,i)-mean_im);
        variance=variance + (diff.*diff);
    end;
end;

variance=variance./size(stack,3);

