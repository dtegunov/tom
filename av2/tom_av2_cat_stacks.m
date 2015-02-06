function [big]=tom_av2_cat_stacks(stack1,stack2)


big=zeros(size(stack2,1),size(stack1,2),size(stack1,3)+size(stack2,3));


big(:,:,1:size(stack1,3))=stack1;
big(:,:,size(stack1,3)+1:size(stack1,3)+size(stack2,3))=stack2;


