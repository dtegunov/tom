function classes=tom_os3_classify_templStack(stack,templ,num_of_classes,bin_al,bin_cl,programm,num_of_rot,num_of_al_iter)


if (nargin < 8)
    num_of_al_iter=1;
end;

sz=size(stack);

mask=tom_spheremask(ones(sz(1),sz(2)),round(sz(1)./2)-1,0);
eigenEnd=10;

stack=tom_av2_process_stack(stack,'','norm','mean0+1std');
% parfor i=1:size(stack,3)
%     stack(:,:,i)=tom_norm(tom_xmipp_normalize(stack(:,:,i),'Ramp'),'mean0+1std');
% end;

if (strcmp(programm,'align and classify') || strcmp(programm,'align rotate and classify'))
    if (isempty(templ)==0)
          [stack all_tmpl cc_out]= tom_os3_alignStack2(stack,templ,mask,[bin_al 0],'mean0+1std',[num_of_al_iter 4],'default');
    end;
end;

stack_samp=stack;

if strcmp(programm,'align rotate and classify')
   stack_samp=zeros(size(stack,1),size(stack,2),size(stack,3).*num_of_rot);
    incre=360./num_of_rot;
    zz=1;
    for i=1:size(stack,3)
        ang=round(rand(1).*360);
        for ii=1:num_of_rot
            part=tom_rotate(stack(:,:,i),ang);
            stack_samp(:,:,zz)=part;
            ang=ang+incre;
            zz=zz+1;
        end;
        
    end;
    stack=stack_samp;
end;

stack_org=stack;


stack=stack+0.001*rand(size(stack));
stack=tom_reshape_stack_memory(stack,bin_cl,mask,'mean0+1std');
  
[scores,coefs, eigenvalues]=tom_calc_pca(double(stack'),eigenEnd,'pca','','','','','','',0);

stack=stack_samp;

for i=1:20
    try
        [classifications,centroids,distances,interClusterDistances] =kmeans(scores(1:eigenEnd,:)',num_of_classes);
        break;
    catch Me
        disp(Me.message);
    end;
    if (i>10)
        num_of_classes=num_of_classes-1;
    end;
        
end;

num_per_class=zeros(num_of_classes,1);
classes=zeros(sz,sz,num_of_classes);
for i=1:length(classifications);
    classes(:,:,classifications(i))=classes(:,:,classifications(i))+stack_org(:,:,i);
    num_per_class(classifications(i))=num_per_class(classifications(i))+1;
end;

%norm the classes
for i=1:num_of_classes
     classes(:,:,i)=classes(:,:,i)./num_per_class(i);
end;

