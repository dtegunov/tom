function tom_classify(stack_in,output_dir,align_st,align_out,pre_processing,classifyer,binning,cut,Filter,Mask,Norm)
%TOM_CLASSIFY creates ...
%
%   tom_classify(stack_in,output_dir,align_st,align_out,pre_processing,classifyer,binning,cut,Filter,Mask,Norm)
%
%PARAMETERS
%
%  INPUT
%   stack_in            ...
%   output_dir          ...
%   align_st            ...
%   align_out           ...
%   pre_processing      ...
%   classifyer          ...
%   binning             ...
%   cut                 ...
%   Filter              ...
%   Mask                ...
%   Norm                ...
%  
%  OUTPUT
%
%EXAMPLE
%   tom_classify(...);
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

error(nargchk(0, 11, nargin, 'struct'))

if (isdir(stack_in)==0)
    h=tom_reademheader(stack_in);
    sz_org=h.Header.Size;
    dim_flag='2d';
else
    d=dir(stack_in);
    h=tom_reademheader([stack_in d(3).name]);
    %sz_org=h.Header.Size./(2.^binning);
    sz_org=h.Header.Size;
    dim_flag='3d';
end;

%parse inputs
if (isempty(align_st)==1)
    align_st(1,1).auto_class=-1;
else
    load(align_st);
    align_st=align2d;
    clear align2d;
end;

%transfer often used
iter_num=size(align_st,1);
num_of_classes=classifyer.num_of_classes;


%process stack
%st=tom_emread(stack_in)

tom_build_classify_file_struct(output_dir,iter_num,num_of_classes,dim_flag);

disp('reshaping...');
tom_reshape_stack(stack_in,[output_dir '/run_' num2str(iter_num)  '/r_stack.em'],1,binning);
tom_reshape_stack(stack_in,[output_dir '/run_' num2str(iter_num)  '/r_stack_org.em'],1,binning);

%preprocessing
if isempty(pre_processing)==0
    if (strcmp(pre_processing.Type,'pca'))
        disp('pca...'); %perform pca
        st=tom_emread([output_dir '/run_' num2str(iter_num) '/r_stack.em']);
        [coefs,scores,variances,t2] = princomp(st.Value);
        %reduce the number of dimensions
        sc_tmp=scores(:,1:pre_processing.num_of_eigenvectors);
        tom_emwrite([output_dir '/run_' num2str(iter_num) '/r_stack.em'],sc_tmp);
        tom_emwrite([output_dir '/run_' num2str(iter_num) '/r_coefs.em'],coefs);
        tom_emwrite([output_dir '/run_' num2str(iter_num) '/r_scores.em'],scores);

    end;
end;


%classifying
if (strcmp(classifyer.Type,'k-means')==1)

    %to be replaced by a paraell approach
    st=tom_emread([output_dir '/run_' num2str(iter_num) '/r_stack.em']);
    st=st.Value;

    disp('clustering...');
    [classes centriod sum_dist dists] = kmeans(st,num_of_classes);

end;

%bookkeeping
num=get_class_numbers(classes,num_of_classes);
class_count=zeros(num_of_classes,1);

if (strcmp(dim_flag,'2d'))
    avg=zeros(sz_org(1),sz_org(2),num_of_classes);
    for i=1:num_of_classes
        dir_out=[output_dir '/run_' num2str(iter_num) '/class_' num2str(i) '/stack.em' ];
        tom_emwritec(dir_out,[sz_org(1) sz_org(2) num(i)],'new');
    end;
else
    d=dir(stack_in);
    avg=zeros(sz_org(1),sz_org(2),sz_org(3),num_of_classes);
end;



for i=1:size(classes,1)

    if (strcmp(dim_flag,'2d'))
        im=tom_emread([stack_in],'subregion',[1 1 i],[sz_org(1)-1 sz_org(2)-1 0]);
    else
        im=tom_emread([stack_in '/' d(i+2).name]);
    end;


    im=im.Value;
    %update align struct
    class_count(classes(i))=class_count(classes(i))+1;
    align_st(iter_num+1,i).auto_class=classes(i);

    if (strcmp(dim_flag,'2d'))
        dir_out=[output_dir '/run_' num2str(iter_num) '/class_' num2str(classes(i)) '/stack.em' ];
        tom_emwritec(dir_out,im,'subregion',[1 1 class_count(classes(i))],[sz_org(1) sz_org(2) 1]);
        avg(:,:,classes(i))=avg(:,:,classes(i))+im;
    else
        dir_out=[output_dir '/run_' num2str(iter_num) '/class_' num2str(classes(i)) '/part_' num2str(class_count(classes(i))) '.em' ];
        tom_emwrite(dir_out,im);
        avg(:,:,:,classes(i))=avg(:,:,:,classes(i))+im;
    end;

    class_statistic(classes(i)).members(class_count(classes(i)))=i;
end;


for i=1:num_of_classes
    class_statistic(i).class_dist=sum_dist(i);
    class_statistic(i).num_of_members=num(i);
    class_statistic(i).centriod=centriod(i);
    dir_out=[output_dir '/run_' num2str(iter_num) '/class_' num2str(i) '/avg.em' ];

    if (strcmp(dim_flag,'2d')==1)
        tom_emwrite(dir_out,avg(:,:,i));
    else
        tom_emwrite(dir_out,avg(:,:,:,i));
    end;


    for ii=1:size(class_statistic(i).members,2)
        align_tmp(1,ii)=align_st(iter_num,class_statistic(i).members(ii));
        align_tmp(1,ii).auto_class=0;
    end;
    dir_out=[output_dir '/run_' num2str(iter_num) '/class_' num2str(i) '/stack.mat' ];
    save(dir_out,'align_tmp');
end;

if (strcmp(dim_flag,'2d'))
    tom_emwrite([output_dir '/run_' num2str(iter_num) '/class_avg.em'],avg);
else
    for i=1:num_of_classes
        tom_emwrite([output_dir '/run_' num2str(iter_num) '/avg/avg_' num2str(i) '.em'],avg(:,:,:,i));
    end;
end;

save([output_dir '/run_' num2str(iter_num) '/' align_out],'align_st');
save([output_dir '/run_' num2str(iter_num) '/class_statistic.mat'],'class_statistic');




%validation


function num=get_class_numbers(classes,num_of_classes)

num=zeros(num_of_classes,1);

for i=1:size(classes,1)
    ref_nr=classes(i);
    num(ref_nr)=num(ref_nr)+1;
end;







