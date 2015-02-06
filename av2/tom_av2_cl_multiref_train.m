function train_st=tom_av2_cl_multiref_train(stack_good,stack_bad,options,verbose)
%  TOM_AV2_CL_MULTIREF_TRAIN trains a multiref classifier 
%  
%     train_st=tom_av2_cl_multiref_train(stack_good,stack_bad,options)
%  
%  PARAMETERS
%  
%    INPUT
%     stack_good          stack containing the good images
%     stack_bad           stack containing the bad images
%     options             options struct
%     verbose             (1) verbose flag use 0 2 switch off
%    
%    OUTPUT
%    
%     train_st            struct containing the trained classifier
%
%
%  EXAMPLE
%  
%  options.nr_good_train=280; 
%  options.nr_bad_train=1120;
%  options.eigs=[1 20];
%  options.bin=1;
%  options.num_cl_good=3;
%  options.num_cl_bad=3;
%  options.num_al=3;
%  options.cc_search_range=[0.22:0.02:0.42];
%  options.test_good_bad_ratio=0.25;
%  options.test_max_num=2000;
%
%  train_st=tom_av2_cl_multiref_train(stack_good,stack_bad,options,verbose);
%  
%
%  
%  REFERENCES
%  
%  SEE ALSO
%     tom_av2_cl_multiref_classify
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

if (nargin < 3)
    nr_good_train=[100 200 400];
    nr_bad_train=[400 800 1600];
    cent_st.eigs=[1 20];
    cent_st.num_cl_good=3;
    cent_st.num_cl_bad=3;
    cent_st.num_al=5;
    cent_st.cc_search_range=0.24:0.02:0.42;
    test_good_bad_ratio=0.25;
    
else
    nr_good_train=options.nr_good_train;
    nr_bad_train=options.nr_bad_train;
    cent_st.eigs=options.eigs;
    cent_st.bin=options.bin;
    cent_st.num_cl_good=options.num_cl_good;
    cent_st.num_cl_bad=options.num_cl_bad;
    cent_st.num_al=options.num_al;
    cent_st.cc_search_range=options.cc_search_range;
    test_good_bad_ratio=options.test_good_bad_ratio;
    test_max_num=options.test_max_num;
end;

if (nargin < 4)
    verbos=0;
end;


for i=1:size(stack_good,3);
    stack_good(:,:,i)=tom_filter(stack_good(:,:,i),4);
end;

for i=1:size(stack_bad,3);
    stack_bad(:,:,i)=tom_filter(stack_bad(:,:,i),4);
end;



tmp_ind_good=randperm(size(stack_good,3));
tmp_ind_bad=randperm(size(stack_bad,3));
%tmp_ind_good=1:size(stack_good,3);
%tmp_ind_bad=1:size(stack_bad,3);

max_num_good_test=test_max_num.*test_good_bad_ratio;
max_num_bad_test=test_max_num.*(1-test_good_bad_ratio);



sz_good_test=length(tmp_ind_good)-max(nr_good_train);
sz_bad_test=length(tmp_ind_bad)-max(nr_bad_train);

if (sz_bad_test > max_num_bad_test)
    sz_bad_test=max_num_bad_test;
end

if (sz_good_test > max_num_good_test)
    sz_good_test=max_num_good_test;
end

% if (sz_bad_test*test_good_bad_ratio > sz_good_test)
%     sz_bad_test=round((1./test_good_bad_ratio).*sz_good_test);  
% else
%     sz_good_test=round(test_good_bad_ratio.*sz_bad_test);  
% end;


if (sz_bad_test ~= ((sz_good_test/test_good_bad_ratio)-sz_good_test) )
    if ( ((sz_bad_test*test_good_bad_ratio)/(1-test_good_bad_ratio)) > sz_good_test)
        sz_bad_test=round((sz_good_test/test_good_bad_ratio)-sz_good_test);
    else
        %round(test_good_bad_ratio.*sz_bad_test);
        sz_good_test=round(((sz_bad_test*test_good_bad_ratio)/(1-test_good_bad_ratio)));
    end;
end;



tmp_test_good=tmp_ind_good(end-sz_good_test+1:end);
% if (length(tmp_test_good) > round(test_max_num*test_good_bad_ratio))
%     tmp_test_good=tmp_test_good(1:round(test_max_num*test_good_bad_ratio));
%     
% end;

tmp_test_bad=tmp_ind_bad(end-sz_bad_test+1:end);
% if (length(tmp_test_bad) > round(test_max_num*(1-test_good_bad_ratio)))
%     tmp_test_bad=tmp_test_bad(1:round(test_max_num*(1-test_good_bad_ratio)));
% end;

good_test=stack_good(:,:,tmp_test_good);
bad_test=stack_bad(:,:,tmp_test_bad);

data_stack_test=cat(3,good_test,bad_test);
groups_test=cat(1,ones(sz_good_test,1),zeros(sz_bad_test,1));

disp(['Nr. good 4 test: ' num2str(size(good_test,3)) ]);
disp(['Nr. bad 4 test: ' num2str(size(bad_test,3)) ]);

if (sz_good_test < 0 || sz_bad_test < 0)
    error('input stack smaller tan number of particles requested in training');
end;

for c_part_nr=1:length(nr_good_train)
    
    
    
    sz_good=nr_good_train(c_part_nr);
    sz_bad=nr_bad_train(c_part_nr);
   
    
    
    good_train=stack_good(:,:,tmp_ind_good(1:sz_good));
    bad_train=stack_bad(:,:,tmp_ind_bad(1:sz_bad));
    
    data_stack_train=cat(3,good_train,bad_train);
    groups_train=cat(1,ones(sz_good,1),zeros(sz_bad,1));
    
    disp(['Nr. good 4 train: ' num2str(sz_good)]);
    disp(['Nr. bad 4 train: ' num2str(sz_bad)]);
    
    [all_cent_st{c_part_nr},cl_rate(c_part_nr)]=cent_train(data_stack_train,groups_train,data_stack_test,groups_test,cent_st);
    
    disp(' ');
    
    
end;

[val pos]=min(cl_rate);

train_st=all_cent_st{pos};

train_st.test_good_bad_ratio=options.test_good_bad_ratio;



function [cent_st final_cl_rate]=cent_train(stack,groups,data_stack_test,groups_test,cent_st)

if (nargin<3)
    eigs=[1 20];
    bin=1;
    num_cl_good=3;
    num_cl_bad=3;
    num_al=4;
    cc_search_range=[0.2:0.02:0.46];
else
    eigs=cent_st.eigs;
    bin=cent_st.bin;
    num_al=cent_st.num_al;
    cc_search_range=cent_st.cc_search_range;
end;


for c_num_cl=1:length(cent_st.num_cl_good)
   
    num_cl_good=cent_st.num_cl_good(c_num_cl);
    num_cl_bad=cent_st.num_cl_bad(c_num_cl);
    
    cent_st.opt_num_cl_good=num_cl_good;
    cent_st.opt_num_cl_bad=num_cl_bad;
    
%     idx_all_good=find(groups);
%     idx_all_bad=find(groups==0);
%     
%     l4mult_good=round(idx_all_good).*(1-percent_opt_val));
%     l4mult_bad=round(length(idx_all_bad).*(1-percent_opt_val));
    
    idx_good=find(groups);
    idx_bad=find(groups==0);
    
    idx_test_good=find(groups_test);
    idx_test_bad=find(groups_test==0);
    
    
    nr_good_test=length(idx_test_good);
    
    
    disp(['nr cl good: ' num2str(num_cl_good)]);
    disp(['nr cl bad: ' num2str(num_cl_bad)]);
    
    %data_stack_test=stack(:,:,cat(1,idx_test_good,idx_test_bad));
  %  groups_test=cat(1,ones(length(idx_test_good),1),zeros(length(idx_test_bad),1));
    
   % mask=tom_spheremask(ones(64,64),5);
    
    if (num_cl_good>1)
        %tmp_stack=tom_os3_alignStack2(stack(:,:,idx_good),stack(:,:,idx_good(1)),'default',[0 0],'mean0+1std',[num_al+1 2],mask);
        tmp_stack=stack(:,:,idx_good);
        [idx,c,sumd]=tom_pca_and_k_means(tmp_stack,num_cl_good,eigs,1);
        u_idx=unique(idx);
        for ii=1:length(u_idx)
            cents_good(:,:,ii)=sum(tmp_stack(:,:,find(idx==u_idx(ii))),3)./length(find(idx==u_idx(ii)));
        end;
        [stack_out all_tmpl cc_out cc_sum]=tom_os3_alignStack2(stack(:,:,idx_good),cents_good,'default',[0 0],'mean0+1std',[num_al 1],'default');
        cents_good=all_tmpl{num_al};
    else
        tmp_stack=tom_os3_alignStack2(stack(:,:,idx_good),stack(:,:,idx_good(1)),'default',[0 0],'mean0+1std',[num_al 2]);
        cents_good=sum(tmp_stack,3);
    end;
    
    if (num_cl_bad>1)
        %tmp_stack=tom_os3_alignStack2(stack(:,:,idx_bad),stack(:,:,idx_bad(1)),'default',[0 0],'mean0+1std',[num_al+1 2],mask);
        tmp_stack=stack(:,:,idx_bad);
        [idx,c,sumd]=tom_pca_and_k_means(tmp_stack,num_cl_bad,eigs,1);
        u_idx=unique(idx);
        for ii=1:length(u_idx)
            cents_bad(:,:,ii)=sum(tmp_stack(:,:,find(idx==u_idx(ii))),3)./length(find(idx==u_idx(ii)));
        end;
        [stack_out all_tmpl cc_out cc_sum ref]=tom_os3_alignStack2(stack(:,:,idx_bad),cents_bad,'default',[0 0],'mean0+1std',[num_al 1],'default');
        cents_bad=all_tmpl{num_al};
    else
        tmp_stack=tom_os3_alignStack2(stack(:,:,idx_bad),stack(:,:,idx_bad(1)),'default',[0 0],'mean0+1std',[num_al 2]);
        cents_bad=sum(tmp_stack,3);
    end;
    
 
    cent_st.cents=cat(3,cents_good,cents_bad);
    
    disp('optimizing cc-thresh');
    
    cl_rate=zeros(length(cc_search_range),1);
    fp=zeros(length(cc_search_range),1);
    fn=zeros(length(cc_search_range),1);
    overlap=zeros(length(cc_search_range),1);
    zz=0;
    for i=cc_search_range
        zz=zz+1;
        cent_st.cc_opt=i;
        classes_test=tom_av2_cl_multiref_classify(cent_st,data_stack_test);
        v=(groups_test==classes_test);
        v=v(1:nr_good_test);
        overlap(zz)=length(find(v))./nr_good_test;
        fp(zz)=(length(find((classes_test==1)))-length(find(v)))./nr_good_test;
        fn(zz)=1-overlap(zz);
        cl_rate(zz)=fp(zz)+fn(zz);
        disp(['opt step ' num2str(zz) ' of ' num2str(length(cc_search_range)) ' done']);
        disp(['cc-thresh: ' num2str(cent_st.cc_opt) ' cl-rate: ' num2str(cl_rate(zz))]);
    end;
    
    disp(' ');
    disp(num2str(cc_search_range));
    disp(num2str(cl_rate'));
    disp(' ');
    
    
    [val pos]=min(cl_rate);
    pos=min(pos);
    
    disp(['opt found: ']);
    disp(['overlap is: ' num2str(overlap(pos))]);
    disp(['false positive is: ' num2str(fp(pos))]);
    disp(['false negative is: ' num2str(fn(pos))]);
    disp(['cl_rate ' num2str(fn(pos)+fp(pos))]);
    disp(['cc thresh ' num2str(cc_search_range(pos)) ]);
    
    cent_st.cc_opt=cc_search_range(pos);
    all_cl_rate(c_num_cl)=fn(pos)+fp(pos);
    
    all_cent_st{c_num_cl}=cent_st;
    
    save('cents.mat','cent_st');
    
end;

[val pos]=min(all_cl_rate);
 
final_cl_rate=val;
cent_st=all_cent_st{pos};



function classes=cent_classify(cent_st,stack,num_cl_good,num_cl_bad)


if (isfield(cent_st,'mult_alg')==0)
    cent_st.mult_alg=1;
end;
classes=zeros(size(stack,3),1);
v=cat(1,ones(num_cl_good,1),zeros(num_cl_bad,1));
ccc=zeros(size(v,1),1);


if (cent_st.mult_alg==1)
     [stack_out all_tmpl cc_out cc_sum m_pos]=tom_os3_alignStack2(stack,cent_st.cents,'default',[0 0],'mean0+1std',[1 2],'default');
     classes=(m_pos<=num_cl_good) .* (cc_out > cent_st.thresh);  
else
    
    for i=1:size(stack,3)
        
        for ii=1:length(v)
             ccc(ii)=tom_ccc(stack(:,:,i),cent_st.cents(:,:,ii),'norm');
        end;
        
        [m_val m_pos]=max(ccc);
        classes(i)=m_pos<=num_cl_good;
        ccc=zeros(size(v,1),1);
    end;
    

end;



%  if (isempty(find(classes_test)) )
%         disp('No good particles in test set!');
%     else
%          figure; tom_dspcub(data_stack_test(:,:,find(classes_test) )); set(gcf,'Name','good');
%        % figure; tom_imagesc(sum(data_stack_test(:,:,find(classes_test)) ,3)); set(gcf,'Name','good')
%     end
%     pause(1);
%     if (isempty(find(classes_test==0)) )
%         disp('No bad particles in test set!');
%     else
%         figure; tom_dspcub(data_stack_test(:,:,find(classes_test==0) )); set(gcf,'Name','bad');
%         %figure; tom_imagesc(sum(data_stack_test(:,:,find(classes_test==0)) ,3)); set(gcf,'Name','bad')
%     end;

%     me=mean(cc_out(:,num_al));
%     st=std(cc_out(:,num_al));
%     idx_bg=find(cc_out(:,num_al) > (me-st));
%     tmp=zeros(64,64,length(unique(ref)));
%     for iiii=1:length(idx_bg)
%         tmp(:,:,ref(idx_bg(iiii)))=tmp(:,:,ref(idx_bg(iiii)))+stack_out(:,:,idx_bg(iiii));
%     end;
%     all_tmpl{num_al}=tmp;



% if (size(cents_good,3)==1)
%     figure; tom_imagesc(cents_good); set(gcf,'Name','good');
% else
%     figure; tom_dspcub(cents_good); set(gcf,'Name','good');
% end;
%     
% pause(0.8);
% if (size(cents_good,3)==1)
%     figure; tom_imagesc(cents_bad); set(gcf,'Name','bad');
% else
%     figure; tom_dspcub(cents_bad); set(gcf,'Name','bad');
% end;
% pause(0.8);
% 
% 



% disp(['reshaphing test stack good']);
% good_test_rs=tom_reshape_stack_memory(good_test)';
% disp(['reshaphing test stack bad']);
% bad_test_rs=tom_reshape_stack_memory(bad_test)';
% disp('reshaping done');
% data_test=cat(1,good_test_rs,bad_test_rs);



