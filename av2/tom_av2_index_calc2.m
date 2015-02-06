function tom_av2_index_calc2(filename_stack,filename_index_stack,mask_st,tree_split,tmp_filename,verbose)
%TOM_AV2_INDEX_CALC calculates a index stack for a given image stack ...to
%                   speed up multi-ref alignment
%
%   tom_av2_index_calc(filename_stack,filename_index_stack,flag,tmp_filename,verbose,lookup)
%PARAMETERS
%
%  INPUT
%   filename_stack       filename of input stack
%   filename_index_stack filename (output) of index stack
%   mast_st              mask structure 
%   tree_split             (classify) structure determines how the original stack is devided '2_halfs' or 'classify'
%                          tree_split.Method= '2_halfs' or 'classify
%                          tree_split.Values.Binning  only needed if Method is classify 
%                          tree_split.Values.Eignevalues only needed Method is classify determines 
%                          which Eigenvalues are used to project the data 
%
%
%   tmp_filename         (./idx_dump+datestr(now)) filename for a tmp directory  
%   verbose              (0) flag for output
%   
%  OUTPUT
%   -                  
%
%EXAMPLE
%     
%tree_split.Method='classify'; tree_split.Values.Binnig=1;
%tree_split.Values.Eigenvalues=[1 10];
%
%tom_av2_index_calc('stack_ref.em','stack_ref_index.em',tree_split,'data_av2')
%   
%
%REFERENCES
%
%SEE ALSO
%   tom_av2_index_search, tom_av2_index_bintree_not_2_index.m,
%   tom_av2_index_plot_indexstack.m
%
%   created by fb(eckster)
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

if nargin < 2
    rest=strtok(filename_stack,'.');
    filename_index_stack=[rest '_index.em'];
end;

if nargin < 3
    mask_st.mask.Apply=0;
end;

if nargin < 4
    tree_split.Method='classify';
    tree_split.Values.Binning=1;
    tree_split.Values.Eigenvalues=[0 10];
end;

if nargin < 5
    tmp_filename=['idx_dump_' datestr(now)];
end;

if (nargin < 5)
    verbose=0;
end;


%check for file I/O
try
    stack=tom_emread(filename_stack); sz=stack.Header.Size; stack=stack.Value;
catch
    disp([filename_stack ' not readable !']);
    rethrow(lasterror);
end;

%apply mask... and (to do filter) to stack
mask=tom_create_mask(mask_st);
for i=1:size(stack,3)
    stack(:,:,i)=stack(:,:,i).*mask;
end;

if isdir(tmp_filename)
    delete([tmp_filename '/idx_*.em']);
else
    mkdir(tmp_filename);
end;



%calculate index stack
lookup=tom_av2_index_calc_recurs2(stack,tmp_filename,[1:size(stack,3)],[],tree_split,verbose,[]);


clear('stack');  %free da memory !


% ... to do in recursion 

% build index Stack out of tmp directory 
% left/right branch is 0 1 encoded 
% 0..n encoding is used for calculating a  position in the index stack 

dir_st=dir([tmp_filename '/*.em']);

for i=1:size(dir_st,1)
    [a b]=strtok(dir_st(i).name,'.');
    [a b]=strtok(a,'_');
    a=strrep(b,'_','');
    ind_deg=tom_av2_index_bintree_not_2_index(a,tree_split.num_of_classes);
end;

max_ind=max(ind_deg);


%allocate memory for index stack
idx_stack=zeros(sz(1),sz(2),max_ind+1);

for i=1:size(dir_st,1)
    tmp=tom_emread([tmp_filename '/' dir_st(i).name]);
    [a b]=strtok(dir_st(i).name,'.');
    [a b]=strtok(a,'_');
    a=strrep(b,'_','');
    ind_deg=tom_av2_index_bintree_not_2_index(a,tree_split.num_of_classes);
    idx_stack(:,:,ind_deg)=tmp.Value;
end;

tom_emwrite(filename_index_stack,idx_stack);
rest=strtok(filename_index_stack,'.');
save([rest '_lookup.mat'],'lookup');

%clean up folder
rmdir(tmp_filename,'s');

function lookup=tom_av2_index_calc_recurs2(stack_in,path_out,part,name,tree_split,verbose,lookup)


if strcmp(tree_split.Method,'2_halfs')==1 
    if (length(part) >= tree_split.num_of_classes)
        num_of_cl=tree_split.num_of_classes;
    else
        num_of_cl=length(part);
    end;
    N = floor(length(part) /num_of_cl);
    start=1;
    for i=1:num_of_cl
        if (i==tree_split.num_of_classes)
            stop=length(part);
        else
            stop=start+N-1;
        end
        part_list{i,:} = part(start:stop);
        start=start+N;   
    end;
end;

if strcmp(tree_split.Method,'classify')==1 || strcmp(tree_split.Method,'classify balanced')==1 ||  strcmp(tree_split.Method,'multiref')==1 
   
    %tree_split.Method='multiref';
    
    if (strcmp(tree_split.Method,'classify balanced') )
        equal_flag=1;
    else
        equal_flag=0;
    end;
    
    if (length(part) >= tree_split.num_of_classes)
        num_of_cl=tree_split.num_of_classes;
    else
        num_of_cl=length(part);
    end;
    
    if  strcmp(tree_split.Method,'classify')==1 || strcmp(tree_split.Method,'classify balanced')==1
        
        for iai=1:20
            try
                idx=tom_pca_and_k_means(stack_in(:,:,part),num_of_cl,tree_split.Values.Eigenvalues,tree_split.Values.Binning,equal_flag,0);
                break;
            catch
                idx=tom_pca_and_k_means(stack_in(:,:,part),num_of_cl,tree_split.Values.Eigenvalues,tree_split.Values.Binning,equal_flag,0);
            end;
            
        end;
        
        
        for i=1:tree_split.num_of_classes
            part_list{i,:}=part(idx==i);
        end;
        centros=zeros(size(stack_in,1),size(stack_in,2),num_of_cl);
        for i=1:num_of_cl
            tmp_list=part_list{i,:};
            for ii=1:length(tmp_list)
                centros(:,:,i)=centros(:,:,i)+stack_in(:,:,tmp_list(ii));
            end;
            centros(:,:,i)= centros(:,:,i)./length(tmp_list);
        end;
    end;
   
   if strcmp(tree_split.Method,'multiref')==1
       disp('multiref');
       [idx centros]=tom_av2_multi_ref_classify(stack_in(:,:,part),num_of_cl,4);
       for i=1:tree_split.num_of_classes
            part_list{i,:}=part(idx==i);
        end;
   end;
   
   
end;



for i=1:num_of_cl

    
    if (isempty(part_list{i,:}))
    else
        offset=2-i;
        if (offset<0)
            offset=0;
        end;
        name(length(name)+offset)=i-1;
        if (name(1)==2)
            disp('');
        end;
        tmp_list=part_list{i,:};
        disp(['name: ' num2str(name) '   members: '  num2str(tmp_list)  ]);
        mean_tmp=zeros(size(stack_in,1),size(stack_in,2));
        mean_tmp=tom_emheader(mean_tmp);
        mean_tmp.Value=centros(:,:,i);
        mean_tmp.Header.Comment=char(strrep(num2str(tmp_list),' ',''));
        t_name=strrep(num2str(name),' ','');
        tom_emwrite([path_out '/idx_' t_name '.em'],mean_tmp);

        
        if (length(tmp_list) > 1)
            lookup=tom_av2_index_calc_recurs2(stack_in,path_out,tmp_list,name,tree_split,verbose,lookup);
        end;
        lookup{tom_av2_index_bintree_not_2_index(strrep(num2str(name),' ',''),tree_split.num_of_classes),2}=tmp_list;
        lookup{tom_av2_index_bintree_not_2_index(strrep(num2str(name),' ',''),tree_split.num_of_classes),1}=strrep(num2str(name),' ','');
        lookup{tom_av2_index_bintree_not_2_index(strrep(num2str(name),' ',''),tree_split.num_of_classes),3}=num_of_cl;
    end;
   

    
 end;

