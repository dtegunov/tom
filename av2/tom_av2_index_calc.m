function tom_av2_index_calc(filename_stack,filename_index_stack,mask_st,tree_split,tmp_filename,verbose)
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
lookup=tom_av2_index_calc_recurs(stack,tmp_filename,[1:size(stack,3)],[],tree_split,verbose,[]);


clear('stack');  %free da memory !


% ... to do in recursion 

% build index Stack out of tmp directory 
% left/right branch is 0 1 encoded 
% 0 1 encoding is used for calculating a  position in the index stack 

dir_st=dir([tmp_filename '/*.em']);

for i=1:size(dir_st,1)
    [a b]=strtok(dir_st(i).name,'.');
    [a b]=strtok(a,'_');
    a=strrep(b,'_','');
    ind_deg=tom_av2_index_bintree_not_2_index(a);
end;

max_ind=max(ind_deg);


%allocate memory for index stack
idx_stack=zeros(sz(1),sz(2),max_ind+1);

for i=1:size(dir_st,1)
    tmp=tom_emread([tmp_filename '/' dir_st(i).name]);
    [a b]=strtok(dir_st(i).name,'.');
    [a b]=strtok(a,'_');
    a=strrep(b,'_','');
    ind_deg=tom_av2_index_bintree_not_2_index(a);
    idx_stack(:,:,ind_deg)=tmp.Value;
end;

tom_emwrite(filename_index_stack,idx_stack);
rest=strtok(filename_index_stack,'.');
save([rest '_lookup.mat'],'lookup');

%clean up folder
rmdir(tmp_filename,'s');

function lookup=tom_av2_index_calc_recurs(stack_in,path_out,part,name,tree_split,verbose,lookup)


if strcmp(tree_split.Method,'2_halfs')==1 
    N = ceil(length(part) / 2);
    part_left = part(1:N);
    part_right = part(1+N:end);
end;

if strcmp(tree_split.Method,'classify')==1 
    idx=tom_pca_and_k_means(stack_in(:,:,part),2,tree_split.Values.Eigenvalues,tree_split.Values.Binning,0,0);
    part_left=part(idx==1);
    part_right=part(idx==2);
end;


if (isempty(part_left))
else
    name(length(name)+1)=0;
    disp(['members: ' num2str(part_left) ' name: ' num2str(name) ]);
    mean_tmp=zeros(size(stack_in,1),size(stack_in,2));
    mean_tmp=tom_emheader(mean_tmp);
    for i=1:length(part_left)
        mean_tmp.Value=mean_tmp.Value+stack_in(:,:,part_left(i));
    end;
    mean_tmp.Value=mean_tmp.Value./length(part_left);
    
    mean_tmp.Header.Comment=char(strrep(num2str(part_left),' ',''));
    t_name=strrep(num2str(name),' ','');
    tom_emwrite([path_out '/idx_' t_name '.em'],mean_tmp);    
    
    if (length(part_left) > 1)
        lookup=tom_av2_index_calc_recurs(stack_in,path_out,part_left,name,tree_split,verbose,lookup);
    end;
    lookup{tom_av2_index_bintree_not_2_index(strrep(num2str(name),' ','')),2}=part_left;
    lookup{tom_av2_index_bintree_not_2_index(strrep(num2str(name),' ','')),1}=strrep(num2str(name),' ','');
end;


if (isempty(part_right))
else
    name(length(name))=1;
    disp(['members ' num2str(part_right) ' name: ' num2str(name) ]);
    if (strcmp(strrep(num2str(name),' ','') ,'01'))
        disp('');
    end;
    mean_tmp=zeros(size(stack_in,1),size(stack_in,2));  
    mean_tmp=tom_emheader(mean_tmp);
    for i=1:length(part_right)
        mean_tmp.Value=mean_tmp.Value+stack_in(:,:,part_right(i));
    end;
    mean_tmp.Value=mean_tmp.Value./length(part_right);
    mean_tmp.Header.Comment=char(strrep(num2str(part_right),' ',''));
    t_name=strrep(num2str(name),' ','');
    tom_emwrite([path_out '/idx_' t_name '.em'],mean_tmp);   
    
    if (length(part_right) > 1)
        lookup=tom_av2_index_calc_recurs(stack_in,path_out,part_right,name,tree_split,verbose,lookup);
    end;
    
    lookup{tom_av2_index_bintree_not_2_index(strrep(num2str(name),' ','')),2}=part_right;
    lookup{tom_av2_index_bintree_not_2_index(strrep(num2str(name),' ','')),1}=strrep(num2str(name),' ','');
end;



