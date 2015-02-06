function [stack_comp align2d]=tom_av2_compress_stack(stack,num_of_classes,equal_flag,pre_alg_th,eigs,binning,cluster_method,outputfold,verbose)
%tom_av2_compress_stack compresses a stack by clustering  
%
%   [stack_comp align2d]=tom_av2_compress_stack(stack,num_of_classes,equal_flag,pre_alg_th,eigs,binning,cluster_method,outputfold,verbose)
%PARAMETERS
%
%  INPUT
%   stack                particle stack or tom-picklist
%   num_of_classes       number of classes    
%   equal_flag           ('not_equal')
%   pre_alg_th           (0.1) cc thresh for pre alignment use -2 to switch
%                               of prealg
%
%   eigs                 (1,20)  vector which determines the used eigenvectors   
%   binning              (0) the align structure
%   cluster_method       (k-means) 
%   outputfold           (./comp_stack) use '' for no output writtein 2 hd
%   verbose              (0) output on/off 1 2 Arpack output levels
%
%  OUTPUT
%   stack_comp           compressed stack
%   align2d              align2d struct for stackbrowser  
%   
%
%
%EXAMPLE
%   
%   st_comp=tom_av2_compress_stack('pickLists/pick_high_1.em.mat',50,'not_equal',0.1,1:10,1,'k-means','out',0);
%
%REFERENCES
%
%SEE ALSO
%   tom_av2_stackbrowser,tom_av2_uncompress_stack
%
%   created by fb (eckster)
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


if (ischar(stack))
    pl_flag=0;
    try
        load(stack);
        stack=tom_av2_xmipp_picklist2stack(align2d,'','','','gradient&mean0+1std',1,[64 64]);
        pl_flag=1;
    catch
        pl_flag=0;
    end;
    
    if (pl_flag==0)
        
        if (exist(stack,'file') )
            stack=tom_emread(stack);
            stack=stack.Value;
        else
            stack=tom_av2_xmipp_picklist2stack(stack,'','','','gradient&mean0+1std',1,[64 64]);
        end;
    end;
    
end;

if (isempty(outputfold)==0)
    warning off; mkdir(outputfold); warning on;
    tom_emwritec([outputfold '/stack_org.em'],stack);
end;




normflag=1;
sz=size(stack);


if (pre_alg_th ~= -2)
     disp('Aligning Stack ...')
     m_tmp=tom_spheremask(ones(sz(1),sz(2)),round(sz(1)./7),5);
     [stack all_tmpl cc_out cc_sum]=tom_os3_alignStack2(stack,sum(stack,3),'default',[1 0],'mean0+1std','',[6 3],m_tmp);
     v_sel=cc_out(:,end)>pre_alg_th;   
     disp('done!');
end;

if (isempty(outputfold)==0)
    warning off; mkdir(outputfold); warning on;
    tom_emwritec([outputfold '/stack_alg.em'],stack);
end;


if ischar(stack)
    dimension=3;
else
    dimension=2;
end;

if (dimension==2)
    disp('Reshaping ...')
    reshapedStack=zeros(sz(3),(sz(1)./2^binning).*(sz(2)./2^binning) );
    for i=1:sz(3)
        im = stack(:,:,i);
        im = tom_bin(im,binning);
        im = reshape(im,1,[]);
        if (normflag == 1)
            im = tom_norm((im+100).*2,'phase');
        end
        reshapedStack(i,:)=im;
    end;
    disp('done!');
else
    %3d
    disp('not implemented!');
end;


%calc pca
if (max(eigs) > 0)
    disp('Calc Eignevectors ...');
    eigs_in=max(eigs);
    [scores,coefs,eigenValues]=tom_calc_pca(double(reshapedStack),eigs_in,'pca','','','','','','',verbose);
    disp('done!')
else
    scores=reshapedStack;
end;

if (strcmp(cluster_method,'k-means'))
    disp('Clustering');
    warning off;
    if ((size(scores,2)./20) >  num_of_classes)
        [idx centriod sum_dist dists] = kmeans(scores',num_of_classes,'Start','Cluster', 'Display','Iter','MaxIter',40,'EmptyAction','drop');
    else
        [idx centriod sum_dist dists] = kmeans(scores',num_of_classes,'Display','Iter','MaxIter',40,'EmptyAction','drop');
    end;
    warning on;
    disp('done!');
end;

if (pre_alg_th ~= -2)
    h_cl=max(idx)+1;
    idx=(idx.*(v_sel==1)) + ((v_sel==0) .*  h_cl);
end;


idx_unique=unique(idx);

stack_comp=zeros(sz(1),sz(2),length(idx_unique));

align2d=tom_av2_create_alignfromstack([sz(1) sz(2) length(idx_unique)]);

for i=1:length(idx_unique)
    idx_cl=find(idx==idx_unique(i));
    stack_comp(:,:,i)=tom_norm(sum(stack(:,:,idx_cl),3)./length(idx_cl),'mean0+1std');
    align2d(1,i).dataset=idx_cl;
    align2d(1,i).ref_class=length(idx_cl);
    if (exist('cc_out','var'))
        align2d(1,i).ccc=mean(cc_out(idx_cl,end));
    end;
end;


if (isempty(outputfold)==0)
    warning off; mkdir(outputfold); warning on;
    save([outputfold '/comp_stack.mat'],'align2d');
    tom_emwrite([outputfold '/comp_stack.em'],stack_comp);
end;










% if ((equal_flag==1) && (num_of_classes==2))
%     
%     
%     if length(find(idx==1)) > length(find(idx==2))
%         diff=round((length(find(idx==1))-length(find(idx==2)))./2);
%         for i=1:diff
%             [a idx_min]=max(dists(:,1));
%             dists(idx_min,1)=min(dists(:,1))-1;
%             idx(idx_min)=2;
%         end;
%         return;
%     end;
%     
%     if length(find(idx==1)) < length(find(idx==2))
%         diff=round((length(find(idx==2))-length(find(idx==1)))./2);
%         for i=1:diff
%             [a idx_min]=min(dists(:,2));
%             dists(idx_min,2)=min(dists(:,2))-1;
%             idx(idx_min)=1;
%         end;
%         return;
%     end;
%     
% end;
%     
    %     for i=1:num_of_classes
    %
    %     end;
    
    
    
    
    
    
    
