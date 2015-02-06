function [part_stack_alg tmp_iddx]=tom_av2_xmipp_doc2em(doc_filename,class_in,alg_flag,output_sel,find_what,replace_with,bin)
%tom_av2_xmipp_doc2em extracts a class from the xmipp .doc file and creates
% an EM stack
%
%  part_stack_alg=tom_av2_xmipp_doc2em(doc_filename,class_in,alg_flag)
%
%   generates aligned 2d stack from doc file (ML2d,xmipp)
%
%PARAMETERS
%
%  INPUT
%   doc_filename  name of the xmipp *.doc file from ml2d ml3d,ProJ match
%   class_in      ('') class nr use '' for all classes
%   alg_flag      (1) flag for alignment 0 for no alignment  ... or 2 for
%                 centering only
%   output_sel    outpt sel  
%   find_what     to be replaced
%   replace_with  replacement string
%   bin           binning for part_stack
%
%  OUTPUT
%   part_stack_alg            aligned stack
%   tmp_iddx                  index of the particel in respect of the input doc-file     
%
%EXAMPLE
%    matlabpool open local 8;
%    extract class nr 8 from the xmipp doc-file:
%    stack=tom_av2_xmipp_doc2em('ml2d_it000008.doc',3);
%
%REFERENCES
%
%   xmipp
%
%SEE ALSO
%   tom_xmippsellread, tom_xmippdocread
%
%   created by FB  150909
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

if (nargin < 2)
    class_in='';
end;
    
if (nargin < 3)
    alg_flag=1;
end;


%     sym='C1';
%     numc=1;
%     sym_angs(1)=0;


if (nargin<4)
    output_sel='';
end;
    
if (nargin<7)
    bin=0;
end;


if (isempty(output_sel)==0)
    fp=fopen(output_sel,'wt');
end;

% if (strcmp(sym,'C1')==0)
%     numc=str2double(sym(2));
%     sym_angs(1)=0;
%     for i=2:numc;
%         sym_angs(i)=360./numc;
%     end;
% end;

if (isstruct(doc_filename))
    st=doc_filename;
else
    fprintf('\n%s ', ['Reading  ' doc_filename ':']);
    st=tom_xmippdocread(doc_filename);
    fprintf('%s \n', ['...done!']);
    disp('');
end;

fprintf('%s ', ['Converting Angles: ' ]); 
tmp=zeros(size(st,1),1);
proj_angles=zeros(size(st,1),3);
tmp_angles=zeros(size(st,1),3);

for i=1:size(st,1)
    tmp(i)=st(i).ref;
    [xx,angles] = tom_eulerconvert_xmipp(st(i).rot, st(i).tilt, st(i).psi);
    proj_angles(i,:) = angles';
    tmp_angles(i,:)=[st(i).rot st(i).tilt st(i).psi];
    if (mod(i,2000)==0)    
        fprintf('%s','.');
    end;
end;
fprintf('%s \n','done!');



fprintf('%s ', ['Allocating Memory: ' ]); 
num_of_ref=max(tmp);
try
    im_tmp=tom_spiderread([st(1).name]);
    path_flag='abs_man';
catch
    im_tmp=tom_spiderread(['./' st(1).name]);
    path_flag='rel_man';
end


for i=1:num_of_ref
    indx(i).filename{1}='';
end;

fprintf('%s \n', ['...done! ' ]); 


if (exist('f_alignment','var') && isempty(f_alignment)==0 )
    load(f_alignment);
end;


fprintf('%s ', ['Building Classes: ' ]); 

vect=[st(1:end).ref];

if (isempty(class_in))
    tmp_iddx=1:length(st);
else
    %tmp_iddx=find(vect==class_in);
    tmp_iddx=find(ismember(vect,class_in));
end;


if (nargout>0)
    part_stack_alg=zeros(size(im_tmp.Value,1)./(2^bin),size(im_tmp.Value,2)./(2^bin),length(tmp_iddx));
end;

% try
%     parfor i=tmp_iddx
%         
%       
%         if (strcmp(path_flag,'rel_man'))
%             tmp_name=['./' st(i).name];
%         else
%             tmp_name=[st(i).name];
%         end;
%         im=tom_spiderread(tmp_name);
%         
%         
%         
%         if (st(i).flip==0)
%             if (alg_flag>0)
%                 if (alg_flag==1)
%                     im_tmp_alg=tom_rotate(tom_shift(im.Value,[st(i).xoff st(i).yoff]),tmp_angles(i,3));
%                 end;
%                 if (alg_flag==2)
%                     im_tmp_alg=tom_shift(im.Value,[st(i).xoff st(i).yoff]);
%                 end;
%                 
%             else
%                 im_tmp_alg=im.Value;
%             end;
%             
%         else
%             if (alg_flag>0)
%                 im.Value=tom_mirror(im.Value,'x');
%                 im.Value=im.Value;
%                  if (alg_flag==1)
%                     im_tmp_alg=tom_rotate(tom_shift(im.Value,[-st(i).xoff st(i).yoff]),tmp_angles(i,3));
%                  end;
%                  if (alg_flag==2)
%                      im_tmp_alg=tom_shift(im.Value,[-st(i).xoff st(i).yoff]);
%                  end
%             else
%                 im_tmp_alg=im.Value;
%             end;
%             
%         end;
%         
%         if (nargout>0)
%             part_stack_alg(:,:,i)=im_tmp_alg;
%         end;
%         
%         if (isempty(output_sel)==0)
%             new_name=strrep(tmp_name,find_what,replace_with);
%             folder_n=fileparts(new_name);
%             if (exist(folder_n,'dir')==0 )
%                 mkdir(folder_n);
%             end;
%             if (strcmp(new_name,tmp_name))
%                 error('name is the same after find&replace ...check pattern');
%             end;
%             
%             all_new_names{i}=new_name;
%             tom_spiderwrite(new_name,im_tmp_alg);
%         end;
%         
%         if (mod(i,50)==0)
%             fprintf('%s','.');
%         end;
%     end;
% catch Me
    %disp(Me.message);
    for i=1:length(tmp_iddx)
        
        if (strcmp(path_flag,'rel_man'))
            tmp_name=['./' st(tmp_iddx(i)).name];
        else
            tmp_name=[st(tmp_iddx(i)).name];
        end;
        im=tom_spiderread(tmp_name);
        all_new_names{i}=tmp_name;  
        if (st(tmp_iddx(i)).flip==0)
            if (alg_flag==1)
                im_tmp_alg=tom_rotate(tom_shift(im.Value,[st(tmp_iddx(i)).xoff st(tmp_iddx(i)).yoff]),tmp_angles(tmp_iddx(i),3));
            else
                im_tmp_alg=im.Value;
            end;
            part_stack_alg(:,:,i)=tom_bin(im_tmp_alg,bin);
        else
            if (alg_flag==1)
                im.Value=tom_mirror(im.Value,'x');
                im.Value=im.Value;
                im_tmp_alg=tom_rotate(tom_shift(im.Value,[-st(tmp_iddx(i)).xoff st(tmp_iddx(i)).yoff]),tmp_angles(tmp_iddx(i),3));
            else
                im_tmp_alg=im.Value;
            end;
            part_stack_alg(:,:,i)=tom_bin(im_tmp_alg,bin);
        end;
      
        if (mod(i,50)==0)
            fprintf('%s','.');
        end;
    end;
%end;


fprintf('%s \n', ['...done! ' ]); 






%gen match
if (isempty(output_sel)==0)
    for i=1:length(all_new_names)
        fprintf(fp,'%s 1\n',all_new_names{i});
    end;
    fclose(fp);
end;








