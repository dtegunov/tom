function [new_doc,idx1,idx2]=tom_av2_xmipp_mergedoc(doc1,doc2,merge_by,use_from2,fmerged_doc,idx1,idx2)
%TOM_AV2_MERGEDOC merges 2 doc xmipp .doc files 
%    m_doc=tom_av2_xmipp_mergedoc(doc1,doc2,merge_by,use_from2,fmerged_doc)
%
%PARAMETERS
%
%  INPUT
%   doc1                 first doc file (filename or struct)
%   doc2                 second doc file (filename or struct ) 
%   merge_by             flag for mergeing 
%                        "name" ...slow !!!
%                        "simple_index" 
%                        "double_index" 
%   use_from2            matlb cell  use{1}='name'
%                                    use{2}='xoff'   
%                                    use{3}='yoff'
%   fmerged_doc          filename for merged doc 
%   idx1                 index build for doc1 use precalc index ...speed!       
%   idx2                 index build for doc2 use precalc index ...speed!   
%
%
%  OUTPUT
%   new_doc            merged doc in memory  
%   idx1               index build for doc1
%   idx2               index build for doc2
%
%EXAMPLE
%   
% tom_av2_xmipp_mergedoc('Iter_12_current_angles_1.doc','sym_doc280_pcent.doc','double_index',{'name','xoff','yoff'});
%
%REFERENCES
%
%SEE ALSO
%   tom_xmippdocread
%
%   created by FB 10/08/2012
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

if (nargin < 5)
    fmerged_doc='';
end;

if (nargin < 6)
    idx1=[];
end;

if (nargin < 7)
    idx2=[];
end;




if (isstruct(doc1)==0)
    disp(['Reading ' doc1]);
    doc1=tom_xmippdocread(doc1);
    disp('done!');
end;

if (isstruct(doc2)==0)
    disp(['Reading ' doc2]);
    doc2=tom_xmippdocread(doc2);
    disp('done!');
end;

if (isempty(idx1))
    disp('building idx for doc1');
    idx1=build_idx(doc1,merge_by);
    disp('done!');
else
    if (length(idx1)~=length(doc1) )
        error(['length idx1 differs form length doc1']);
    end;
    if (length(idx1)~=length(unique(idx1)))
        error('idx1 is not unique');
    end;
    disp('using pre calc idx1');
end;

if (isempty(idx2))
    disp('building idx for doc2');
    idx2=build_idx(doc2,merge_by);
    disp('done!');
else
    if (length(idx2)~=length(doc2) )
        error(['length idx2 differs form length doc2']);
    end;
    if (length(idx2)~=length(unique(idx2)))
        error('idx2 is not unique');
    end;
    disp('using pre calc idx2 for doc2');
end;

[in_both idx_both]=ismember(idx1,idx2);
in_both=find(in_both);
idx_both=idx_both(in_both);
new_doc=doc1(in_both);

disp(['found ' num2str(length(in_both)) ' of ' num2str(length(doc1)) ' particles form doc1 in doc2' ]);

tic;
for i=1:length(idx_both)
    for ii=1:length(use_from2) %check for idx!!
        eval(['new_doc(i).' use_from2{ii} '=doc2(idx_both(i)).' use_from2{ii} ';']);
        new_doc(i).run_num=i;
    end;
    if (mod(i,100000)==0)
        disp([num2str(i) ' entries merged!']);
        toc;
        tic;
    end;
end;
toc;
disp('done!');
new_doc(1).header=doc1(1).header;
new_doc(1).part_idx_unique=0;
disp(' ');
if (isempty(fmerged_doc)==0)
    disp(['Writing: ' fmerged_doc]);
    tom_xmippdocwrite(fmerged_doc,new_doc);
    disp('done!');
end;




function idx=build_idx(doc,merge_by)

if (strcmp(merge_by,'name')==0)
    idx=zeros(length(doc),1);
end;

for n = 1:length(doc)
    
    [a b c]=fileparts(doc(n).name);
    if (strcmp(merge_by,'name'))
        idx{n}=[b c];
    end;
    
    if (strcmp(merge_by,'double_index'))
        tmp=b;
        [ctok rest]=strtok(tmp,'_');
        num_ind=ismember(ctok,['0' '1' '2' '3' '4' '5' '6' '7' '8' '9' ]);
        tmp_ind1=str2double(ctok(num_ind==1));
        [ctok rest]=strtok(rest,'_');
        tmp_ind2=str2double(strrep(rest,'_',''));
        idx(n)=(tmp_ind1.*1000000000) + tmp_ind2;
        if (isnan(idx(n)))
            error(['idx build error on part: ' b]);
        end;
    end;
    if (strcmp(merge_by,'single_index'))
        idx(n)=tom_get_max_numeric({b});
    end;
    
    if (mod(n,500000)==0)
        disp(num2str(n));
    end
end;

if (length(idx)~=length(unique(idx)))
    error('index build error');
end;
    

