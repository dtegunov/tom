function tom_av2_xmipp_mergesel2(in_sels,out_sel,extr_flag,operator,use_path)
%  tom_av2_xmipp_mergesel2 merges n .sel filese
%
%    
%       tom_av2_xmipp_mergesel2(in_sel,out_sel,extr_flag,operator)
%    
%    PARAMETERS
%    
%      INPUT
%       in_sels        input sel of sel-files
%       out_sel        filenme of output sel
%       extr_flag      extraction flag can be 
%                      'full' full path is used 2 merge
%                      'part' particle name is used 2 merge (if root path differs) 
%                      'num'  only particle number is used 2 merge         
% 
%       operator       operator for mergeing
%                      '&&'  use only particles which are in every set
%                      '||'  use all particels of all sets (== cat and unique )
%       use_path       use path of in_sel nr works only with && operator       
%    
%                       
%      OUTPUT
%      
%    
%    EXAMPLE
%        
%      tom_av2_xmipp_mergesel2({'my_sel1.sel';'my_sel2.sel'},'merge.sel','full','&&');
%      
%      %merge with particle name (path ignored) use path of my_sel2.sel for  output
%      tom_av2_xmipp_mergesel2({'my_sel1.sel';'my_sel2.sel'},'merge.sel','part','&&',1); 
%
%    REFERENCES
%    
%    NOTE:
%   
%    
%    
%  
%    SEE ALSO
%      
%    
%       created by FB 15/02/10
%    
%       Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%       Journal of Structural Biology, 149 (2005), 227-234.
%    
%       Copyright (c) 2004-2007
%       TOM toolbox for Electron Tomography
%       Max-Planck-Institute of Biochemistry
%       Dept. Molecular Structural Biology
%       82152 Martinsried, Germany
%       http://www.biochem.mpg.de/tom

if (nargin < 5)
    use_path=0;
end;


if (use_path > 0 && strcmp(operator,'||'))
    error('operator || and use_path > 0 not defined');
end;

fid=fopen(out_sel,'wt');

for i=1:length(in_sels)
    tmp_sel=importdata(in_sels{i});
    for ii=1:length(tmp_sel.textdata)
        if (strcmp(extr_flag,'full'))
            proc{ii}=tmp_sel.textdata{ii};
        end;
        if (strcmp(extr_flag,'part'))
            [a b c]=fileparts(tmp_sel.textdata{ii});
            proc{ii}=[b c];
         end;
        if (strcmp(extr_flag,'num'))
            tmp=tom_get_max_numeric({tmp_sel.textdata{ii}});
            proc{ii}=num2str(tmp);
        end;
    end;
    
    if (i>1)
        if (strcmp(operator,'&&'))
            [res idxa idxb]=intersect(res,proc);
            if (i==use_path)
                idx4path=idxb;             
            end;
            if (i-1==use_path)
                idx4path=idxa;             
            end;
        end;
        if (strcmp(operator,'||'))
            res=cat(2,res,proc);
            res=unique(res);
        end;
    else
        res=proc;
    end;
    clear('proc');
    if (i==use_path)
        sel4path=tmp_sel.textdata;
    end;
end;


if (use_path > 0 )
    res=sel4path(idx4path);
end;
    
for i=1:length(res)
  fprintf(fid,'%s 1 \n',res{i});  
end;

fclose(fid);










