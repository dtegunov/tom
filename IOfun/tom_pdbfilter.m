function [pdb_out idx]=tom_pdbfilter(pdb,operator,value,fieldname,abs_of_field,comb_flag,f_output_pdb)
% tom_pdbfilter filters a pdb 
% 
% Input
%  pdb                 struct or filename
%  operator            filter operator (< > == and )
%  value               filter value
%  fieldname           pdb fieldname (segID,resSeq ...)
%  abs_of_fieldname    (0) use abs of field (0/1)
%  comb_flag           (&&) comb flag for more the one condition (&&,||)
%  f_output_pdb        (opt.) name of output pdb
%  
% Output:
%  pdb_out             filtered pdb 
%  idx                 index of used atoms   
% 
% Example:
% 
% matlabpool open local 8;
%
% %filter 4 one subunit
% pdb_out=tom_pdbfilter(pdb,{'strcmp'},{'B6_2'},{'segID'});
% tom_vol2chimera(pdb,pdb_out);
%
% %filter 1 one Alpha ring 
% pdb_out=tom_pdbfilter(pdb,{'regexp'},{'^A.*_1'},{'segID'});
% tom_vol2chimera(pdb,pdb_out);
%
% %filter 1 one Alpha ring and by coordinates 
% pdb_out=tom_pdbfilter(pdb,{'regexp','<'},{'^A.*_1','245'},{'segID','X'}); 
% tom_vol2chimera(pdb,pdb_out);
%
% %get only C-Alpha from Beta-ring
% pdb_out=tom_pdbfilter(pdb,{'regexp','strcmp'},{'^B.*_1','CA'},{'segID','AtomName'}); 
% tom_vol2chimera(pdb,pdb_out);
%
% %get only C-Alpha from Alpha ring out of org 1ryp pdb 
% pref_pdb_tmpl=tom_pdbfilter(org_pdb_tmpl,{'ismember','strcmp'},{{'A','B','C','D','E','F','G'},'CA'},{'chainID','AtomName'}); 
%
% References:
% 
%
%SEE ALSO
%   tom_vol2chimera
%
%   created by fb
%
%   Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%   Journal of Structural Biology, 149 (2005), 227-234.
%
%   Copyright (c) 2004-2010
%   TOM toolbox for Electron Tomography
%   Max-Planck-Institute of Biochemistry
%   Dept. Molecular Structural Biology
%   82152 Martinsried, Germany
%   http://www.biochem.mpg.de/tom

if (nargin < 5)
    abs_of_field=zeros(length(fieldname),1);
end;

if (nargin < 6)
    comb_flag{1}='&&';
end;

if (nargin < 7)
    f_output_pdb='';
end;

if (ischar(pdb))
    pdb=tom_pbdread(pdb);
end;

%alloc memory 
atom=pdb.Model.Atom;

in=zeros(length(atom),length(fieldname));


for ii=1:length(fieldname)
    
    parfor i=1:length(atom)
        in(i,ii)=filter_st(atom(i),fieldname{ii},operator{ii},value{ii},abs_of_field(ii));
        if (mod(i,10000)==0)
            disp(num2str(i));
        end;
    end;
    disp(' ');
end;

call='';
for i=1:length(fieldname)
    
    call=[call 'in(:,' num2str(i) ') ' ];
    if (i==length(fieldname))
        break;
    end;
    if (strcmp(comb_flag{i},'&&'))
        call=[call '.*'];
    else
        call=[call '+'];
    end;
end;

idx=find(eval(call)>0);

pdb_out=pdb;
pdb_out.Model=rmfield(pdb.Model,'Atom');
pdb_out.Model.Atom=pdb.Model.Atom(idx);


disp(['org. pdb has length of: ' num2str(length(pdb.Model.Atom))]);
disp(['filt. pdb has length of: ' num2str(length(pdb_out.Model.Atom))]);



if (isempty(f_output_pdb)==0)
    tom_pdbwrite(f_output_pdb,pdb_out);
end;



function in=filter_st(atom,fieldname,operator,value,abs_of_field)

in=0;
act_val=getfield(atom,fieldname);
if (abs_of_field)
    act_val=abs(act_val);
end;

if (strcmp(operator,'ismember'))
    if(ismember(act_val,value))
        in=1;
    end;
end;
if (strcmp(operator,'strcmp'))
    if(strcmp(act_val,value))
        in=1;
    end;
end;
if (strcmp(operator,'regexp'))
    if(regexp(act_val,value))
        in=1;
    end;
end;

if (ismember(operator,{'>','<','=='}))
    tmp_val=value;
    if isnumeric(value)
        tmp_val=num2str(value);
    end;
    cond=[num2str(act_val) ' ' operator ' ' tmp_val];
    
    if (eval(cond))
        in=1;
    end;
end;










%old code 
% act_val=getfield(atom(i),fieldname{ii});
%         if (abs_of_field(ii))
%             act_val=abs(act_val);
%         end;
%         
%         if (strcmp(operator{ii},'ismember'))
%             if(ismember(act_val,value))
%                 in(i,ii)=1;
%             end;
%         end;
%         if (strcmp(operator{ii},'strcmp'))
%             if(strcmp(act_val,value{ii}))
%                 in(i,ii)=1;
%             end;
%         end;
%         if (strcmp(operator{ii},'regexp'))
%             if(regexp(act_val,value{ii}))
%                 in(i,ii)=1;
%             end;
%         end;
%         
%         if (ismember(operator{ii},{'>','<','=='}))
%             tmp_val=value{ii};
%             if isnumeric(value{ii})
%                 tmp_val=num2str(value{ii});
%             end;
%             cond=[num2str(act_val) ' ' operator{ii} ' ' tmp_val];
%             
%             if (eval(cond))
%                 in(i,ii)=1;
%             end;
%         end;
