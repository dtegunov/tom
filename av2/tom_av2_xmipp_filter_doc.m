function tom_av2_xmipp_filter_doc(f_doc,f_outputdoc,operator,value,fieldname,abs_of_field,comb_flag,f_output_sel)
%TOM_AV2_XMIPP_FILTER_DOC filters xmipp doc file
%   
%
%  TOM_AV2_XMIPP_FILTER_DOC(f_doc,f_outputdoc,operator,value)
%
%  TOM_AV2_XMIPP_FILTER_DOC filters xmipp doc file according 2 given operators 
%  
%
%PARAMETERS
%
%  INPUT
%   f_doc               *.doc filename use abs filename 
%   f_outputdoc         name of output doc
%   operator            filter operator (< > == and )
%   value               filter value
%   fieldname           doc fieldname (ccf,xoff ...)
%   abs_of_field        (0) use abs of field (0/1)
%   comb_flag           comb flag for more the one condition (&&,||)
%   f_output_sel        (opt) sel outputname
%
%  OUTPUT
%
%EXAMPLE
%      
%   tom_av2_xmipp_filter_doc('Iter_17_current_angles.doc','f_doc.doc',{'<';'<'}
%   ,[15 15],{'xoff';'yoff'},[1 1],{'&&'},'f_sel.sel');
%   %filter shifts smaller than 15 in x and y
%  
%  
%   tom_av2_xmipp_filter_doc('Iter_27_current_angles.doc','f_doc.doc',{'>'},
%  [72],{'tilt'},[1],{'&&'},'f_sel.sel');
%  %filter out top views from 20s reconstruction ...corresponds 2 a low
%  %tilt angle of 0 and high tilt angle of 18 in xmipp_projmatch restrict
%  %tilt angles
%
% tom_av2_xmipp_filter_doc('sym_doc.doc','sym_doc.doc',{'=='},{'mod_7_low.sel'},{'name'},[0],{'&&'});
% %filter by names in mod_7_low.sel
%
% tom_av2_xmipp_filter_doc('Iter_22_current_angles.doc','f_doc.doc',{'ismember'},1:1000,{'ref'},[1],{'&&'},'f_sel.sel');
% %get only particles with a class reference which is in 1:1000
% 
%
%REFERENCES
%
%NOTE 
%  
% fieldnames are:
%     flip
%     ref
%     yoff
%     xoff
%     psi
%     tilt
%     rot
%     part_idx
%     name
%     run_num
%     maxCC
%
%
%
%SEE ALSO
%   tom_av2_em_classify3d
%
%   created by FB 08/09/09
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

try
    doc=tom_xmippdocread(f_doc);
catch ME
    disp(['Error reading: ' f_doc]);
    error(ME.message);
end;

%alloc memory 
in=zeros(length(doc),length(fieldname));

cell_flag=0;
for ii=1:length(fieldname)
    if (iscell(value))
        cell_flag=1;
        tmp_sel=importdata(value{ii});
    end;    
    for i=1:length(doc)
        act_val=getfield(doc(i),fieldname{ii});
        if (abs_of_field(ii))
            act_val=abs(act_val);
        end;
        if (cell_flag==0)
            if (strcmp(operator{ii},'ismember'))
                if(ismember(act_val,value))
                    in(i,ii)=1;
                end;   
            else
                cond=[num2str(act_val) ' ' operator{ii} ' ' num2str(value(ii))];
                if (eval(cond))
                    in(i,ii)=1;
                end;    
            end;
           
        else
            if ( isempty(ismember(act_val,tmp_sel.textdata))==0 );
                in(i,ii)=1;
            end;
        end;
        
    end;
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

disp(['org. doc has length of: ' num2str(length(doc))]);
disp(['filt. doc has length of: ' num2str(length(idx))]);


tmp_doc=doc(idx);

tmp_doc(1).header=doc(1).header;
tmp_doc(1).part_idx_unique=doc(1).part_idx_unique;

tom_xmippdocwrite(f_outputdoc,tmp_doc);

fid=fopen(f_output_sel,'wt');
for i=1:length(idx)
    fprintf(fid,'%s 1  \n',doc(idx(i)).name);
end;
fclose(fid);




















