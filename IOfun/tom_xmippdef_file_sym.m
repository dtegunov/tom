function tom_xmippdef_file_sym(input_def,input_symed_doc,sym,f_def_sym)
%tom_xmippdef_file_sym syms a defocus part file according 2 given sym doc
%   
%
%  tom_xmippdef_file_sym(input_def,input_symed_doc,output_def,basepath,sym)
%
%  TOM_XMIPPDEF_FILE_SYM syms a defocus file
%                  
%                  
%PARAMETERS
%
%  INPUT
%   input_def               input defocus file
%   input_symed_doc         name of output em stack
%   sym_fact                factor for symmetry for e.g c2 is 2
%                                               for e.g d2 is 4
%   f_def_sym               filename for symmed def file
%
%  OUTPUT
%
%EXAMPLE
%    
%!cat p280_parts_defocus.txt | fgrep -f tmp_cl1.tmp > p280_parts_defocus_cl1_cut.txt
% tom_xmippdef_file_sym('p280_parts_defocus_cl1_cut.txt','p280_low_cl_1_cut.doc',2)
%
% 
%
%REFERENCES
%
%SEE ALSO
%   tom_xmippdocsym,tom_av2_xmipp_ctf_gen_defocus_partlist
%
%   created by FB 10/26/12
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


disp(['Reading ' input_symed_doc]);
doc_in=tom_xmippdocread(input_symed_doc);
disp('done!');


disp(['Reading ' input_def]);
tmp_fid=fopen(input_def,'r');
C = textscan(tmp_fid, '%s %f\n');
parts_and_def.textdata= C{1,1};
parts_and_def.data = C{1,2};
disp('done!');

l_org_doc=round(length(doc_in)./sym);

[is_in idx_def]=ismember({doc_in(:).name},parts_and_def.textdata);

for i=1:l_org_doc
    def{i}=[parts_and_def.textdata{idx_def(i)} ' ' num2str(parts_and_def.data(idx_def(i)))];
end;

zz=l_org_doc;
for i=1:sym-1
    idx=(i.*l_org_doc);
    for ii=1:l_org_doc
        zz=zz+1;
        def{zz}=[doc_in(idx+ii).name ' ' num2str(parts_and_def.data(idx_def(ii)) )];
    end;
end;

fid=fopen(f_def_sym,'wt');
for i=1:length(def)
    fprintf(fid,'%s\n',def{i});
end;
fclose(fid);



