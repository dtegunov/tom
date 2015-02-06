function tom_av2_xmipp_filter_doc_by_sel(input_doc,output_doc,sel_file,use_names,find_what,replace_with)
%tom_av2_xmipp_filter_doc_by_sel filters a doc by a given sel file to have
%only the parts from the sel file
%
%   tom_av2_xmipp_filter_doc_by_sel(input_doc,output_doc,sel_file,use_names,find_what,replace_with)
%
%  
%
%PARAMETERS
%
%  INPUT
%   input_doc         filename of the input sel
%   output_doc        string in the sel which should be replaced
%   sel_file          new string (according 2 find and replace in a text editor)
%   use_names         (1) use filenames for matching only
%   find_what         ('') find what  
%   replace_with      ('') replace with               
%
%EXAMPLE
%  
%
%  tom_av2_xmipp_filter_doc_by_sel('Iter_3_current_angles_1.doc','p_cl11.doc','parts_11.sel',1,'_high_','_low_');  
%
%REFERENCES
%
%SEE ALSO
%   
% tom_av2_xmipp_filter_doc
%  
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

if (nargin<4)
    use_names=1;
end;

if (nargin<5)
    find_what='';
end;

%generate template 4 fgrep -f
if (use_names==1)
    if (isempty(find_what))
        call=['cat ' sel_file  ' | awk ''{split($1,a,"/"); print a[length(a)]}'' > template.sel '];
    else
        call=['cat ' sel_file  ' | awk ''{split($1,a,"/"); print a[length(a)]}'' | awk ''{gsub("' find_what '","' replace_with '"); print$0; }''  > template.sel '];
    end;
else
    if (isempty(find_what))
        call=['cat ' sel_file  ' | awk ''{print $1}'' > template.sel '];
    else
        call=['cat ' sel_file  ' | awk ''{print $1}'' | awk ''{gsub("' find_what '","' replace_with '"); print$0; }''  > template.sel '];
    end;
end;

disp(call);
unix(call);

%Generate empty doc
call=['head -1 ' input_doc ' > ' output_doc];
disp(call);
unix(call);

%match doc files
call=['cat ' input_doc ' | fgrep -A1 -f template.sel | grep -v "^\--" >>  ' output_doc];
disp(call);
unix(call);

%remove template
disp('rm template.sel');
unix('rm template.sel');

