function align2d=tom_av2_pickList_extend(pickList,num_of_frames,num_lead,find_what,replace_with,add,f_pickListOut)
%  tom_av2_pickList_extend extends a pickList for sumed images 2 frames
%  
%     mod_list=tom_av2_pickList_extend(pickList,num_of_frames,find_what,replace_with,add,f_pickListOut)
%  
%  PARAMETERS
%  
%    INPUT
%     pickList                 input picklist
%     num_of_frames            number of frames 4 extension   
%     num_lead                 number of leading zeros
%     find_what                part 2 be replaced 4 the frame name
%     replace_with             replacement 
%     add                      added part usually extension (.em,.mrc ...)
%     f_pickListOut            (opt.) filename 4 output of the picklist  
%     
%    OUTPUT
%     mod_pickList             extended pickList
% 
%  EXAMPLE
%  
%  mod_pickList=tom_av2_pickList_extend(align2d,40,2,'.em','-em/K_','.em')
%
%  REFERENCES
%  
%  SEE ALSO
%     ...
%  
%     created by FB 25/04/13
%  
%     Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%     Journal of Structural Biology, 149 (2005), 227-234.
%  
%     Copyright (c) 2004-2007
%     TOM toolbox for Electron Tomography
%     Max-Planck-Institute of Biochemistry
%     Dept. Molecular Structural Biology
%     82152 Martinsried, Germany
%     http://www.biochem.mpg.de/tom

if (nargin < 7)
    f_pickListOut='';
end;

if (isstruct(pickList))
    align2d=pickList;
else
    load(pickList);
end;

f_string=['%0' num2str(num_lead+1) 'd'];

org_filenames={align2d(1,:).filename};
align2d=rmfield(align2d,'filename');
for i=1:size(align2d,2)
    for ii=1:num_of_frames
        [a b c]=fileparts(org_filenames{i});
        new_names{ii}=[a filesep strrep([b c],find_what,replace_with) sprintf(f_string,ii) add];
     end;
     align2d(1,i).Orgfilename=org_filenames{i};
     align2d(1,i).filename=new_names;
    
    if (mod(i,2000)==0)
        disp([num2str(i) ' processed!']);
    end;
end;

if (isempty(f_pickListOut)==0)
    save(f_pickListOut,'align2d');
end;


