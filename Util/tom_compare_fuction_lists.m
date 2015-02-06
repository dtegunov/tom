function tom_compare_fuction_lists(old_folder,new_folder,wildcard,output)
%  TOM_COMPARE_FUCTIONS_LISTS compares functions of 2 folders
%                             useful for release management
%  
%     
%     tom_compare_fuctions_lists(folder1,folder2,prefix,ouputfolder)
%  
%  
%  PARAMETERS
%  
%    INPUT
%     old_folder     folder containing files
%     new_folder     folder containing files
%     wildcard       wildcard for searching
%     output         (out_lists) outputfolder for cp scripts and lists     
%                         
%  
%    
%    OUTPUT
%  
%    NOTE:
%    outputfolder/missing_functions_**         contains the missing m-fucntions
%    outputfolder/missing_functions_calls_**   contains the call of missing m-fucntions 
%
%  
%  EXAMPLE
%     
%  tom_compare_fuction_lists('TOM_Release_2008','/fs/pool/pool-bmsan-apps/utils/av3/','tom_','out2');
%
%  
%  REFERENCES
%  
%  SEE ALSO
%     
%  
%     created by FB 12/04/04
%     updated by FB 07/10/05
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


warning off;
mkdir(output);
mkdir([output '/steps']);
warning on;


outNameNew=strrep(new_folder,'/','_');
if (strcmp(outNameNew(end),'_'))
    outNameNew=outNameNew(1:end-1);
end;
if (strcmp(outNameNew(1),'_'))
    outNameNew=outNameNew(2:end);
end;
outNameOld=strrep(old_folder,'/','_');
if (strcmp(outNameOld(end),'_'))
    outNameOld=outNameOld(1:end-1);
end;
if (strcmp(outNameOld(1),'_'))
    outNameOld=outNameOld(2:end);
end;

call=['ls -R ' old_folder ' | grep ' wildcard ' | sort > ' output '/list_' outNameOld '.txt'];
unix(call);

call=['ls -R ' new_folder ' | grep ' wildcard ' | sort > ' output '/list_' outNameNew '.txt'];
unix(call);

call=['find ' new_folder ' -name "*.m" -exec grep "' wildcard '" {} \; > ' output '/steps/list_with_comment.txt'];
unix(call);

call=['cat ' output '/steps/list_with_comment.txt | grep -v "^%" | grep "' wildcard '" > ' output '/steps/list_without_comment.txt' ];
unix(call);

ll=importdata([output '/steps/list_without_comment.txt']);

wt=[strrep(wildcard,'*','') ''];

zz=1;
for i=1:length(ll)
    lin=ll{i};
    if (strfind(lin,'function')==1)
        continue;
    end;
    lin_clean=lin;
    p_com=strfind(lin,'%');
    if (isempty(p_com)==0)
        lin_clean=lin(1:p_com);
    end;
    p_wt=strfind(lin_clean,wt);
    lin_clean=lin(p_wt:end);
    p_br=strfind(lin_clean,'(');
    lin_clean=lin_clean(1:p_br-1);
    if (isempty(lin_clean)==0)
        all_lins{zz}=lin_clean;
        zz=zz+1;
    end;

end;

all_lins=unique(all_lins);



fid=fopen([output '/calls_' outNameNew '.txt'],'wt');
for i=1:length(all_lins)
    fprintf(fid,'%s.m\n',all_lins{i});
end;
fclose(fid);

call=['cat ' output '/list_' outNameOld '.txt | grep -f ' output '/calls_' outNameNew '.txt > ' output '/steps/in_both.txt'];
unix(call);

call=['cat ' output '/calls_' outNameNew '.txt | grep -vf ' output '/steps/in_both.txt > ' output '/steps/missing_all.txt'];
unix(call);

call=['cat ' output '/steps/missing_all.txt | grep -vf ' output '/list_' outNameNew '.txt  > ' output '/steps/missing_all_without_internal.txt'];
unix(call);


missing=importdata([output '/steps/missing_all_without_internal.txt']);

zz=1;
disp('checking for interanl function ');
for i=1:length(missing)
    mm=strtok(missing{i},'.');
    call=['find ' new_folder '/ -name "*.m" -exec grep -H ' mm ' {} \; | grep :function' ];
    [a b]=unix(call);
    if (isempty(b))
        clean_missing{zz}=mm;
        zz=zz+1;
    else
        disp('found internal function:');
        disp(b);
        disp(' ');
        
    end;
end;
disp('done !');


disp([num2str(zz) ' functions remain unresolved']);

%find missing in files
fid=fopen([output '/missing_functions_calls_' outNameNew '.txt'],'wt');
fid2=fopen([output '/missing_functions_' outNameNew '.txt'],'wt');

for i=1:length(clean_missing)
    mm=clean_missing{i};
    disp('*******************************************************************************');
    disp(['************** ' mm '*******************************************************']);
    call=['find ' new_folder '/ -name "*.m" -exec grep -H ' mm ' {} \; | grep -v ":%"' ];
    [a b]=unix(call);
    fprintf(fid,'%s\n',b);
    fprintf(fid2,'%s\n',clean_missing{i});
    disp(b);
end;
fclose(fid);    
fclose(fid2);

disp(' ');
disp('Missing m-files: ');

unix(['cat ' output '/missing_functions_' outNameNew '.txt']);






