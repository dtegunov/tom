function [result jobs]= tom_os3_checkMissingResults(folder,prefix,numberFiles,jobs)


files = dir([folder '*mat']);


files = {files.name};
preLength = length(prefix);

result = zeros(1,numberFiles);

for i = 1:length(files)
    
    name = files{i};
    
    name = name(preLength+1:end);
    
    pos = strfind(name,'-');
    
    index = str2double(name(1:pos-1));
    
    result(index)=1;
end;

j = {};

if(exist('jobs'))
    for i=1:numberFiles
        if(result(i) == 0)
            j{length(j)+1} = jobs{i};
        end;
    end;
    jobs = j;
else
    jobs ={};
end;
