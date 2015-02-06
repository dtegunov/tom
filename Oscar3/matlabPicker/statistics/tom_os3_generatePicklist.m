function picklist = tom_os3_generatePicklist(options)


resultDirectory = options.job.resultDirectory;

if(strcmp(resultDirectory(end),'/'))
    slash = '';
else
    slash = '/';
end;

files = dir(resultDirectory);

files = {files.name};

picklist = {};

for i=3:length(files)
    
    load([resultDirectory slash files{i}]);
    
    picklist{length(picklist)+1} = tom_os3_returnPicklist(tom_os3_peakValue(result.ccc,result.psr,result.autoc),result.angles,result.job.templateSize,50,options);
    
    clear result;
end;
