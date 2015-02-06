function h = tom_os3_generateHistogram(options)


resultDirectory = options.job.resultDirectory;

if(strcmp(resultDirectory(end),'/'))
    slash = '';
else
    slash = '/';
end;

files = dir(resultDirectory);

files = {files.name};

h = [];

bins = -1:0.1:1;

for i=3:length(files)
    
    load([resultDirectory slash files{i}]);
    
    if(isempty(h))
        h = hist(result.ccc,bins) ;
    else
        h = h + hist(result.ccc,bins) ;
    end;
    
    clear result;
end;

h = h / (length(files)-2);