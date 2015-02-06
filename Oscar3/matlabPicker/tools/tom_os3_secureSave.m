function [written resultName resultDirectory]= tom_os3_secureSave(result,job,options,logFile,doOnline)




%%  save the results to file in specified path
%volumeName
posOfSlash = strfind(job.volumeFile,'/');
volumeName = job.volumeFile(posOfSlash(end)+1:end);
posOfSlash = strfind(job.templateFile,'/');
if(~ isempty(posOfSlash))
    templateName = job.templateFile(posOfSlash(end)+1:end);
else
    templateName = job.templateFile;
end;

resultDirectory = options.job.resultDirectory;
 
if(strcmp(options.job.jobType,'2d-online') && doOnline)    
    if(~exist([options.job.resultDirectory '/' volumeName],'dir'))
        mkdir([options.job.resultDirectory '/' volumeName]);
        system(['chmod -R 777 ' options.job.resultDirectory '/' volumeName ]);
    end;    
    
    resultDirectory = [ resultDirectory '/' volumeName ];
end;

resultName = [ resultDirectory '/' volumeName '-result-' num2str(job.id) '-' templateName '.mat' ];
save(resultName,'result');

result2= result;
clear result;

%% open the logFile
logOpendHere = false;

if(~exist('logFile') || logFile == -1)

    hostname = tom_os3_getHostname;
    logFile = fopen([ options.job.resultDirectory '/logFile-' hostname '.txt'],'a');
    logOpendHere = true;

end;


%%  clean up

try
    %try to load the file
    load(resultName);
    written = true;
catch
    %if the mat file could not be opened, retry 10 times, else skip and report to log file
    i =10;
    written =false;
    system([ 'rm ' resultName ]);
    result = result2;
    save(resultName,'result');
    clear result;

    fprintf(logFile,[ 'Problems when reading result of job no: ' num2str(job.id) ' Retry ' num2str(i) '\n'  ]);
    fclose(logFile);

    while(i>0 && ~written)

        try
            load([ resultName ]);
            written =true;
        catch
            logFile = fopen([ options.job.resultDirectory '/logFile.txt'],'a');
            fprintf(logFile,[ 'Problems when reading result of job no: ' num2str(job.id) ' Retry ' num2str(i) '\n' ]);
            fclose(logFile);
            system([ 'rm ' resultName ]);
            result = result2;
            save(resultName,'result');
            clear result;
            i = i-1;
        end;
    end;

end;

if(written)
    fprintf(logFile,[ 'Result of job no: ' num2str(job.id) ' written to disk \n' ]);
else
    fprintf(logFile,[ 'Could not ensure the result of job no: ' num2str(job.id) ' written to disk is read able. \n' ]);
end;


if(logOpendHere)
    fclose(logFile);
end;
