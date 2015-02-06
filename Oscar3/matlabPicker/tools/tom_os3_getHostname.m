function hostname = tom_os3_getHostname


[a hostname]=system('hostname');

if(strcmp(hostname,''))
    [a hostname]=system('hostname -i');
end;

hostname = hostname(1:end-1);

