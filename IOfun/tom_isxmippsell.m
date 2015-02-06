function out=tom_isxmippsell(filename)

flag=0;

[path name ext]=fileparts(filename);

if (strcmp(ext,'.sel')==0)
    out=0;
    return;
end;




fid=fopen(filename);
line=fgetl(fid);
fclose(fid);

line=strtok(line,' ');

if (fid==-1)
    out=0;
    return;
end;

try 
    im=tom_spiderread(line);
    flag=1;
catch
    disp([line ' not found. Change directory.']);
end;
 

try
    im=tom_spiderread(['.' line]);
    flag=2;    
catch

end;

if (flag==0)
    out=0;
    return;
end;


out=1;


