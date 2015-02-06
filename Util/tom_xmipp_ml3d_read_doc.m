function [head,s,f]=tom_xmipp_ml3d_read_doc(filename)

% [head,s,f]=tom_xmipp_ml3d_read_doc(filename)
% reads xmipp ML3D doc file and gives back the header
% the filenames and the alignment data, see header
%
[fp,msg]=fopen(filename,'r');
if fp==-1
    error(['Open File: ' msg]); return;
end;
head = fgetl(fp);
i=1;
while feof(fp)~=1
    s = fgetl(fp);
    i=i+1;
end
istep=round(i./20);
fseek(fp,0,-1);
head = fgetl(fp);
f=zeros(14,i./2,'single');
disp([filename ', total of: ' num2str((i-1)./2) ' particles.']);
i=1;
clear s;
while feof(fp)~=1
    tmp=fgetl(fp);
    s{:,i} = tmp(4:end);
    f(:,i) = fscanf(fp,'%f',14);
    c = fscanf(fp,'%c',2);
    i=i+1;
    if i./istep==fix(i./istep)
        p=sprintf('%2.2f %%', 10.*i./istep);     
        disp(p);
    end;
end
fclose(fp);
disp('Statistics of LogLikelihood:');
tom_dev(f(10,:));

