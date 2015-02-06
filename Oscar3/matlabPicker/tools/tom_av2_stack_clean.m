function tom_av2_stack_clean(outfilestruct,flag)

if (nargin < 2)
    flag='no_outfilestruct';
end;

%get folders
st_fold=dir([outfilestruct.path '*']);

baspath=fileparts(outfilestruct.path);
warning off;

for i=1:size(st_fold,1)
    if  (st_fold(i).isdir==1)
        unix(['rm -f ' baspath '/' st_fold(i).name '/*.*']);
        unix(['rm -Rf ' baspath '/' st_fold(i).name]);
    end;
end;


if (strcmp(flag,'outfilestruct'))
    try
        unix(['rm -f ' outfilestruct.path '.mat']);
    catch
    end;
end;

warning on;