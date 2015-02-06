function [st,align2d]=tom_av2_check_written(outfilestruct,align2d,stop)


if nargin < 2
    align2d='';
end;

if nargin < 3
    stop='no_stop';
end;

%basepath
base_path=outfilestruct.path(1:(max(strfind(outfilestruct.path,'/'))));

%get folders
st_fold=dir([outfilestruct.path '*']);
zz=1;

vect=zeros(1,10000000);


fprintf('\n %s','checking folder: ');
for i=1:size(st_fold,1)
    %get files
    if (st_fold(i).isdir==1)
        st_files=dir([base_path st_fold(i).name '/' outfilestruct.filename '*']);
        for ii=1:size(st_files,1);
            num=tom_get_max_numeric({st_files(ii).name});
            if isempty(num)
                continue;
            end;
            vect(zz)=num;
            zz=zz+1;
            if (nargin > 1)
                align2d(1,i).iswritten=1;
            end;
        end;

        if (strcmp(stop,'no_stop')==0)
            if (zz > stop)
                st.written_index=sort(vect);
                st.max_abs_num=max(vect);
                st.num_of_files=zz;
                return;
            end;
        end;
        if mod(i,10)~=0
            fprintf('%d ',tom_get_max_numeric({st_fold(i).name}));
        else
            fprintf('%d \n',tom_get_max_numeric({st_fold(i).name}));
        end;
    end;
    
end;
fprintf('%s \n',' done!');

vect=vect(1:zz-1);

st.written_index=sort(vect);
st.max_abs_num=max(vect);
st.num_of_files=zz;

