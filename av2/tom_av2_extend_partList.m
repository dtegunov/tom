function counts=tom_av2_extend_partList(list)

dat=importdata(list);

counts=zeros(length(dat.textdata),1);
tic;
for i=1:length(dat.textdata)-1
    try
        v = tom_emreadc3(dat.textdata{i+1,1}, [dat.data(i,2)-128,dat.data(i,3)-128,0 256,256,1], [8,8,1], [1,1,1]);
    catch ME
        disp([ME.message]);
        v = tom_emreadc3(dat.textdata{i+1,1}, [dat.data(i,2)-32,dat.data(i,3)-32,0 64,64,1], [2,2,1], [1,1,1]);
    end;
    
    counts(i)=mean2(v.Value);
    if (mod(i,500)==0)
        toc;
        disp(num2str(i));
        tic;
    end;
end;
toc;

tom_dev(counts)

fp=fopen([list '.ext'],'wt');

for idx=1:length(dat.textdata)-1
    fprintf(fp,'%s %d %s %d %d %d %f %f %f %f %f  \n',dat.textdata{idx+1,1},str2double(dat.textdata{idx+1,2}),dat.textdata{idx+1,3},dat.data(idx,1),dat.data(idx,2),dat.data(idx,3),dat.data(idx,4),dat.data(idx,5),dat.data(idx,6),dat.data(idx,7),counts(idx));
    if (mod(idx,1000)==0)
        disp(num2str(idx));
    end;
end;

fclose(fp);