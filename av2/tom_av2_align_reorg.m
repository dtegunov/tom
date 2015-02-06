function tom_av2_align_reorg(f_align,new_f_align)

tmp=load(f_align);
align2d=tmp.align2d;
clear('tmp');

for i=1:size(align2d,2) 
    tmp{i}=align2d(1,i).filename; 
    [a b c]=fileparts(tmp{i});
    tmp_path{i}=a;
end;

datasets=unique(tmp_path);
new_al=align2d;
zz=1;

h = waitbar(0,'Please wait...');

for i=1:length(datasets)
   idx_ds=ismember(tmp_path,datasets{i});
   tmp_names=unique({tmp{idx_ds}});
   tmp_names=tom_sort_nat(tmp_names);
   for ii=1:length(tmp_names)
       waitbar(ii./length(tmp_names),h,['Dataset: ' num2str(i) ' of ' num2str(length(datasets)) ]);
       idx_al=find(ismember(tmp,tmp_names{ii}));
       for iii=1:length(idx_al)
            new_al(1,zz)=align2d(1,idx_al(iii));
            zz=zz+1;
       end;
    end; 
end;

align2d=new_al;

clear('new_al');
save(new_f_align,'align2d');

disp('datasets:');
for i=1:length(datasets)
    disp(datasets{i});
end;

try
close(h);
catch
end;

