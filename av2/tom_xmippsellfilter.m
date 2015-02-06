function tom_xmippsellfilter(doc_filename_in,sell_filename_out_include,sell_filename_out_exclude,thres)

st=tom_xmippdocread(doc_filename_in);

for i=1:size(st,2); 
    tmp(i)=st(i).pmax_sump; 
end;

idx_exclude=find(tmp <= thres);
idx_include=find(tmp > thres);


fid=fopen(sell_filename_out_exclude,'w');
for i=1:size(idx_exclude,2)
    fprintf(fid,'%s 1 \n',st(idx_exclude(i)).name);
end;
fclose(fid);


fid=fopen(sell_filename_out_include,'w');
for i=1:size(idx_include,2)
    fprintf(fid,'%s 1 \n',st(idx_include(i)).name );
end;
fclose(fid);






