function st_out=tom_xmippdoc2classstruct(filename)

st=tom_xmippdocread(filename);


%determine number of classes !
for i=1:size(st,2)
    tmp_ref(i)=st(i).ref;
    disp(num2str(i));
end;

num_of_classes=max(tmp_ref);



for i=1:size(st,2)
    ref_nr=st(i).ref;
    st_out(ref_nr)=st(i);
end;

